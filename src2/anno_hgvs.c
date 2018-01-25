#include "utils.h"
#include "hgvs.h"
#include "anno_hgvs.h"
#include "anno_col.h"
#include "variant_type.h"
#include "htslib/vcf.h"
#include "stack_lite.h"
#include "name_list.h"

static int anno_hgvs_update_buffer(struct anno_hgvs_file *f, bcf_hdr_t *hdr, bcf1_t *line)
{
    int i;
    // clear buffer
    for ( i = 0; i < f->n_allele; ++i ) hgvs_destroy(f->files[i]);
    if ( f->n_allele ) free(f->files);
    f->n_allele = line->n_allele-1; // emit ref
    f->files = malloc(f->n_allele*sizeof(struct hgvs));
    for ( i = 0; i < f->n_allele; ++i ) 
        f->files[i] = hgvs_init(bcf_seqname(hdr, line), line->pos+1, line->pos+line->rlen, line->d.allele[0], line->d.allele[i+1]);
    return 0;
}
// generate variant string in annovar format
// OR4F5:NM_001005484:exon1:c.T809A:p.V270E
static char *generate_annovar_name(struct hgvs *h)
{
    kstring_t str = {0,0,0};
    int i;
    for ( i = 0; i < h->n_tran; ++i ) {
        if ( i ) kputc('|', &str);
        struct hgvs_type *type = &h->trans[i].type;
        struct hgvs_inf  *inf  = &h->trans[i].inf;
        // for noncoding donot export annovar name
        if ( type->func1 != func_region_cds || inf->offset != 0 ) {
            kputc('.', &str);
            continue;
        }
                
        char *ref, *alt;
        if ( inf->strand == '+' ) {
            ref = h->ref ? strdup(h->ref) : NULL;
            alt = h->alt ? strdup(h->alt) : NULL;
        }
        else {
            ref = h->ref ? rev_seqs(h->ref, strlen(h->ref)) : NULL;
            alt = h->alt ? rev_seqs(h->alt, strlen(h->alt)) : NULL;
        }


        
        ksprintf(&str,"%s:%s:", inf->gene, inf->transcript);
        if ( inf->offset != 0 ) ksprintf(&str, "intron%d:", type->count);
        else ksprintf(&str, "exon%d:", type->count);
        if ( type->func1 == func_region_noncoding ) kputs("n.", &str);
        else if ( type->func1 == func_region_cds  ) kputs("c.", &str);
        else if ( type->func1 == func_region_utr5 ) kputs("c.-", &str);
        else if ( type->func1 == func_region_utr3 ) kputs("c.*", &str);

        if ( inf->dup_offset > 0 ) {
            assert(inf->loc != 0 );
            int len = strlen(h->alt);
            if ( len > 1 ) ksprintf(&str, "%d_%ddup%s", inf->loc+inf->dup_offset-len+1, inf->loc+inf->dup_offset+1,h->alt);
            else ksprintf(&str, "%ddup%s", inf->loc+inf->dup_offset, h->alt);
        }
        else if ( inf->dup_offset < 0 ) {
            assert(inf->loc != 0 );
            int len = strlen(h->alt);
            if ( len > 1 ) ksprintf(&str, "%d_%ddup%s", inf->loc+inf->dup_offset+1, inf->loc, h->alt);
            else {
                assert(inf->dup_offset == -1);
                ksprintf(&str, "%ddup%s", inf->loc, h->alt);
            }
        }
        else {
            if ( h->type == var_type_snp ) {
                if ( ref ) ksprintf(&str, "%s", ref);            
                if ( inf->loc != 0 ) ksprintf(&str,"%d", inf->loc);
                else kputc('?', &str);
                
                if ( inf->offset > 0 ) ksprintf(&str, "+%d", inf->offset);
                else if ( inf->offset < 0 ) ksprintf(&str, "%d", inf->offset);
                ksprintf(&str, "%s", alt);
            }
            else {
                if ( inf->loc != 0 ) ksprintf(&str,"%d", inf->loc);
                else kputc('?', &str);
                
                if ( inf->offset > 0 ) ksprintf(&str, "+%d", inf->offset);
                else if ( inf->offset < 0 ) ksprintf(&str, "%d", inf->offset);

                if ( h->start != h->end ) {
                    kputc('_', &str);
                    if ( type->func2 == func_region_utr5 ) kputc('-', &str);
                    else if ( type->func2 == func_region_utr3 ) kputc('*', &str);
                    ksprintf(&str, "%d", inf->end_loc);
                    if ( inf->end_offset > 0 ) ksprintf(&str,"+%d", inf->end_offset);
                    else if ( inf->end_offset < 0 ) ksprintf(&str,"%d", inf->end_offset);
                }
                if ( h->type == var_type_del ) kputs("del", &str);
                else if ( h->type == var_type_ins) ksprintf(&str, "ins%s", alt);
                else if ( h->type == var_type_delins) ksprintf(&str, "delins%s", alt);
                else error("Variant type inconsistant. %d", h->type);           
            }
        }

        if ( ref ) free(ref);
        if ( alt ) free(alt);

        if ( type->loc_amino > 0 && h->type == var_type_snp) ksprintf(&str, ":p.%s%d%s", codon_short_names[type->ori_amino], type->loc_amino, codon_short_names[type->mut_amino]);
        else {
            int i;
            if ( type->vartype == var_is_inframe_insertion ) {
                assert(type->loc_end_amino -1 == type->loc_amino);
                ksprintf(&str, ":p.%d_%dins", type->loc_amino, type->loc_end_amino);
                for (i = 0; i < type->n; ++i) kputs(codon_short_names[type->aminos[i]], &str);
            }
            else if ( type->vartype == var_is_inframe_deletion ) {
                assert(type->loc_end_amino > 0);
                if ( type->loc_end_amino == type->loc_amino ) ksprintf(&str, ":p.%ddel%s", type->loc_amino, ref);
                else ksprintf(&str, ":p.%d_%ddel", type->loc_amino, type->loc_end_amino);
                //for (i = 0; i < type->n; ++i) kputs(codon_short_names[type->aminos[i]], &str);
            }
            else if ( type->vartype == var_is_inframe_delins ) {
                assert(type->loc_end_amino > 0);
                if ( type->loc_end_amino == type->loc_amino ) ksprintf(&str, ":p.%ddelins", type->loc_amino);
                else ksprintf(&str, ":p.%d_%ddelins", type->loc_amino, type->loc_end_amino); 
                for (i = 0; i < type->n; ++i) kputs(codon_short_names[type->aminos[i]], &str);
            }
            else if ( type->vartype == var_is_frameshift ) {
                ksprintf(&str, ":p.%s%d",codon_short_names[type->ori_amino], type->loc_amino);
                if ( type->fs > 0 ) kputs("fs", &str);
            }
        }
    }
    return str.s;
}
static char *generate_hgvsnom_string(struct hgvs *h)
{
    kstring_t str = {0,0,0};
    int i;
    for ( i = 0; i < h->n_tran; ++i ) {
        if ( i ) kputc('|', &str);
        struct hgvs_type *type = &h->trans[i].type;
        struct hgvs_inf *inf = &h->trans[i].inf;
        ksprintf(&str, "%s", inf->transcript);
        int l = strlen (inf->transcript);
        int k;
        for ( k = 0; k < l; ++k )
            if ( inf->transcript[k] == '.') break;
        if ( k == l && inf->version > 0 ) ksprintf(&str, ".%d", inf->version);
        kputc(':', &str);
        if ( type->func1 == func_region_noncoding ) kputs("n.", &str);
        else if ( type->func1 == func_region_cds ) kputs("c.", &str);
        else if ( type->func1 == func_region_utr5 ) kputs("c.-", &str);
        else if ( type->func1 == func_region_utr3 ) kputs("c.*", &str);

        // for dup 
        if ( inf->dup_offset > 0 ) {           
            assert(inf->loc != 0);
            int len = strlen(h->alt);
            if ( len > 1 ) ksprintf(&str, "%d_%ddup", inf->loc+inf->dup_offset-len+1, inf->loc+inf->dup_offset+1);
            else ksprintf(&str, "%ddup", inf->loc+inf->dup_offset);
        }
        else if ( inf->dup_offset < 0 ) {
            assert(inf->loc != 0 );
            int len = strlen(h->alt);
            if ( len > 1 ) ksprintf(&str, "%d_%ddup", inf->loc+inf->dup_offset+1, inf->loc);
            else {
                assert(inf->dup_offset == -1);
                ksprintf(&str, "%ddup", inf->loc);
            }
        }
        else {
            // for nondup, inf->dup_offset will always be 0
            if ( inf->loc != 0 ) ksprintf(&str, "%d", inf->loc);
            else kputc('?', &str);
            
            if ( inf->offset > 0 ) ksprintf(&str, "+%d", inf->offset);
            else if ( inf->offset < 0 ) ksprintf(&str, "%d", inf->offset);

            while ( h->start != h->end) {
                //if ( inf->dup_offset && h->end - h->start == 1 ) break;
                kputc('_', &str);
                if ( type->func2 == func_region_utr5 ) kputc('-', &str);
                else if ( type->func2 == func_region_utr3 ) kputc('*', &str);
                ksprintf(&str, "%d", inf->end_loc);
                if ( inf->end_offset > 0 ) ksprintf(&str,"+%d", inf->end_offset);
                else if ( inf->end_offset < 0 ) ksprintf(&str,"%d", inf->end_offset);
                break;
            }
            char *ref, *alt;
            if ( inf->strand == '+' ) {
                ref = h->ref ? strdup(h->ref) : NULL;
                alt = h->alt ? strdup(h->alt) : NULL;
            }
            else {
                ref = h->ref ? rev_seqs(h->ref, strlen(h->ref)) : NULL;
                alt = h->alt ? rev_seqs(h->alt, strlen(h->alt)) : NULL;
            }
            if ( h->type == var_type_snp ) ksprintf(&str, "%s>%s", ref, alt);
            else if ( h->type == var_type_del ) ksprintf(&str, "del%s", ref);
            else if ( h->type == var_type_ins ) ksprintf(&str, "ins%s", alt);
            else if ( h->type == var_type_delins ) ksprintf(&str, "delins%s", alt);
            else error("Failed to parse HGVS nom.");

            if ( ref ) free(ref);
            if ( alt ) free(alt);
        }

        // protein code
        if ( type->loc_amino > 0 && h->type == var_type_snp ) 
            ksprintf(&str, "(p.%s%d%s/p.%s%d%s)", codon_names[type->ori_amino], type->loc_amino, codon_names[type->mut_amino], codon_short_names[type->ori_amino], type->loc_amino, codon_short_names[type->mut_amino]);
        // indels
        else {
            int i;
            if ( type->vartype == var_is_inframe_insertion ) {
                assert(type->loc_amino +1 == type->loc_end_amino);
                ksprintf(&str, "(p.%s%d_%s%dins",codon_names[type->ori_amino], type->loc_amino, codon_names[type->ori_end_amino], type->loc_end_amino);
                for (i = 0; i < type->n; ++i) kputs(codon_names[type->aminos[i]], &str);
                kputc(')', &str);
            }
            else if ( type->vartype == var_is_inframe_deletion ) {
                assert(type->loc_end_amino != 0);
                if ( type->loc_end_amino == type->loc_amino ) ksprintf(&str, "(p.%s%ddel",codon_names[type->ori_amino], type->loc_amino);
                else  ksprintf(&str, "(p.%s%d_%s%ddel",codon_names[type->ori_amino], type->loc_amino, codon_names[type->ori_end_amino], type->loc_end_amino);
                for (i = 0; i < type->n; ++i) kputs(codon_names[type->aminos[i]], &str);
                kputc(')', &str);                                         
            }
            else if ( type->vartype == var_is_inframe_delins ) {
                if ( type->loc_end_amino == type->loc_amino ) ksprintf(&str, "(p.%s%ddelins",codon_names[type->ori_amino], type->loc_amino);
                else ksprintf(&str, "(p.%s%d_%s%ddelins",codon_names[type->ori_amino], type->loc_amino, codon_names[type->ori_end_amino], type->loc_end_amino); 
                for (i = 0; i < type->n; ++i) kputs(codon_names[type->aminos[i]], &str);
                kputc(')', &str);                                         
            }
            else if ( type->vartype == var_is_frameshift ) {
                ksprintf(&str, "(p.%s%d%s",codon_names[type->ori_amino], type->loc_amino, codon_names[type->mut_amino]);
                if ( type->fs > 0 ) ksprintf(&str, "fs*%d", type->fs);
                kputc(')', &str);                                         
            }
        }
    }
    
    return str.s;   
}
static char *generate_gene_string(struct hgvs *h)
{
    kstring_t str = {0,0,0};
    //struct anno_stack *s = anno_stack_init();
    int i;
    for ( i = 0; i < h->n_tran; ++i ) {
     if ( i ) kputc('|', &str);
        struct hgvs_inf *inf = &h->trans[i].inf;
        if ( inf->gene ) kputs(inf->gene, &str);
        else kputc('.', &str);
    }   
    return str.s;
}
static char *generate_transcript_string(struct hgvs *h)
{
    kstring_t str = {0,0,0};
    int i, k, l;
    for ( i = 0; i < h->n_tran; ++i ) {
        if ( i ) kputc('|', &str);
        struct hgvs_inf *inf = &h->trans[i].inf;
        kputs(inf->transcript, &str);
        l = strlen(inf->transcript);
        for ( k = 0; k < l; ++k )
            if ( inf->transcript[k] == '.') break;
        if ( k == l && inf->version > 0 ) ksprintf(&str, ".%d", inf->version);
    }
    return str.s;
}
static char *generate_vartype_string(struct hgvs *h)
{
    kstring_t str = {0,0,0};
    int i;
    for ( i = 0; i < h->n_tran; ++i ) {
        if ( i ) kputc('|', &str);
        struct hgvs_type *type = &h->trans[i].type;
        if ( type->vartype2 != var_is_not_splice ) {
            kputs(var_type_splice_string(type->vartype2), &str);
            if (type->vartype != var_is_intron && type->vartype != var_is_reference)
                ksprintf(&str, "(%s)", var_func_type_string(type->vartype));
        }
        else kputs(var_func_type_string(type->vartype), &str);
    }
    return str.s;    
}
static char *generate_exonintron_string(struct hgvs *h)
{
    kstring_t str = {0,0,0};
    int i;
    for ( i = 0; i < h->n_tran; ++i ) {
        if ( i ) kputc('|', &str);
        struct hgvs_type *type = &h->trans[i].type;
        struct hgvs_inf *inf = &h->trans[i].inf;
        if ( inf->offset != 0 )
            ksprintf(&str, "I%d", type->count);
        else {
            ksprintf(&str, "E%d", type->count);
            if ( type->count2 != 0 )
                ksprintf(&str, "/C%d", type->count2);
        }
    }
    return str.s;
}
static char *generate_ivsnom_string(struct hgvs *h)
{
    kstring_t str = {0,0,0};
    int i;
    int empty = 1;
    for ( i = 0; i < h->n_tran; ++i ) {        
        struct hgvs_type *type = &h->trans[i].type;
        struct hgvs_inf *inf = &h->trans[i].inf;
        if ( i ) kputc('|', &str);
        if ( inf->offset == 0 ) { kputc('.', &str); continue; }
        ksprintf(&str, "c.IVS%d", type->count);
        if ( inf->offset > 0 ) ksprintf(&str, "+%d", inf->offset);
        else if ( inf->offset < 0 ) ksprintf(&str, "%d", inf->offset);
        if ( h->start != h->end ) {
            kputc('_', &str);
            if (inf->offset == 0 || inf->end_offset == 0 ) kputc('?', &str);
            else if ( inf->end_offset > 0 ) ksprintf(&str, "+%d", inf->end_offset);
            else if ( inf->end_offset < 0 ) ksprintf(&str, "%d", inf->end_offset);
        }
        if ( h->type == var_type_del ) {
            assert(h->ref);
            if ( inf->strand == '+' ) ksprintf(&str, "del%s", h->ref);
            else {
                char *rev = rev_seqs(h->ref, strlen(h->ref));
                ksprintf(&str, "del%s", rev);
                free(rev);
            }
        }
        else if ( h->type == var_type_ins ) {
            assert(h->alt);
            if ( inf->strand == '+' ) ksprintf(&str, "ins%s", h->alt);
            else {
                char *rev = rev_seqs(h->alt, strlen(h->alt));
                ksprintf(&str, "ins%s", rev);
                free(rev);
            }
        }
        else {
            assert(h->ref && h->alt);
            if ( inf->strand == '+' ) ksprintf(&str, "%s>%s", h->ref, h->alt);
            else {
                char *ref = h->ref == NULL ? NULL : rev_seqs(h->ref, strlen(h->ref));
                char *alt = h->alt == NULL ? NULL : rev_seqs(h->alt, strlen(h->alt));
                ksprintf(&str, "%s>%s", ref, alt);
                free(ref); free(alt);
            }
        }
        empty = 0;
    }
    if ( empty == 1 ) {
        free(str.s);
        return NULL;
    }
    return str.s;
}
static char *generate_oldnom_string(struct hgvs *h)    
{
    kstring_t str = {0,0,0};
    int i;
    for ( i = 0; i < h->n_tran; ++i ) {
        if ( i ) kputc('|', &str);
        //struct hgvs_type *type = &h->trans[i].type;
        struct hgvs_inf *inf = &h->trans[i].inf;
        ksprintf(&str, "%s", inf->transcript);
        int l = strlen (inf->transcript);
        int k;
        for ( k = 0; k < l; ++k )
            if ( inf->transcript[k] == '.') break;
        if ( k == l && inf->version > 0 ) ksprintf(&str, ".%d", inf->version);
        kputc(':', &str);
        if ( inf->pos != 0 ) ksprintf(&str, "n.%d", inf->pos);
        else kputs("n.?", &str);

        if ( inf->offset > 0 ) ksprintf(&str, "+%d", inf->offset);
        else if ( inf->offset < 0 ) ksprintf(&str, "%d", inf->offset);

        if ( h->start != h->end) {
            kputc('_', &str);
            ksprintf(&str, "%d", inf->end_pos);
            if ( inf->end_offset > 0 ) ksprintf(&str,"+%d", inf->end_offset);
            else if ( inf->end_offset < 0 ) ksprintf(&str,"%d", inf->end_offset);
        }
        char *ref, *alt;
        if ( inf->strand == '+' ) {
            ref = h->ref ? strdup(h->ref) : NULL;
            alt = h->alt ? strdup(h->alt) : NULL;
        }
        else {
            ref = h->ref ? rev_seqs(h->ref, strlen(h->ref)) : NULL;
            alt = h->alt ? rev_seqs(h->alt, strlen(h->alt)) : NULL;
        }
        if ( h->type == var_type_snp ) ksprintf(&str, "%s>%s", ref, alt);
        else if ( h->type == var_type_del ) ksprintf(&str, "del%s", ref);
        else if ( h->type == var_type_ins ) ksprintf(&str, "ins%s", alt);        
        else if ( h->type == var_type_delins ) ksprintf(&str, "%s>%s", ref, alt);
        else {
            error("Failed to parse HGVS nom.");
        }
        if ( ref ) free(ref);
        if ( alt ) free(alt);
    }
    return str.s;

}
static char *generate_aalength_string(struct hgvs *h)
{
    kstring_t str = {0,0,0};
    int i;
    for ( i = 0; i < h->n_tran; ++i ) {
        if ( i ) kputc('|', &str);
        kputw(h->trans[i].inf.aa_length, &str);
    }
    return str.s;
}
//static int anno_hgvs_setter_info(struct anno_hgvs_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col)
static int anno_hgvs_setter_hgvsnom(struct anno_hgvs_file *file, bcf_hdr_t *hdr, bcf1_t *line)
{
    int i, j;
    struct hgvs_handler *h = file->h;    
    int empty = 1;
    kstring_t *str = malloc(file->n_col*sizeof(kstring_t));
    for ( i = 0; i < file->n_col; ++i ) memset(&str[i], 0, sizeof(kstring_t));
    
    for ( i = 0; i < file->n_allele; ++i ) {
        struct hgvs *f = file->files[i];

        if ( i > 0 ) 
            for ( j = 0; j < file->n_col; ++j ) kputc(',', &str[j]);
        
        if ( hgvs_anno_trans(f, h) == 0 ) {
            for ( j = 0; j < file->n_col; ++j ) kputc('.', &str[j]);                
            continue;
        }
        for ( j = 0; j < file->n_col; ++j ) {
            struct anno_col *col = &file->cols[j];
            if ( strcmp(col->hdr_key, "HGVSnom") == 0 ) {
                char *name = generate_hgvsnom_string(f);
                kputs(name, &str[j]);
                if (name) free(name);
            }
            else if ( strcmp(col->hdr_key, "Gene") == 0 ) {
                char *name = generate_gene_string(f);
                kputs(name, &str[j]);
                if (name) free(name);
            }
            else if ( strcmp(col->hdr_key, "Transcript") == 0 ) {
                char *name = generate_transcript_string(f);
                kputs(name, &str[j]);
                if (name) free(name);                
            }
            else if ( strcmp(col->hdr_key, "VarType") == 0 ) {
                char *name = generate_vartype_string(f);
                kputs(name, &str[j]);
                if (name) free(name);
            }
            else if ( strcmp(col->hdr_key, "ExonIntron") == 0 ) {
                char *name = generate_exonintron_string(f);
                kputs(name, &str[j]);
                if (name) free(name);
            }
            else if ( strcmp(col->hdr_key, "IVSnom") == 0 ) {
                char *name = generate_ivsnom_string(f);
                if ( name ) {
                    kputs(name, &str[j]);
                    free(name);
                }
            }
            else if ( strcmp(col->hdr_key, "Oldnom") == 0 ) {
                char *name = generate_oldnom_string(f);
                kputs(name, &str[j]);
                if (name) free(name);
            }            
            else if ( strcmp(col->hdr_key, "AAlength") == 0 ) {
                char *name = generate_aalength_string(f);
                kputs(name, &str[j]);
                if (name) free(name);
            }
            else if ( strcmp(col->hdr_key, "ANNOVARname") == 0 ) {
                char *name = generate_annovar_name(f);
                kputs(name, &str[j]);
                if ( name ) free(name);
            }
                
        }
        empty = 0;
    }
    if ( empty == 1 ) {
        for ( j = 0; j < file->n_col; ++j ) free(str[j].s);
        free(str);
        return 1;
    }
    // update INFO
    for ( i = 0; i < file->n_col; ++i ) {
        struct anno_col *col = &file->cols[i];
        if ( col->replace == REPLACE_MISSING ) {
            int ret = bcf_get_info_string(hdr, line, col->hdr_key, &file->tmps, &file->mtmps);
            if ( ret > 0 && (file->tmps[0]!= '.' || file->tmps[1] != 0 ) ) continue;
        }
        bcf_update_info_string_fixed(hdr, line, col->hdr_key, str[i].s);
    }
    
    for ( i = 0; i < file->n_col; ++i ) free(str[i].s);
    free(str);
    return 0;
}
struct anno_hgvs_file *anno_hgvs_file_duplicate(struct anno_hgvs_file *f)
{
    struct anno_hgvs_file *d = malloc(sizeof(*d));
    d->h = hgvs_handler_duplicate(f->h);
    d->n_allele = 0;
    d->files = NULL;
    d->n_col = f->n_col;
    d->cols = malloc(d->n_col*sizeof(struct anno_col));
    int i;
    for ( i = 0; i < d->n_col; ++i) anno_col_copy(&f->cols[i], &d->cols[i]);
    return d;
}
void anno_hgvs_file_destroy(struct anno_hgvs_file *f)
{
    int i;
    // clear buffer
    for ( i = 0; i < f->n_allele; ++i ) hgvs_destroy(f->files[i]);
    if ( f->n_allele ) free(f->files);
    hgvs_handler_destroy(f->h);
    //if ( f->tmps) free(f->tmps);
    for ( i = 0; i < f->n_col; ++i ) free(f->cols[i].hdr_key);
    free(f);
}
struct anno_hgvs_file *anno_hgvs_file_init(bcf_hdr_t *hdr, const char *column, const char *data, const char *rna, const char *reference)
{
    if ( data == NULL || rna == NULL ) return NULL;
        
    struct anno_hgvs_file *f = malloc(sizeof(*f));
    memset(f, 0, sizeof(*f));
    
    if ( column == NULL ) {
      full_column:
        f->n_col = 8;
        f->cols = malloc(8*sizeof(struct anno_col));
        memset(f->cols, 0, 8*sizeof(struct anno_col));
        f->cols[0].hdr_key = strdup("HGVSnom");
        f->cols[1].hdr_key = strdup("Gene");
        f->cols[2].hdr_key = strdup("Transcript");
        f->cols[3].hdr_key = strdup("VarType");
        f->cols[4].hdr_key = strdup("ExonIntron");
        f->cols[5].hdr_key = strdup("IVSnom");
        f->cols[6].hdr_key = strdup("Oldnom");
        f->cols[7].hdr_key = strdup("AAlength");        
    }
    else {
        kstring_t str = {0,0,0};
        kputs(column, &str);
        int i, n;
        int *s = ksplit(&str, ',', &n);
        if ( n == 0 ) {
            warnings("Failed to parse columns, try to annotate all tags. %s.", column);
            goto full_column;            
        }
        f->cols = malloc(n*sizeof(struct anno_col));
        for ( i = 0; i < n; ++i ) {
            char *ss = str.s + s[i];
            struct anno_col *col = &f->cols[f->n_col];
            memset(col, 0, sizeof(struct anno_col));
            //col->icol = -1;
            col->replace = REPLACE_MISSING;
            if ( *ss == '+' ) ss++;
            else if (*ss == '-') { col->replace = REPLACE_EXISTING; ss++; }
            if ( ss[0] == '\0') continue;
            if ( strncmp(ss, "INFO/", 5) == 0 ) ss += 5;
            if ( strcmp(ss, "HGVSnom")  && strcmp(ss, "Gene") && strcmp(ss, "Transcript") && strcmp(ss, "VarType") && strcmp(ss, "ExonIntron") && strcmp(ss, "IVSnom")
                 && strcmp(ss, "Oldnom") && strcmp(ss, "AAlength") && strcmp(ss, "ANNOVARname") ){
                warnings("Do NOT support tag %s.", ss);
                continue;
            }
            col->hdr_key = strdup(ss);            
            f->n_col++;
        }
        free(s);
    }    

    // update header
#define BRANCH(_key, _description) do {                                 \
        int id;                                                         \
        id = bcf_hdr_id2int(hdr, BCF_DT_ID, _key);                      \
        if (id == -1) {                                                 \
            bcf_hdr_append(hdr, _description);                          \
            bcf_hdr_sync(hdr);                                          \
            id = bcf_hdr_id2int(hdr, BCF_DT_ID, _key);                  \
            assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));        \
        }                                                               \
    } while(0)
    int i;
    for ( i = 0; i < f->n_col; ++i ) {
        struct anno_col *col = &f->cols[i];
        if ( strcmp(col->hdr_key, "HGVSnom") == 0 ) 
            BRANCH("HGVSnom", "##INFO=<ID=HGVSnom,Number=A,Type=String,Description=\"HGVS nomenclature for the description of DNA sequence variants\">");
        else if ( strcmp(col->hdr_key, "ANNOVARname") == 0 )
            BRANCH("ANNOVARname", "##INFO=<ID=ANNOVARname,Number=A,Type=String,Description=\"Variant description in ANNOVAR format.\">");
        else if ( strcmp(col->hdr_key, "Gene") == 0 ) 
            BRANCH("Gene","##INFO=<ID=Gene,Number=A,Type=String,Description=\"Gene names\">");
        else if ( strcmp(col->hdr_key, "Transcript") == 0 ) 
            BRANCH("Transcript","##INFO=<ID=Transcript,Number=A,Type=String,Description=\"Transcript names\">");
        else if ( strcmp(col->hdr_key, "VarType") == 0 ) 
            BRANCH("VarType", "##INFO=<ID=VarType,Number=A,Type=String,Description=\"Variant type.\">");
        else if ( strcmp(col->hdr_key, "ExonIntron") == 0 ) 
            BRANCH("ExonIntron","##INFO=<ID=ExonIntron,Number=A,Type=String,Description=\"Exon/CDS or intron id on transcripts.\">");
        else if ( strcmp(col->hdr_key, "IVSnom") == 0 ) 
            BRANCH("IVSnom", "##INFO=<ID=IVSnom,Number=A,Type=String,Description=\"Old style nomenclature for the description of intron variants. Not recommand to use it.\">");
        else if ( strcmp(col->hdr_key, "Oldnom") == 0 ) 
            BRANCH("Oldnom", "##INFO=<ID=Oldnom,Number=A,Type=String,Description=\"Old style nomenclature, compared with HGVSnom use gene position instead of UTR/coding position.\">");
        else if ( strcmp(col->hdr_key, "AAlength") == 0 ) 
            BRANCH("AAlength", "##INFO=<ID=AAlength,Number=A,Type=String,Description=\"Amino acid length for each transcript. 0 for noncoding transcript.\">");
    }
#undef BRANCH

    f->h = hgvs_handler_init(rna, data, reference);
    return f;
}
void anno_hgvs_core(struct anno_hgvs_file *f, bcf_hdr_t *hdr, bcf1_t *line)
{
    anno_hgvs_update_buffer(f, hdr, line);
    anno_hgvs_setter_hgvsnom(f, hdr, line);    
}


#ifdef ANNO_HGVS_MAIN
#include "anno_thread_pool.h"
#include "anno_pool.h"
#include "number.h"
#include <unistd.h>

int usage()
{
    fprintf(stderr, "anno_hgvs [options] in.vcf\n");
    fprintf(stderr, " -data <genepred_plus.gz>   GenePredPlus database.\n");
    fprintf(stderr, " -rna  <rna.fa>             RNA sequence in fasta format.\n");
    fprintf(stderr, " -ref  <reference.fa>       Genome reference sequence in fasta format.\n");
    fprintf(stderr, " -tag <tag,tag>             Specify tags.\n");    
    fprintf(stderr, " -t [1]                     Threads.\n");
    fprintf(stderr, " -O <u|v|b|z>               Output format.\n");
    fprintf(stderr, " -o <output.vcf>            Output file.\n");
    fprintf(stderr, " -r [1000]                  Record per thread per time.\n");
    fprintf(stderr, " -gene                      Selected gene list for annotation.\n");
    fprintf(stderr, " -transcript                Selected transcript list for annotation.\n");
    return 1;
}
static const char *hts_bcf_wmode(int file_type)
{
    if ( file_type == FT_BCF ) return "wbu";    // uncompressed BCF
    if ( file_type & FT_BCF ) return "wb";      // compressed BCF
    if ( file_type & FT_GZ ) return "wz";       // compressed VCF
    return "w";                                 // uncompressed VCF
}
struct args {
    const char *input_fname;
    const char *output_fname;
    const char *data_fname;
    const char *rna_fname;
    const char *reference_fname;
    int n_thread;
    htsFile *fp_input;
    htsFile *fp_out;
    bcf_hdr_t *hdr_out;
    int n_record;
    void    *gene_hash;
    void    *trans_hash;
    struct anno_hgvs_file **files;
} args = {
    .input_fname  = NULL,
    .output_fname = NULL,
    .data_fname   = NULL,
    .rna_fname    = NULL,
    .reference_fname = NULL,
    .n_thread     = 1,
    .files        = NULL,
    .fp_input     = NULL,
    .fp_out       = NULL,
    .hdr_out      = NULL,
    .n_record     = 100,
    .gene_hash    = NULL,
    .trans_hash   = NULL,
};

void *anno_hgvs(void *arg, int idx) {
    struct anno_pool *pool = (struct anno_pool*)arg;
    struct args *args = (struct args*)pool->arg;
    struct anno_hgvs_file *f = args->files[idx];
    int i;
    for ( i = 0; i < pool->n_reader; ++i) {
        bcf_unpack(pool->readers[i], BCF_UN_INFO);
        anno_hgvs_core(f, args->hdr_out, pool->readers[i]);
    }
    return pool;
}

int parse_args(int argc, char **argv)
{
    int i;
    if ( argc == 1 )
        return usage();

    const char *tags        = 0;
    const char *thread      = 0;
    const char *record      = 0;
    const char *output_type = 0;
    const char *gene_list_fname = 0;
    const char *trans_list_fname = 0;
    
    for ( i = 1; i < argc; ) {
        
        const char *a = argv[i++];
        const char **var = 0;
        
        if ( strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0 ) return usage();
        else if ( strcmp(a, "-data") == 0 ) var = &args.data_fname;
        else if ( strcmp(a, "-rna") == 0 ) var = &args.rna_fname;
        else if ( strcmp(a, "-tag") == 0 || strcmp(a, "-column") == 0) var = &tags;
        else if ( strcmp(a, "-t") == 0 ) var = &thread;
        else if ( strcmp(a, "-o") == 0 ) var = &args.output_fname;
        else if ( strcmp(a, "-O") == 0 ) var = &output_type;
        else if ( strcmp(a, "-r") == 0 ) var = &record;
        else if ( strcmp(a, "-gene") == 0 || strcmp(a, "-data") == 0 ) var = &gene_list_fname;
        else if ( strcmp(a, "-trans") == 0 ) var = &trans_list_fname;
        else if ( strcmp(a, "-ref") == 0 ) var = &args.reference_fname;
        if ( var != 0 ) {
            if ( i == argc) error("Missing an argument after %s", a);
            *var = argv[i++];
            continue;
        }        
        if ( args.input_fname != 0 ) error("Unknown argument, %s", a);
        if ( a[0] == '-' && a[1] ) error("Unknown parameter, %s", a);
        args.input_fname = a;
    }
    
    if ( args.data_fname == 0 ) error("Please specify bed database with -data.");
    if ( args.reference_fname == 0 ) error("Please specify genome reference database with -ref.");
    if ( args.rna_fname == 0 ) error("Please specify rna reference database with -rna.");
    if ( args.input_fname == 0 && (!isatty(fileno(stdin))) ) args.input_fname = "-";
    if ( args.input_fname == 0 ) error("No input file.");
    args.fp_input = hts_open(args.input_fname, "r");
    if ( args.fp_input == NULL ) error("%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(args.fp_input);
    if ( type.format != vcf && type.format != bcf ) error("Unsupport input format, only accept VCF/BCF.");

    int out_type = FT_VCF;
    if ( output_type != 0 ) {
        switch (output_type[0]) {
            case 'b': out_type = FT_BCF_GZ; break;
            case 'u': out_type = FT_BCF; break;
            case 'z': out_type = FT_VCF_GZ; break;
            case 'v': out_type = FT_VCF; break;
            default: error("The output type \"%d\" unrecognised", out_type);
        }
    }
    args.fp_out = args.output_fname == 0 ? hts_open("-", hts_bcf_wmode(out_type)) : hts_open(args.output_fname, hts_bcf_wmode(out_type));
    
    if ( thread ) args.n_thread = str2int((char*)thread);
    if ( record ) args.n_record = str2int((char*)record);
    if ( args.n_thread < 1) args.n_thread = 1;
    if ( args.n_record < 1 ) args.n_record = 100;

    args.hdr_out = bcf_hdr_read(args.fp_input);

    if ( args.hdr_out == NULL ) error("Failed to parse header of input.");

    if ( gene_list_fname ) {
        args.gene_hash = name_hash_init(gene_list_fname);
        if ( args.gene_hash == NULL ) warnings("Failed to load gene list. %s", gene_list_fname);
    }
    if ( trans_list_fname ) {
        args.trans_hash = name_hash_init(trans_list_fname);
        if ( args.trans_hash == NULL ) warnings("Failed to load transcript list. %s", trans_list_fname);
    }
    
    args.files = malloc(args.n_thread*sizeof(void*));
    args.files[0] = anno_hgvs_file_init(args.hdr_out, tags, args.data_fname, args.rna_fname, args.reference_fname);
    //args.files[0]->gene_hash = args.gene_hash;
    //args.files[0]->trans_hash = args.trans_hash;
    
    for ( i = 1; i < args.n_thread; ++i ) args.files[i] = anno_hgvs_file_duplicate(args.files[0]);
    bcf_hdr_write(args.fp_out, args.hdr_out);
    return 0;
}

int bcfanno_hgvs()
{
    struct thread_pool *p = thread_pool_init(args.n_thread);
    struct thread_pool_process *q = thread_pool_process_init(p, args.n_thread*2, 0);
    struct thread_pool_result  *r;
    
    for ( ;; ) {
        struct anno_pool *arg = anno_reader(args.fp_input, args.hdr_out, args.n_record);
        if ( arg->n_reader == 0 )
            break;
        arg->arg = &args;
        int block;
        do {
            block = thread_pool_dispatch2(p, q, anno_hgvs, arg, 1);
            if ( ( r = thread_pool_next_result(q) ) ) {
                // generate output
                struct anno_pool *data = (struct anno_pool*)r->data;
                int i;
                for ( i = 0; i < data->n_reader; ++i ) {
                    bcf_write1(args.fp_out, args.hdr_out, data->readers[i]);
                    bcf_destroy(data->readers[i]);
                }
                free(data->readers);
                thread_pool_delete_result(r, 1);
            }
            // flush output
        } while ( block == -1);
    }

    thread_pool_process_flush(q);
    while ( (r = thread_pool_next_result(q)) ) {
        // generate output
        struct anno_pool *data = (struct anno_pool*)r->data;
        int i;
        for ( i = 0; i < data->n_reader; ++i ) {
            bcf_write1(args.fp_out, args.hdr_out, data->readers[i]);
            bcf_destroy(data->readers[i]);
        }
        free(data->readers);
        thread_pool_delete_result(r, 1);
    }
    thread_pool_process_destroy(q);
    thread_pool_destroy(p);
    return 0;
}

void release_memory()
{
    hts_close(args.fp_input);
    hts_close(args.fp_out);
    if ( args.gene_hash ) name_hash_destroy(args.gene_hash);
    if ( args.trans_hash ) name_hash_destroy(args.trans_hash);
    bcf_hdr_destroy(args.hdr_out);
    int i;
    for ( i = 0; i < args.n_thread; ++i )
        anno_hgvs_file_destroy(args.files[i]);
    free(args.files);
}

int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    bcfanno_hgvs();

    release_memory();
    return 0;
}

#endif
