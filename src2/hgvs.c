#include "utils.h"
#include "hgvs.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/vcf.h"

struct hgvs_handler *hgvs_handler_init(const char *rna_fname, const char *data_fname, const char *reference_fname)
{
    struct hgvs_handler *h = malloc(sizeof(*h));
    memset(h, 0, sizeof(*h));
    h->reference_fname = reference_fname;
    h->rna_fname = rna_fname;
    h->data_fname = data_fname;
    h->rna_fai = fai_load(rna_fname);
    
    if ( h->rna_fai == NULL ) error("Failed to load index of %s : %s.", rna_fname, strerror(errno));

    h->fp_idx = hts_open(data_fname, "r");
    if ( h->fp_idx == NULL ) error("%s : %s.", data_fname, strerror(errno));

    h->idx = tbx_index_load(data_fname);
    if ( h->idx == NULL ) error("Failed to load index of %s : %s.", data_fname, strerror(errno));

    set_format_genepredPlus();
    
    return h;
}

struct hgvs_handler *hgvs_handler_duplicate(struct hgvs_handler *h)
{
    struct hgvs_handler *d = malloc(sizeof(*d));
    memset(d, 0, sizeof(*d));
    d->rna_fname = h->rna_fname;
    d->data_fname = h->data_fname;
    d->rna_fai = fai_load(d->rna_fname);
    d->idx = tbx_index_load(d->data_fname);
    d->fp_idx = hts_open(d->data_fname, "r");
    d->gene_hash = h->gene_hash;
    d->trans_hash = h->trans_hash;    
    return d;
}
void hgvs_handler_destroy(struct hgvs_handler *h)
{
    fai_destroy(h->rna_fai);
    tbx_destroy(h->idx);
    hts_close(h->fp_idx);
    int i;
    for ( i = 0; i < h->n_gene; ++i ) genepred_line_destroy(h->gls[i]);
    if ( h->n_gene) free(h->gls);
    free(h);
}

struct hgvs *hgvs_init(const char *chrom, int start, int end, char *ref, char *alt)
{
    struct hgvs *h = malloc(sizeof(*h));
    memset(h, 0, sizeof(*h));
    h->chr = chrom;
    h->start = start;
    h->end   = end;
    // if not consider alleles, only convert locations
    if ( ref == NULL && alt == NULL ) return h;
    // trim capped and tails
    int lr = strlen(ref);
    int la = strlen(alt);

    if ( ref && alt ) {
        // for GATK users, treat <NON_REF> as reference
        if ( strcmp(alt, "<NON_REF>") == 0 ) { h->type = var_type_nonref; return h;}

        for ( ;; ) {
            if ( ref == NULL || alt == NULL ) break;
            if (*ref == *alt) {
                ref++; alt++; lr--; la--;
                h->start++;
            }
            else if ( *(ref+lr-1) == *(alt+la-1) && lr > 0 && la > 0) {
                lr--, la--;
                h->end--;
            }
            else break;
        }
        //assert(lr >= 0 && la >= 0);
        if ( lr > 0) {
            h->ref = strdup(ref);
            h->ref[lr] = '\0';
        }
        if ( la > 0 ) {
            h->alt = strdup(alt);
            h->alt[la] = '\0';
        }
        if ( la == 0 ) {
            if ( lr == 0 ) h->type = var_type_ref;
            else h->type = var_type_del;
        } 
        else if ( lr == 0 ) h->type = var_type_ins;
        else {
            if ( la == 1 && lr == 1 ) h->type = var_type_snp;
            else h->type = var_type_delins;
        }
    }
    else if ( ref == NULL ) {
        h->type = var_type_ins;
        h->alt = strdup(alt);
    }
    else if ( alt == NULL ) {
        h->type = var_type_del;
        h->ref = strdup(ref);
    }
    // for insert, end smaller than start.
    if ( lr == 0 && h->end < h->start ) {
        h->end = h->start;
        h->start--;
    }
    return h;
}
void hgvs_destroy(struct hgvs *h)
{
    if ( h->ref ) free(h->ref);
    if ( h->alt ) free(h->alt);
    int i;
    for ( i = 0; i < h->n_tran; ++i) {
        struct hgvs_inf *inf = &h->trans[i].inf;
        struct hgvs_type *type = &h->trans[i].type;
        if ( inf->transcript ) free(inf->transcript);
        if ( inf->gene ) free(inf->gene);
        if ( type->aminos ) free(type->aminos);        
    }
    free(h->trans);
    free(h);
}
static int hgvs_handler_fill_buffer(struct hgvs *n, struct hgvs_handler *h)
{
    extern int name_hash_key_exists(void *hash, char *key);
    
    //struct genepred_line *gl = genepred_retrieve_region(h->idx, n->chr, n->start, n->end);
    //if ( gl == NULL )
    //return 0;
    
    int i;
    // clean buffer
    if ( h->n_gene > 0 ) {
        for ( i = 0; i < h->n_gene; ++i ) genepred_line_destroy(h->gls[i]);
        free(h->gls);
        h->gls = NULL;
    }
    h->n_gene = 0;
    
    int id;
    id = tbx_name2id(h->idx, n->chr);
    if ( id == -1 ) return 0;

    // retrieve genepred records from database
    hts_itr_t *itr = tbx_itr_queryi(h->idx, id, n->start-1, n->end);
    kstring_t string = {0,0,0};
    struct genepred_line *head = NULL;
    struct genepred_line *temp = NULL;
    while ( tbx_itr_next(h->fp_idx, h->idx, itr, &string) >= 0 ) {
        struct genepred_line *line = genepred_line_create();
        if ( parse_line(&string, line) ) continue;
        if ( head == NULL ) head = line;
        if ( temp ) temp->next = line;
        temp = line;
        string.l = 0;
    }  
    free(string.s);
    tbx_itr_destroy(itr);
    int l = count_list(head);
    if ( l == 0 || head == NULL ) return 0; // if no record retrieved    
    h->gls = malloc(l*sizeof(struct genepred_line));
    memset(h->gls, 0, sizeof(struct genepred_line)*l);
    for ( i = 0; i < l; ++i ) {
        assert(head != NULL);
        
        if ( h->gene_hash || h->trans_hash ) {
            // if defined transcript or gene list, only export predefined records, these usually used to emit redundant transcripts (XM_ likes)
            if ( h->gene_hash && name_hash_key_exists(h->gene_hash, head->name2) ) goto refresh_buffer;
            if ( h->trans_hash && name_hash_key_exists(h->trans_hash, head->name1) ) goto refresh_buffer;
            struct genepred_line *old = head;
            head = head->next;
            genepred_line_destroy(old);
            continue;
        }
      refresh_buffer:
        h->gls[h->n_gene++] = head;
        head = head->next;
    }

    // release overhang memory
    if ( h->n_gene == 0 && l > 0 ) free(h->gls);
    
    return h->n_gene;
}
static int find_the_block(struct genepred_line *l, int *s, int *e, int pos)
{
    *s = 0; *e = 0;
    int i;
    for ( i = 0; i < l->exon_count; ++i ) {
        int start = l->exons[BLOCK_START][i];
        int end = l->exons[BLOCK_END][i];
        if ( pos < start) { *e = i; break; }
        else {
            *s = i;
            if ( pos <= end ) { *e = i; break; }
        }
    }
    // out of range
    if ( i == l->exon_count && *s > 0 && *e == 0 ) return 1;
    return 0;
}
// 0 on success
// 1 on out of range
static int find_locate(struct genepred_line *l, int *pos, int *offset, int start, int *id)
{
    int b1 = -1, b2 = -1;
    *pos = 0;
    *offset = 0;
    // for location out of range set pos and offset to 0, annotate as '?'
    if ( find_the_block(l, &b1, &b2, start ) ) { *pos = 0, *offset = 0; return 1;}
    
    if ( b1 == b2 ) {
        if ( l->strand == '+' ) *pos = l->loc[BLOCK_START][b1] + start - l->exons[BLOCK_START][b1];
        else *pos = l->loc[BLOCK_END][b1] + (l->exons[BLOCK_END][b1] - start);
        *offset = 0;
    }
    else {
        int upstream =  start - l->exons[BLOCK_END][b1];
        int downstream = l->exons[BLOCK_START][b2]-start;
        if ( upstream > downstream ) {
            *pos = l->loc[BLOCK_START][b2];
            *offset = l->strand == '+' ? -downstream : downstream;
        }
        else {
            *pos = l->loc[BLOCK_END][b1];
            *offset = l->strand == '+' ? upstream : -upstream;
        }
    }
    // exon ID
    *id = b1;

    // adjust position if located in exon and there are realignments in the transcript sequence
    if ( *offset != 0 ) return 0;
    if ( l->n_cigar == 0 ) return 0;
    if ( l->n_cigar == 1 && (l->cigars[0] & GENEPRED_CIGAR_MATCH_TYPE) ) return 0;

    int adjust = 0;
    int i;
    int match = 0;
    int ins, del;
    for ( i = 0; i < l->n_cigar; ++i ) {            
        if ( l->cigars[i] & GENEPRED_CIGAR_MATCH_TYPE ) match += l->cigars[i] >> GENEPRED_CIGAR_PACKED_FIELD;            
        else if ( l->cigars[i] & GENEPRED_CIGAR_DELETE_TYPE ) {
            del = l->cigars[i] >> GENEPRED_CIGAR_PACKED_FIELD;
            // there is no need to check match and *pos, becase for plus strand match will be always smaller than *pos
            // check if this deletion in the target block
            if ( l->loc[l->strand == '+' ? BLOCK_START : BLOCK_END][i] <= match && *pos > match) {
                adjust -= del;            
                // if this variant located in the deletion, just put pos to the edge of this gap
                if ( *pos < match +  del ) { *pos += 1; return 0; }
            }
            match += del;
        }
        else if ( l->cigars[i] & GENEPRED_CIGAR_INSERT_TYPE ) {
            ins = l->cigars[i] >> GENEPRED_CIGAR_PACKED_FIELD;
            if ( l->loc[l->strand == '+' ? BLOCK_START : BLOCK_END][b1] <= match  && *pos > match ) adjust += ins;
        }
        
        if ( l->strand == '+' && match >= *pos ) break;            
        else if ( l->strand == '-' && match < *pos ) break;
    }
    *pos += adjust;
    return 0;
}
//static int check_func_vartype(struct hgvs_handler *hand, struct hgvs_inf *inf, struct hgvs_type *type, struct genepred_line *gl)
// n : iterate of transcript
static int check_func_vartype(struct hgvs_handler *h, struct hgvs *hgvs, int n, struct genepred_line *gl)
{
    struct hgvs_inf *inf = &hgvs->trans[n].inf;
    struct hgvs_type *type = &hgvs->trans[n].type;
    memset(type, 0, sizeof(struct hgvs_type));

    int pos = inf->pos;
    int ref_length = hgvs->ref == NULL ? 0 : strlen(hgvs->ref);
    int alt_length = hgvs->alt == NULL ? 0 : strlen(hgvs->alt);
    char *ref_seq = hgvs->ref == NULL ? NULL : strdup(hgvs->ref);
    char *alt_seq = hgvs->alt == NULL ? NULL : strdup(hgvs->alt);

    // for reverse strand, complent sequence
    if ( inf->strand == '-' ) {        
        if ( ref_seq ) compl_seq(ref_seq, ref_length);
        if ( alt_seq ) compl_seq(alt_seq, alt_length);
    }

    if (gl->cdsstart == gl->cdsend ) type->func1 = type->func2 = func_region_noncoding;
    // if coding transcript, check UTR or coding regions
    else {
        if ( inf->pos <= gl->utr5_length ) type->func1 = func_region_utr5;
        else if ( inf->pos > gl->cds_length ) type->func1 = func_region_utr3;
        else type->func1 = func_region_cds;
        
        if ( inf->pos == inf->end_pos && inf->offset == inf->end_offset) type->func2 = type->func1;
        else if ( type->func2 == func_region_unknown ) { // in case out of range
            if ( inf->end_pos <= gl->utr5_length ) type->func2 = func_region_utr5;
            else if ( inf->end_pos > gl->cds_length ) type->func2 = func_region_utr3;
            else type->func2 = func_region_cds;
        }
    }

    // Convert start and end  location variant from genepred structure and genome coordinate.
    // start
    if ( type->func1 == func_region_utr5 ) inf->loc = gl->utr5_length - inf->pos +1;
    else if ( type->func1 == func_region_utr3 ) inf->loc = inf->pos - gl->cds_length;
    else if ( type->func1 == func_region_cds ) inf->loc = inf->pos - gl->utr5_length;
    else inf->loc = inf->pos;
    //end
    if ( type->func2 == func_region_utr5 ) inf->end_loc = gl->utr5_length - inf->end_pos +1;
    else if ( type->func2 == func_region_utr3 ) inf->end_loc = inf->end_pos - gl->cds_length;
    else if ( type->func2 == func_region_cds ) inf->end_loc = inf->end_pos - gl->utr5_length;
    else inf->end_loc = inf->end_pos;

    // check exon,intron,cds ID, check the splice sites
    int i; // exon / intron ID
    int c; // cds id
    // Variant type should be consider as function variant type (vartype1) and splice site (vartype2);    
    type->vartype = var_is_unknown;
    type->vartype2 = var_is_not_splice;
    
    // Only check the start of variants
    // Count the Exome/Intron and CDS id, for minus strand id should count from backward
    // Only check the start of variants. pos == start of variant 
    if ( gl->strand == '+' ) {
        for ( i = 0, c = 0; i < gl->exon_count; ) {
            if ( gl->loc[BLOCK_END][i] > gl->utr5_length && pos > gl->utr5_length && pos <= gl->cds_length ) c++;
            if ( pos >= gl->loc[BLOCK_START][i] && pos <= gl->loc[BLOCK_END][i]) {
                if ( inf->offset == 0 ) {
                    if ( pos <= gl->loc[BLOCK_START][i] + SPLICE_SITE_EXON_RANGE || pos >= gl->loc[BLOCK_END][i] - SPLICE_SITE_EXON_RANGE)
                        type->vartype2 = var_is_splice_site;
                    else if ( pos + ref_length <= gl->loc[BLOCK_START][i] + SPLICE_SITE_EXON_RANGE || pos + ref_length >= gl->loc[BLOCK_END][i] - SPLICE_SITE_EXON_RANGE )
                        type->vartype2 = var_is_splice_site;
                }
                break;
            }
            ++i;         
        }
        if ( inf->offset == 0 ) { type->count = i+1; type->count2 = c; }
        else { type->count = inf->offset < 0 ? i : i +1; type->count2 = 0; }
    }
    else { // for minus strand, count CDS and EXON id backward
        for ( i = 0, c = 0; i < gl->exon_count; ) {
            int j = gl->exon_count -i -1;
            if ( gl->loc[BLOCK_START][j] > gl->utr5_length && pos > gl->utr5_length && pos <= gl->cds_length) c++;
            if ( pos >= gl->loc[BLOCK_END][j] && pos <= gl->loc[BLOCK_START][j] ) {
                if ( inf->offset == 0 ) {
                    if (pos <= gl->loc[BLOCK_END][j] + SPLICE_SITE_EXON_RANGE || pos >= gl->loc[BLOCK_START][j] - SPLICE_SITE_EXON_RANGE)
                        type->vartype2 = var_is_splice_site;
                    // in case indels cover splice site
                    else if ( pos + ref_length <= gl->loc[BLOCK_END][j] + SPLICE_SITE_EXON_RANGE || pos + ref_length >= gl->loc[BLOCK_START][j] - SPLICE_SITE_EXON_RANGE ) 
                        type->vartype2 = var_is_splice_site;
                }
                break;
            }
            ++i;
        }
        if ( inf->offset == 0 ) { type->count = i+1; type->count2 = c; }
        else { type->count = inf->offset < 0 ? i : i + 1; type->count2 = 0; }
    }

    
    // check if variantion located in the splice sites around the edge of utr and cds regions    
    if ( inf->offset < 0 ) {
        type->vartype = var_is_intron;
        if ( inf->offset > -SPLICE_SITE_EXON_RANGE ) type->vartype2 = var_is_splice_acceptor;
        goto no_amino_code;
    }
    else if ( inf->offset > 0 ) {
        type->vartype = var_is_intron;
        if ( inf->offset < SPLICE_SITE_EXON_RANGE ) type->vartype2 = var_is_splice_donor;
        goto no_amino_code;
    }
    else if ( pos > gl->utr5_length - SPLICE_SITE_EXON_RANGE && pos < gl->utr5_length + SPLICE_SITE_EXON_RANGE )
        type->vartype2 = var_is_splice_site; // init
    else if ( pos > gl->cds_length - SPLICE_SITE_EXON_RANGE && pos < gl->cds_length + SPLICE_SITE_EXON_RANGE )
        type->vartype2 = var_is_splice_site; // stop    
    
#define BRANCH(_type) do {                              \
        if ( type->vartype == var_is_unknown ) {        \
            type->vartype = _type;                      \
        }                                               \
} while(0)

    // Check if noncoding transcript.
    if ( gl->cdsstart == gl->cdsend ) {
        // no cds count for noncoding transcript
        type->count2 = 0;
        BRANCH(var_is_noncoding);
        goto no_amino_code;
    }

    //  deletion variant from UTR5 that changes at least one base of the canonical start codon, trim the UTR part
    // int deletion_variant_cover_utr5_length = 0;
    if ( type->func1 == func_region_utr5 && type->func2 != func_region_utr5 ) {

        // Do not try to predict the amnino acid changes if partly or whole start codon deleted

        BRANCH(var_is_transcript_ablation);
        if ( type->func2 == func_region_cds ) {
            //deletion_variant_cover_utr5_length = ref_length - inf->end_loc;            
            //assert(deletion_variant_cover_utr5_length > 0);
            type->vartype2 = var_is_splice_site; // for variants covered start codon, should always be splice site
            type->vartype  = var_is_start_lost;
        }
        goto no_amino_code;
    }
    else {
        //  Check if UTR regions.  
        if ( pos <= gl->utr5_length -SPLICE_SITE_EXON_RANGE ) { BRANCH(var_is_utr5); goto no_amino_code; }
        if ( pos <= gl->utr5_length ) { type->vartype2 = var_is_splice_site; goto no_amino_code; }  // init
        if ( pos >= gl->cds_length + SPLICE_SITE_EXON_RANGE ) { BRANCH(var_is_utr3); goto no_amino_code;}
        if ( pos >= gl->cds_length ) { type->vartype2 = var_is_splice_site; goto no_amino_code; }
    }
    
    // next we will check amino acid changes
    if ( inf->offset != 0 ) goto no_amino_code;
    
    // For variants in coding region, check the amino acid changes.
    int cds_pos = pos - gl->utr5_length;
    // variant start in the amino codon, 0 based, [0,3)
    int cod = (cds_pos-1) % 3;
    // codon inition position in the transcript, 0 based.
    // NOTICE: One to 3 base(s) will be CAPPED for the start of insertion; remove caps to check amino acid changes.
    //int start = cod == 0 ? gl->utr5_length + 1 : (cds_pos-1)/3*3 + gl->utr5_length;
    int start = (cds_pos-1)/3*3 + gl->utr5_length;
    // retrieve affected sequences and downstream    
    int transcript_retrieve_length = 0;
    char *name = gl->name1;
    // original sequnence from RNA sequence, start round from first aa changed position;
    // start aa location may be changed becase realignment, but ori_seq will be "stable" (reset if duplicate) in this function;
    // ori_seq is a temp sequence, will be free before level this function.
    char *ori_seq = NULL;
    ori_seq = faidx_fetch_seq(h->rna_fai, name, start, start + 10000, &transcript_retrieve_length);
    if ( ori_seq == NULL || transcript_retrieve_length == 0 ) goto failed_check;

    // Sometime trancated transcript records will disturb downstream analysis.
    if ( transcript_retrieve_length < 3 ) {
        warnings("Record %s probably trancated.", name);
        goto failed_check;
    }
        
    
    // Check if insert bases in tandam short repeats region, amino acids will be changed in the downstream;
    // Realign the alternative allele ..
    // check if and only if in coding region !!!
    if ( ref_length == 0 && alt_length > 0) {
        int offset = 0;
        char *ss = ori_seq;
        // remove caps
        ss += cod+1;
        
        while ( same_DNA_seqs(alt_seq, ss, alt_length) == 0 ) { ss += alt_length; offset += alt_length; }

        if ( offset > 0 ) {
            if ( transcript_retrieve_length <= offset ) goto failed_check;
            // update dup_offset for original inf struct, check this value when output hgvs position
            inf->dup_offset = offset;

            cds_pos += offset;
            cod = (cds_pos-1)%3;
            int new_start = (cds_pos-1)/3*3 + gl->utr5_length;
            assert(new_start >= start);
            int start_offset = new_start -start;
            char *offset_seq = strdup(ori_seq+start_offset);
            transcript_retrieve_length -= start_offset;
            offset_seq[transcript_retrieve_length] = '\0';
            // point ori_seq to new memory address, this is not safe!
            free(ori_seq);
            ori_seq = offset_seq;
            start = new_start;
        }
    }
    
    char codon[4];
    memcpy(codon, ori_seq, 3);
    codon[3] = '\0';
    type->ori_amino = codon2aminoid(codon);
    type->loc_amino = (cds_pos-1)/3 + 1;

    // reset type    
    if ( type->n ) { free(type->aminos); type->n = 0; }

    //if ( ref_length == 1 && ref_length == alt_length ) {
    if ( hgvs->type == var_type_snp ) {
        // Check the ref sequence consistant with database or not. If variants are same with transcript sequence, set type to no_call.
        if ( seq2code4(ori_seq[cod]) != seq2code4(ref_seq[0]) ) {
            if (seq2code4(ori_seq[cod]) == seq2code4(alt_seq[0]) ) type->vartype = var_is_no_call;
            else warnings("Inconsistance nucletide: %s,%d,%c vs %s,%d,%c.", hgvs->chr, hgvs->start, ref_seq[0], gl->name1, pos, ori_seq[cod]);
        }

        codon[cod] = *alt_seq;
        type->mut_amino = codon2aminoid(codon);
        if ( type->ori_amino == type->mut_amino ) {
            if ( type->ori_amino == 0 ) type->vartype = var_is_stop_retained;
            else BRANCH(var_is_synonymous);
        }
        else {
            if ( type->ori_amino == 0 ) type->vartype = var_is_stop_lost;
            else if ( type->mut_amino == 0 ) type->vartype = var_is_nonsense;
            else BRANCH(var_is_missense);
        }
    }
    // if insert or delete
    //else {
    else if ( hgvs->type == var_type_ins ) {
        // Insertions
            
        // if insertion break the frame goto check delins
        if ( alt_length%3 ) goto delins;
        // if insertion break the original codon goto check delins
        if ( cod != 2 ) goto delins;
        
        // Check if insertion does NOT change the amino acid, treat as inframe insertion, else go to delins
        // memcpy(codon, ori_seq, 3);

        // ?
        //if ( type->ori_amino != codon2aminoid(codon) ) goto delins;
        
        BRANCH(var_is_inframe_insertion);
        
        // end amino acid
        memcpy(codon, ori_seq+3, 3);
        type->ori_end_amino = codon2aminoid(codon);
        type->loc_end_amino = type->loc_amino + 1;
        
        // free aminos buffer before reallocated it, previous memory leaks
        type->n = alt_length/3;
        type->aminos = malloc(type->n *sizeof(int));

        int i;
        for ( i = 0; i < alt_length/3; ++i ) {
            memcpy(codon, alt_seq+i*3, 3);
            type->aminos[i] = codon2aminoid(codon);
        }
    } 
    // Deletion
    else if ( hgvs->type == var_type_del ) {
        // for an exon-span deletion, the transcript has been trimmed, report stop-lost?
        if ( ref_length >= transcript_retrieve_length ) {
            type->vartype = var_is_stop_lost;
            goto no_amino_code;
        }
        
        if ( ref_length%3 ) goto delins;

        if ( cod != 2 ) goto delins;
        
        BRANCH(var_is_inframe_deletion);
            
        if ( ref_length == 3 ) {
            type->loc_end_amino = type->loc_amino;
            type->ori_end_amino = codon2aminoid(ori_seq);
            //type->mut_amino = codon2aminoid(ori_seq);
        }
        else {
            memcpy(codon, ori_seq + ref_length-3, 3);            
            type->ori_end_amino = codon2aminoid(codon);
            type->loc_end_amino = type->loc_amino + ref_length/3;        
            type->n = ref_length/3;
            type->aminos = malloc(type->n *sizeof(int));
            int i;
            for ( i = 0; i < ref_length/3; ++i ) type->aminos[i] = codon2aminoid(ori_seq+i*3);
        }
    }
    // Delins
    else {
      delins:

        // assume it is a inframe delins first
        BRANCH(var_is_inframe_delins);
        
        type->loc_end_amino = (cds_pos+ref_length-1)/3+1;
        type->ori_end_amino = codon2aminoid(ori_seq+(type->loc_end_amino - type->loc_amino)*3);
        
        if ( alt_length == ref_length) {
            char *new_seq_p = ori_seq;
            
            // if inframe delins break one amino acid frame, try to trim both ends of amino acid sequences and predict alternate aa
            if ( type->loc_end_amino > type->loc_amino ) {
                int i, j;
                char *alt_allele = strndup(ori_seq, ref_length+4);
                for ( i = 0, j = cod; i < ref_length && j < ref_length; ++i ) alt_allele[++j] = alt_seq[i];
                
                type->n = type->loc_end_amino - type->loc_amino + 1;

                // trim ends of amino acid sequences
                int offset_init_codon = 0;
                int offset_tail_codon = 0;
                // check start aa
                for ( ; offset_init_codon < type->n; ) {
                    if ( codon2aminoid(alt_allele+offset_init_codon*3) == codon2aminoid(ori_seq+offset_init_codon*3) ) offset_init_codon++;
                    break;
                }
                if ( offset_init_codon > 0 ) {
                    type->loc_amino += offset_init_codon;
                    type->ori_amino = codon2aminoid(ori_seq+i*3);
                    new_seq_p = ori_seq + offset_init_codon*3;
                    assert ( type->loc_end_amino >= type->loc_amino );
                }
                // check end aa
                for ( ; offset_tail_codon < type->n-offset_init_codon-1; ) {
                    if ( codon2aminoid(alt_allele+(type->n-1-offset_tail_codon)*3) == codon2aminoid(ori_seq+(type->n-1-offset_tail_codon)*3) ) offset_tail_codon++;
                    break;
                }
                if ( offset_tail_codon > 0 ) {
                    type->loc_end_amino -= offset_tail_codon;
                    type->ori_end_amino = codon2aminoid(ori_seq+(type->n-1-offset_tail_codon)*3);
                }

                if ( alt_allele ) free(alt_allele);
            }
            type->n = type->loc_end_amino - type->loc_amino + 1;
            type->aminos = (int*)malloc(sizeof(int)*type->n);
            for ( i = 0; i < type->n; ++i ) type->aminos[i] = codon2aminoid(new_seq_p+i*3);
            // repredict the change type if only one aa changed
            if ( type->n == 1 ) {
                type->mut_amino = type->aminos[0];
                if ( type->ori_amino == type->mut_amino ) {
                    if ( type->ori_amino == X_CODO ) type->vartype = var_is_stop_retained;
                    else type->vartype = var_is_synonymous;
                }
                else {
                    if ( type->ori_amino == X_CODO ) type->vartype = var_is_stop_lost;
                    else if ( type->mut_amino == X_CODO) type->vartype = var_is_nonsense;
                    else type->vartype = var_is_missense;
                }
                free(type->aminos);
                type->aminos = NULL;
                type->n = 0;
            }
        }
        // inframe insertion
        else if ( ref_length == 0 && alt_length %3 == 0 ) {
            // if insert codons without check orignal codon will interpret as inframe_insertion
            // else if orignal codon changed, interpret as inframe_delins
            kstring_t str = {0,0,0};
            kputsn(ori_seq, cod+1, &str);
            kputs(alt_seq, &str);
            kputsn(ori_seq+cod+1, 3-cod-1, &str);
            assert(str.l%3 == 0);
            type->n = str.l/3;
            
            // if inserted base disrupt original aa, 2 nearby aa will be checked
            int i = 0, j;
            if ( cod != 2 && codon2aminoid(ori_seq) == codon2aminoid(str.s) ) {
                // convert p.XXdelinsXX ==> p.XX_XXinsXX
                type->vartype = var_is_inframe_insertion;
                i = 1;
                type->loc_end_amino = type->loc_amino+1;
                type->ori_end_amino = codon2aminoid(ori_seq+3);
            }
            type->n -= i;
            type->aminos = malloc(sizeof(int)*(type->n));            
            for ( j = 0;j < type->n; ++i,++j) type->aminos[j] = codon2aminoid(str.s+3*i);
            if ( str.m ) free(str.s);
        }
        // inframe deletion
        else if ( alt_length == 0 && ref_length %3 == 0 ) {
            
            if ( ref_length < transcript_retrieve_length - cod ) {
                // cod == 0 for first base of codon, cod > 0 indicate deletion breakup original codon, treat it as delins
                if ( cod == 0 ) {
                    type->n = ref_length/3;
                    type->aminos = malloc(type->n*sizeof(int));
                    int i;
                    for (i=0; i<type->n; ++i) type->aminos[i] = codon2aminoid(ori_seq+i*3);
                }
                else {
                    memcpy(codon, ori_seq, cod);
                    memcpy(codon+cod, ori_seq+ref_length+cod, 3-cod);
                    // original nearby codons break, connect trancated codons on both ends into a new codon 
                    type->n = 1;
                    type->aminos = (int*)malloc(sizeof(int));
                    type->aminos[0] = codon2aminoid(codon);
                    // there is no need to update ori_seq in this function
                    // if first codon or last codon same with alternative codon, trim it
                    if ( type->ori_amino == type->aminos[0] ) {
                        type->vartype = var_is_inframe_deletion;
                        type->loc_amino++;
                        type->ori_amino = codon2aminoid(ori_seq+3);
                        type->n = 0;
                        free(type->aminos);
                        type->aminos = NULL;
                    }
                    else if ( type->ori_end_amino == type->aminos[0] ) {
                        type->vartype = var_is_inframe_deletion;
                        type->loc_end_amino--;
                        // !!! type->loc_amino should inited and NOT changed again before enter this function
                        int l = (type->loc_end_amino - type->loc_amino)*3;
                        type->ori_end_amino = codon2aminoid(ori_seq+l);
                        type->n = 0;
                        free(type->aminos);
                        type->aminos = NULL;
                    }
                }
            }
            // if deletion break the end of transcription..
            else {
                goto variant_frameshift;
            }
        }
        // frameshift
        else {
          variant_frameshift:
            type->vartype = var_is_frameshift;
            kstring_t str = {0,0,0};
            int ori_stop = gl->cds_length - type->loc_amino;
            char codon[4];
            memcpy(codon, ori_seq, 3);
            type->ori_amino = codon2aminoid(codon);
                
            if ( cod > 0 ) kputsn(ori_seq, cod, &str);
            if ( alt_length > 0 ) kputsn(alt_seq, alt_length, &str);
            if ( ref_length < transcript_retrieve_length ) {
                kputsn(ori_seq+cod+ref_length, transcript_retrieve_length-cod-ref_length, &str);
                // delins in stop codon, and no UTR 3 in this transcript                
                if (str.l < 3) {
                    faidx_t *fai = fai_load(h->reference_fname);
                    if ( fai ) {
                        char *refseq = NULL;
                        int length = 0;
                        if ( gl->strand == '+' ) 
                            refseq = faidx_fetch_seq(fai, gl->chrom, gl->txend, gl->txend+10, &length);
                        else {
                            refseq = faidx_fetch_seq(fai, gl->chrom, gl->txstart - 10, gl->txstart, &length);
                            compl_seq(refseq, length);
                        }
                        if ( length && refseq ) {
                            kputs(refseq, &str);
                            memcpy(codon, str.s, 3);
                            type->mut_amino = codon2aminoid(codon);
                        }
                        fai_destroy(fai);
                    }
                }
                else {
                    memcpy(codon, str.s, 3);
                    type->mut_amino = codon2aminoid(codon);
                }
                
                int i = 0, j;
                // ajust amnio acid change; for some frameshift variants first aa maybe retained, we should adjust aa forward
                // until the most first aa change position. For example,
                // NM_014638.3:c.3705_3721dup(p.Leu1240Leufs*29) should be interpret as :p.(Ser1241Phefs*42)

                for ( ;; ) {
                    if (type->ori_amino != type->mut_amino || i*3>str.l) break;
                    i++;
                    memcpy(codon, ori_seq+i*3, 3);
                    type->ori_amino = codon2aminoid(codon);
                    memcpy(codon, str.s+i*3, 3);
                    type->mut_amino = codon2aminoid(codon);
                }
                for ( j = i; j < str.l/3; ++j ) {
                    if ( check_is_stop(str.s+j*3) ) break;
                }
                // frameshift -1 for no change,else for termination site
                type->fs = ori_stop == j+1-i ? -1 : j+1-i;
                // adjust location to the nearest variant position
                type->loc_amino += i;
                type->loc_end_amino += i;
                    
                if ( str.m ) free(str.s);
            } // end of frameshift
        } // end of Delins        
    } // end var type
    
#undef BRANCH
    
    if ( ref_seq != NULL ) free(ref_seq);
    if ( alt_seq != NULL ) free(alt_seq);
    if ( ori_seq != NULL ) free(ori_seq);
    return 0;
    
  failed_check:
    if ( ori_seq != NULL ) free(ori_seq);
    type->vartype = var_is_unknown;
    
  no_amino_code:
    if ( ref_seq != NULL ) free(ref_seq);
    if ( alt_seq != NULL ) free(alt_seq);

    type->loc_amino = 0;
    type->ori_amino = 0;
    type->mut_amino = 0;
    type->ori_end_amino = 0;
    type->loc_end_amino = 0;
    type->n = 0;
    type->aminos = 0;
    type->fs = 0;
    return 0;
}
int hgvs_anno_trans(struct hgvs *n, struct hgvs_handler *h)
{
    if ( hgvs_handler_fill_buffer(n, h) == 0 ) return 0;
    n->trans = malloc(h->n_gene*sizeof(struct hgvs_core));
    int i;
    for ( i = 0; i < h->n_gene; ++i ) {
        struct genepred_line *gl = h->gls[i];
        if ( gl->loc_parsed == 0 ) {
            if ( parse_line_locs(gl) ) {
                warnings("Failed to parse locs of %s", gl->name2);
                continue;
            }
        }
        struct hgvs_core *trans = &n->trans[n->n_tran];
        memset(trans, 0, sizeof(*trans));
        // update hgvs locations
        struct hgvs_inf *inf = &trans->inf;
        struct hgvs_type *type = &trans->type;
        inf->aa_length = gl->cdsstart == gl->cdsend ? 0 : (gl->cds_length - gl->utr5_length)/3;
        inf->strand = gl->strand;
        inf->version = gl->name_version;
        
        int ex1 = 0, ex2 = 0;
        if ( find_locate(gl, &inf->pos, &inf->offset, n->start, &ex1) ) continue;            
        if ( n->start != n->end ) {
            // unannotated location will export as '?'
            if ( find_locate(gl, &inf->end_pos, &inf->end_offset, n->end, &ex2) ) type->func2 = func_region_outrange;
            
            if ( gl->strand == '-' ) {
                int t = inf->pos;
                inf->pos = inf->end_pos;
                inf->end_pos = t;
                t = inf->offset;
                inf->offset = inf->end_offset;
                inf->end_offset = t;
            } 
        }
        else {
            inf->end_pos = inf->pos;
            inf->end_offset = inf->offset;
        }
        inf->transcript = strdup(gl->name1);
        inf->gene = strdup(gl->name2);

        // update variant functional types

        check_func_vartype(h, n, n->n_tran, gl);
        n->n_tran++;
    }
    return n->n_tran;    
}
