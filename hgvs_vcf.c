#include "utils.h"
#include "hgvs.h"
#include "hgvs_vcf.h"
#include "genepred.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include <string.h>


static int init_flag = 0;

int hgvs_update_vcf_header(bcf_hdr_t *hdr)
{
    int id;

    id = bcf_hdr_id2int(hdr, BCF_DT_ID, "Gene");
    if (id == -1) {
	bcf_hdr_append(hdr, "##INFO=<ID=Gene,Number=A,Type=String,Description=\"Gene names\">");
	bcf_hdr_sync(hdr);
	id = bcf_hdr_id2int(hdr, BCF_DT_ID, "Gene");
	assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));
    }    

    id = bcf_hdr_id2int(hdr, BCF_DT_ID, "Transcript");
    if (id == -1) {
	bcf_hdr_append(hdr, "##INFO=<ID=Transcript,Number=A,Type=String,Description=\"Transcript names\">");
	bcf_hdr_sync(hdr);
	id = bcf_hdr_id2int(hdr, BCF_DT_ID, "Transcript");
	assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));
    }

    id = bcf_hdr_id2int(hdr, BCF_DT_ID, "HGVSnom");
    if (id == -1) {
	bcf_hdr_append(hdr, "##INFO=<ID=HGVSnom,Number=A,Type=String,Description=\"HGVS nomenclature for the description of DNA sequence variants\">");
	bcf_hdr_sync(hdr);
	id = bcf_hdr_id2int(hdr, BCF_DT_ID, "HGVSnom");
	assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));
    }

    /* id = bcf_hdr_id2int(hdr, BCF_DT_ID, "ExIn_id"); */
    /* if (id == -1) { */
    /*     bcf_hdr_append(hdr, "##INFO=<ID=ExIn_id,Number=A,Type=String,Description=\"Exon or intron id on transcripts.\">"); */
    /*     bcf_hdr_sync(hdr); */
    /*     id = bcf_hdr_id2int(hdr, BCF_DT_ID, "ExIn_id"); */
    /*     assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id)); */
    /* } */

    /* id = bcf_hdr_id2int(hdr, BCF_DT_ID, "FlankSeq"); */
    /* if (id == -1) { */
    /*     bcf_hdr_append(hdr, "##INFO=<ID=FlankSeq,Number=1,Type=String,Description=\"Three nearby bases of current position.\">"); */
    /*     bcf_hdr_sync(hdr); */
    /*     id = bcf_hdr_id2int(hdr, BCF_DT_ID, "FlankSeq"); */
    /*     assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id)); */
    /* } */
    
    id = bcf_hdr_id2int(hdr, BCF_DT_ID, "VarType");
    if (id == -1) {
	bcf_hdr_append(hdr, "##INFO=<ID=VarType,Number=A,Type=String,Description=\"Variant type.\">");
	bcf_hdr_sync(hdr);
	id = bcf_hdr_id2int(hdr, BCF_DT_ID, "VarType");
	assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));
    }
    
    return 0;
}

int init_hgvs_anno(const char *data, const char *fasta, bcf_hdr_t *hdr)
{
    if ( init_hgvs_spec(data, fasta) )
        return 1;

    hgvs_update_vcf_header(hdr);


    // turn init flag on
    init_flag = 1;
    
    return 0;
}

int close_hgvs_anno()
{
    hgvs_spec_destroy();
    
    return 0;
}

static char *retrieve_gene_des(struct hgvs_des *des)
{
    kstring_t string = { 0, 0, 0};

    int i;
    for ( i = 0; i < des->l; ++i ) {
        struct hgvs_name *name = &des->a[i].name;
        if ( i ) kputc('|', &string);
        kputs(name->name2, &string);
    }

    return string.s;
}

static char *retrieve_trans_des(struct hgvs_des *des)
{
    kstring_t string = { 0, 0, 0};

    int i;
    for ( i = 0; i < des->l; ++i ) {
        struct hgvs_name *name = &des->a[i].name;
        if ( i ) kputc('|', &string);
        kputs(name->name1, &string);
    }
    return string.s;
}

static char *retrieve_hgvs_des(struct hgvs_des *des)
{
    kstring_t string = { 0, 0, 0};
    int i;
    for ( i = 0; i < des->l; ++i ) {
        struct hgvs_name *name = &des->a[i].name;
        struct var_func_type *type = &des->a[i].type;
        if ( i ) kputc('|', &string);
        ksprintf(&string, "%s:",name->name1);
        if ( type->func == func_region_noncoding ) {
            kputs("n.", &string);
        } else if ( type->func == func_region_cds ) {
            kputs("c.", &string);
        } else if ( type->func == func_region_utr5 ) {
            kputs("c.-", &string);
        } else if ( type->func == func_region_utr3 ) {
            kputs("c.*", &string);
        }
        ksprintf(&string, "%d", name->loc);
        if ( name->offset > 0 ) {
            ksprintf(&string, "+%d", name->offset);
        } else if (name->offset < 0) {
            ksprintf(&string, "%d", name->offset);
        }
        if ( des->start != des->end ) {
            kputc('_', &string);
            if ( type->func == func_region_utr5 ) {
                kputc('-', &string);
            } else if ( type->func == func_region_utr3 ) {
                kputc('*', &string);
            }
            ksprintf(&string, "%d", name->end_loc);
            if ( name->end_offset > 0 ) {
                ksprintf(&string, "+%d", name->end_offset);
            } else if ( name->end_offset < 0) {
                ksprintf(&string, "%d", name->end_offset);
            }            
        }

        if ( des->ref_length == 0 ) {
            ksprintf(&string, "ins%s", des->alt);
        } else if ( des->alt_length == 0 ) {
            ksprintf(&string, "del%s", des->ref);
        } else {
            ksprintf(&string, "%s>%s", des->ref, des->alt);
        }
        if ( type->loc_amino > 0 && des->type == var_type_snp ) {
                if ( type->ori_amino != type->mut_amino ) {
                    ksprintf(&string, "(p.%s%d%s)", codon_names[type->ori_amino], type->loc_amino, codon_names[type->mut_amino]);
            } else {
                kputs("(p.=)", &string);
            }
        }
    }
    return string.s;
}

static char *retrieve_vartype_des(struct hgvs_des *des)
{
    kstring_t string = { 0, 0, 0};

    int i;
    for ( i = 0; i < des->l; ++i ) {
        struct var_func_type *type = &des->a[i].type;
        if ( i ) kputc('|', &string);
        kputs(var_type_string(type->vartype), &string);
    }
    return string.s;
}
// 0 on success, 1 on failed, 2 on NOT inited.
int setter_hgvs_vcf(bcf_hdr_t *hdr, bcf1_t *line)
{
    // If not inited, skip.
    if ( init_flag == 0 )
        return 2;
    
    int i;
    kstring_t gene = { 0, 0, 0};
    kstring_t transcript = { 0, 0, 0};
    kstring_t hgvs_nom = { 0, 0, 0};
    kstring_t vartype = { 0, 0, 0};

    int is_empty = 1;
    for ( i = 1; i < line->n_allele; ++i ) {
        const char *name = bcf_hdr_id2name(hdr, line->rid);        
        setter_description(name, line->pos+1, line->d.allele[0], line->d.allele[i]);
        struct hgvs_des *des = fill_hgvs_name();
        
        if ( i > 1 ) {
            kputc(',', &gene);
            kputc(',', &transcript);
            kputc(',', &hgvs_nom);
            // kputc(',', &exin);
            // kputc(',', &flank_seq);
            kputc(',', &vartype);
        }
    
        char *gene_string = retrieve_gene_des(des);
        char *trans_string = retrieve_trans_des(des);
        char *hgvs_string = retrieve_hgvs_des(des);
        // char *exin_string = retrieve_exin_des(des);
        char *vartype_string = retrieve_vartype_des(des);

        hgvs_des_clear(des);
        if ( gene_string == NULL ) {
            kputs(".", &gene);
            kputs(".", &transcript);
            kputs(".", &hgvs_nom);
            kputs(".", &vartype);
        } else {
            is_empty = 0;
            kputs(gene_string, &gene);
            kputs(trans_string, &transcript);
            kputs(hgvs_string, &hgvs_nom);
            kputs(vartype_string, &vartype);
        
            free(gene_string);
            free(trans_string);
            free(hgvs_string);
            free(vartype_string);
        }        
    }
    if ( i > 1 && is_empty == 0) {
        bcf_update_info_string(hdr, line, "Gene", gene.s);
        bcf_update_info_string(hdr, line, "Transcript", transcript.s);
        bcf_update_info_string(hdr, line, "HGVSnom", hgvs_nom.s);
        bcf_update_info_string(hdr, line, "VarType", vartype.s);
    }

    free(gene.s);
    free(transcript.s);
    free(hgvs_nom.s);
    free(vartype.s);
    
    return 0;                  
}
