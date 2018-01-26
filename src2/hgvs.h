#ifndef HGVS_HEADER
#define HGVS_HEADER

#include "utils.h"
#include "genepred.h"
#include "variant_type.h"

#define SPLICE_SITE_EXON_RANGE      3
#define SPLICE_SITE_INTRON5_RANGE   3
#define SPLICE_SITE_INTRON3_RANGE   8
#define SPLICE_DONDOR_RANGE         2
#define SPLICE_ACCEPTOR_RANGE       2
 
enum func_region_type {
    _func_region_promote_to_int = -1,
    func_region_unknown = 0,
    //func_region_split_sites,
    func_region_cds,
    func_region_noncoding,
    func_region_intron,
    func_region_utr5,
    func_region_utr3,
    func_region_intergenic,
    // If variants across more than one exon, define large fragment.
    //func_region_large,
    func_region_outrange,
};

struct hgvs_type {
    // start
    enum func_region_type func1;
    // end
    enum func_region_type func2;
    enum var_func_type vartype;
    // update:2017/12/01
    // another type to mark splice site, in case lost variants type in splice sites
    enum var_type_splice vartype2;
    // Exon or intron count. Start from 1.
    int count;
    // CDS count, for noncoding transcript always be 0.
     int count2;
    // Location of amino first influenced. Sometime a indel may change the amino sequence after several condons.
    int loc_amino;
    // Original amino acid and mutated amino acid. Check amino type only if variants happed in cds region.
    int ori_amino;
    // the end_amino only used in delins
    int ori_end_amino;
    int loc_end_amino;
    // fast access for single amino acid change
    int mut_amino;
    // buffer for insertion and deletion
    int n;
    int *aminos;

    // Terminal sequence after loc_amino.
    int fs;
    // Extern codon, for stop-loss
    int ext;
};

struct hgvs_inf {
    // transcript name or Loc name in same case, name1
    char *transcript;
    // transcript version
    int version;
    // gene name
    char *gene;
    
    // Amino acid length, for noncoding transcript should always be 0.
    int aa_length;    
    
    // position on gene coordinate, not account the function regions.
    int pos;
    int end_pos;
    // 0 for unknown "?"
    // location for function region.
    int loc;
    int end_loc;
    // offset in intron
    int offset;
    int end_offset;

    // check if repeat sequences(duxp), only for insterion
    int dup_offset;
    
    // If strand is '-', convert sequence to complement strand. Default is '+'.
    char strand;
};

struct hgvs_core {
    struct hgvs_inf inf;
    struct hgvs_type type;
};
struct hgvs {
    // point to chrom of VCF header, do not free it
    const char *chr;
    int start;
    int end;
    char *ref;
    char *alt;
    // variant type include SNP,INDEL only accunt for allele changes
    enum variant_type type;

    int n_tran;
    struct hgvs_core *trans;
};
struct hgvs_handler {
    const char *rna_fname;
    const char *data_fname;
    const char *reference_fname;
    faidx_t *rna_fai;
    tbx_t *idx;
    htsFile *fp_idx;
    // point to gene hash
    void *gene_hash;
    // point to transcript hash
    void *trans_hash;
    int n_gene;
    struct genepred_line **gls;
};

extern struct hgvs_handler *hgvs_handler_init(const char *rna_fname, const char *data_fname, const char *reference_fname);
extern struct hgvs_handler *hgvs_handler_duplicate(struct hgvs_handler *h);
extern void hgvs_handler_destroy(struct hgvs_handler *h);

extern struct hgvs *hgvs_init(const char *chrom, int start, int end, char *ref, char *alt);
extern void hgvs_destroy(struct hgvs *h);
extern int hgvs_anno_trans(struct hgvs *n, struct hgvs_handler *h);

#endif
