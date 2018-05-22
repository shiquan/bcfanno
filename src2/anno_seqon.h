#ifndef ANNO_SEQON_HEADER
#define ANNO_SEQON_HEADER

#include "utils.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "htslib/faidx.h"
#include "anno_pool.h"
#include "anno_col.h"
#include "variant_type.h"

extern int file_is_GEA(const char *fn);

struct mc_handler {
    const char *rna_fname;
    const char *data_fname;
    const char *reference_fname;
    faidx_t *rna_fai;
    tbx_t *idx;
    htsFile *fp_idx;
    struct gea_hdr *hdr;
    
    void *name_hash;
    
    // sweep record per chunk
    // for Random Access Mode, do not require input records sorted, annotate each record one by one
    int end_pos_for_skip;
    
    int i_record;
    int n_record;
    int m_record;
    // gene and regulatory records
    void **records;

    // point to nearest gene record, used to interupt the up/downstream gene of intergenic variants
    void *last_gene;
    void *next_gene;
};
enum func_region_type {
    _func_region_promote_to_int = -1,
    func_region_unknown = 0,
    func_region_cds,
    func_region_noncoding_intron,
    func_region_noncoding_exon,
    func_region_intron,
    func_region_utr5_exon,
    func_region_utr5_intron,
    func_region_utr3_exon,
    func_region_utr3_intron,
    func_region_intergenic,
    func_region_outrange,
};

// Sequence ontology variant types
enum mol_con {
    mc_unknown = 0,
    mc_intergenic_1KB,
    mc_upstream_1KB,
    mc_upstream_10KB,
    mc_downstream_1KB,
    mc_downstream_10KB,
    mc_intergenic,
    mc_tfbs_ablation,
    mc_tfbs_variant,
    mc_intragenic,
    mc_exon_loss,
    mc_noncoding_intron,
    mc_coding_intron,
    mc_utr3_intron,
    mc_utr5_intron,
    mc_exon_splice_sites,
    mc_splice_donor,
    mc_splice_acceptor,
    mc_donor_region,
    mc_polypyrimidine_trait,
    mc_noncoding_exon,
    mc_noncoding_splice_region,
    mc_utr5_exon,
    mc_utr3_exon,
    mc_utr5_premature_start_codon_gain,
    mc_frameshift_truncate,
    mc_frameshift_elongation,
    mc_inframe_indel,
    mc_disruption_inframe_deletion,
    mc_disruption_inframe_insertion,
    mc_conservation_inframe_deletion,
    mc_conservation_inframe_insertion,
    mc_start_loss,
    mc_start_retained,
    mc_stop_loss,
    mc_stop_retained,
    mc_stop_gained,
    mc_missense,
    mc_synonymous,
    mc_nocall,
    mc_whole_gene,
    mc_coding_variant,
};

struct intergenic_core {
    // Should be molecular consequence in intergenic region only.
    enum mol_con con1; // motifs
    enum mol_con con2; // up/downstream or inside of gene
    char *gene; // up/downstream gene name
    int gap_length; // gap between up/downstream gene and variant
};

struct mc_type {
    // start
    enum func_region_type func1;
    // end
    enum func_region_type func2;

    // for each variant, we support at max two molecular consequences, such as a missense_variant could also be 
    // splice_site_variant : missense_variant+splice_sites_variant
    // con1 used to record the main variant type, con2 used to record the supplymentary type such as splice sites and motifs
    enum mol_con con1;
    enum mol_con con2;
    // Exon or intron count. Start from 1.
    int count;
    // CDS count, for noncoding transcript always be 0.
    int count2;

    // Location of amino acids, end_amino should always be 0 for SNPs, and should be equal or greater than loc_amnio for indels
    int loc_amino;
    int loc_end_amino;
    
    // Original amino acid and mutated amino acid. Check amino type only if variants happed in cds region.
    int ori_amino; // first amino acid of original sequence
    int ori_end_amino; // last amino acid of orignal sequence
    
    // fast access for single amino acid change
    int mut_amino;

    // buffered sequences for mutated amino acids
    int n;
    int *aminos;

    // codon is break ?
    // int disrupt_codon;
    
    // Terminal sequence after loc_amino.
    int fs;
    // Extern codon, for stop-loss
    int ext;
};

struct mc_inf {
    // transcript name or Loc name in same case, name1
    char *transcript;
    // transcript version
    // int version;
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

    // the bases used to cache the realigned insertion or deletion, for repeat sequences after realignment the insert or
    // deleted bases may be changed. Need to be freed at the end.
    char *ref;
    char *alt;
    
    // offset in intron
    int offset;
    int end_offset;

    // check if repeat sequences(dup), only for insterion
    int is_dup;
    int dup_offset; // delete it
    
    // If strand is '-', convert sequence to complement strand. Default is '+'.
    char strand;

    // in frame contains stop codon
    int inframe_stop;
};

struct mc_core {
    struct mc_inf inf;
    struct mc_type type;
};

// molecular consequence description structure
struct mc {
    // point to chrom of VCF header, do NOT free it
    const char *chr;
    int start;
    int end;
    char *ref;
    char *alt;
    // variant type include SNP,INDEL only accunt for allele changes
    enum variant_type type;

    // variant description for each transcript
    int n_tran;
    struct mc_core *trans;

    // for intergenic region
    struct intergenic_core inter;
};

struct anno_mc_file {
    struct mc_handler *h;
    int n_allele; // assume n_allele == 1 
    struct mc **files;
    int n_col;
    struct anno_col *cols;
    char *tmps;
    int mtmps;
};

extern struct anno_mc_file *anno_mc_file_init(bcf_hdr_t *hdr, const char *column, const char *data, const char *rna, const char *reference, const char *name_list);
extern struct anno_mc_file *anno_mc_file_duplicate(struct anno_mc_file *f);
extern void anno_mc_file_destroy(struct anno_mc_file *f, int l);
//extern void anno_mc_core(struct anno_mc_file *f, bcf_hdr_t *hdr, bcf1_t *line);
extern int anno_mc_chunk(struct anno_mc_file *f, bcf_hdr_t *hdr, struct anno_pool *pool);


#endif
