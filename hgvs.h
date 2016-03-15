/* hgvs
 *
 */
#ifndef VCFANNO_HGVS_HEADER
#define VCFANNO_HGVS_HEADER

#include "anno.h"
#include <htslib/khash.h>

enum variant_type {
    is_ref,
    is_del,
    is_ins,
    is_indel,
    is_snp,
    is_cnv, // cnv is also a kind of SV
    is_sv,    
};

enum func_type {
    is_utr5, // ncRNA should not have UTR!
    is_utr3,
    is_cds,
    is_splitsite,
    is_promotor,
    is_intron,
    is_unknown,
};

struct trans {
    char *hgvs_trans_name;
    char *alt_seq; // alternative allele sequences
    int exon_count; // init with -1
    int cds_count; // init with -1
    int cds_pos;
    int cds_shift; // init with 0, no shift
};

struct hgvs_name {
    char *hgvs_gene_name;
    int geno_pos_start;
    int ref_length;
    enum variant_type vtype;
    int n_trans; // there are several transcipts for one gene
    struct trans * trans;
};


#define clear_all(x) (x ＝ 0)
#define clear_flag(x, n) (x &= ~(n))

#define REFGENE_PRASE_NONA   0
#define REFGENE_PRASE_BIN    1
#define REFGENE_PRASE_NAME1  2
#define REFGENE_PRASE_REG    4
#define REFGENE_PRASE_STRAND 8
#define REFGENE_PRASE_EXONS  (1<<4 | 8)
#define REFGENE_PRASE_ALL (REFGENE_PRASE_BIN | REFGENE_PRASE_NAME1 | REFGENE_PRASE_REG| REFGENE_PRASE_EXONS)

struct refgene_entry {
    int prase_flag;
    char *buffers;
    int n_buffer;
    int m_buffer;
    int bin; // see USCS bin scheme for details, be used to check regions
    
    const char *name1; // rna name usually
    const char *name2; // gene name

    int rid; // chromosome id, should be one of contigs in the VCF header, -1
    
    int strand; // strand_plus, strand_minus, strand_unknown
    uint32_t exon_begin; // exon region begin from xxx
    uint32_t exon_end; // exon end to xxx, remember exon contains cds and utr
    uint32_t cds_start; // cds region start position, 5'UTR is exon_begin to cds_start-1.
    uint32_t cds_stop; // 3'UTR is cds_end+1 to exon_end in plus strand
    int total_exons; // exon number
    struct region *exons;
    struct region *cds;
};

struct anno_hgvs_option {
    struct refgene_file refgene;
    struct refrna_file refrna;
    void *name_hash1; // transcript hash
    void *name_hash2; // gene hash
    struct annot_col *cols;
    int ncols;
    struct refgene_entry *buffers;
    int mbuffer, nbuffer;
};

/* retrieve data from local tabix indexed refgene database, see manual for more details
   about databases */
extern void extract_refgene(struct refgene_entry *entry, int type);
extern int local_setter_hgvs_func(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int local_setter_hgvs_names(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);

extern void refgene_fill_buffers_from_local(struct anno_hgvs_option *hgvs_opts, int prase_flag);

/* genes list and transcript list would be checked in this step, empty file will be warned,
 *
 * Note: 
 * 1. only random one transcript would be listed if no transcript and gene list specified.
 * 2. all transcripts would be listed for the gene specified in the gene list.
 * 3. if transcipt list is specified, only these transcipts would be annotated.
 *
 */
extern int init_refgene_gene_trans_list(const char *gene_list, const char *nm_list);

//extern bcf_t * insert_hgvs_tag(bcf_t *line);


#endif
