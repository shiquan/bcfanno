/* hgvs
 *
 */
#ifndef VCFANNO_HGVS_HEADER
#define VCFANNO_HGVS_HEADER

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

struct refgene {
    struct refgene *next;
    int is_root;
    int chr_id;
    char *gene_name;
    char *trans_name;
    int exon_count;
    int strand;
    int exon_start;
    int exon_end;
    int cds_start;
    int cds_end;
    int *exon_starts;
    int *exon_stops;
};

/* INPUT should looks like this
 *  chromo,pos,ref,alt
 * Database should contain trans fasta, regions
 * - fasta should be indexed by samtools faidx
 * - refgene should be sorted and indexed by tabix
 */

/* retrieve data from local tabix indexed refgene database, see manual for more details
   about databases */
extern struct refgene * retrieve_refgene_from_local(const char *fname, int retrieve_rule);

extern void release_refgene_list(struct refgene *root);

/* genes list and transcript list would be checked in this step, empty file will be warned,
 *
 * Note: 
 * 1. only random one transcript would be listed if no transcript and gene list specified.
 * 2. all transcripts would be listed for the gene specified in the gene list.
 * 3. if transcipt list is specified, only these transcipts would be annotated.
 *
 */
extern int init_refgene_gene_trans_list(const char *gene_list, const char *nm_list);

extern bcf_t * insert_hgvs_tag(bcf_t *line);


#endif
