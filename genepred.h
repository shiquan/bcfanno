#ifndef GENEPRED_HEADER
#define GENEPRED_HEADER
#include "htslib/hts.h"
#include "htslib/tbx.h"
#include "htslib/faidx.h"
#include "htslib/bgzf.h"
#include "sort_list.h"

// 
// The dna reference offset for DNA coding and noncoding transcript consist of two parts. the offset value
// and the function region construct a 32-bits value.
// coding/nocoding coordinate
// 32                        4321
// |_________________________||||
// ONLY one-bit below is accept per value
// 
// offset 1:  is_noncoding
// offset 2:  is_coding
// offset 3:  is_utr5
// offset 4:  is_utr3
// 

/* #define REG_NONCODING  1 */
/* #define REG_CODING     2 */
/* #define REG_UTR5       4 */
/* #define REG_UTR3       8 */
/* #define REG_MASK       0xF */
/* #define TYPEBITS       4 */

/* #define read_type(a) ((a) & REG_MASK) */
/* #define read_loc(a)  ((a) >> TYPEBITS) */
/* #define compact_loc(a, t) ((a)<<TYPEBITS | ((t) & REG_MASK)) */

#define BLOCK_START 0
#define BLOCK_END   1

#define read_start(a, l) (a[BLOCK_START][l])
#define read_end(a, l) (a[BLOCK_END][l])

struct genepred_line {
    struct genepred_line *next;
    // Chromosome id.
    char *chrom;
    int txstart;
    int txend;
    char strand;
    // Usually the transcript name.
    char *name1;
    // Gene name.
    char *name2;

    int cdsstart;
    int cdsend;

    // Forward length and backward length are the UTRs in the transcript, for different
    // strand the forward and backward could be UTR5' or UTR3'; and most importantly,
    // for noncoding transcript, no UTRs, so the forward and backward length should
    // always be 0.
    int utr5_length;
    int utr3_length;

    // Length of this transcript.
    int reference_length;

    int exon_count;
    // [start, end]
    int *exons[2];
    int loc_parsed;
    // Offsets location on gene.
    // int *dna_ref_offsets[2];
    // Transcript locations of each block, coding transcripts consist of UTRs and CDS,
    // so loc[0,1] is not the edge of coding sequences. Consistant with the strand.
    int *loc[2];    
};

struct genepred_format {
    int chrom;
    // Usually name1 is transcript, name2 is gene.
    int name1;
    int name2;
    int strand;
    int txstart;
    int txend;
    int cdsstart;
    int cdsend;
    int exon_count;
    int exonstarts;
    int exonends;
};

struct list {
    void *hash;
    int lines;
    char **reads;
};
struct genepred_spec {
    const char *data_fname;
    htsFile *fp;
    const char *fai_fname;
    faidx_t *fai;
    tbx_t *idx;
    struct list *genes;
    struct list *trans;
};

extern void set_format_refgene();
extern void set_format_genepred();
extern void set_format_refflat();

extern struct genepred_line *genepred_line_create();
extern void genepred_line_destroy(void *line);

struct genepred_spec *genepred_spec_init();
void genepred_spec_destroy(struct genepred_spec *spec);
extern int parse_line(kstring_t *string, struct genepred_line *line);
extern int parse_line_locs(struct genepred_line *line);

extern int genepred_load_data(struct genepred_spec *spec, const char *fname);
extern int genepred_load_fasta(struct genepred_spec *spec, const char *fname);
extern int genepred_load_genes(struct genepred_spec *spec, const char *fname);
extern int genepred_load_trans(struct genepred_spec *spec, const char *fname);
extern struct genepred_line *genepred_retrieve_gene(struct genepred_spec *spec, const char *name);
extern struct genepred_line *genepred_retrieve_trans(struct genepred_spec *spec, const char *name);
extern struct genepred_line *genepred_retrieve_region(struct genepred_spec *spec, char *name, int start, int end);
extern void generate_dbref_database(struct genepred_line *line);
#endif
