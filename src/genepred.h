#ifndef GENEPRED_HEADER
#define GENEPRED_HEADER
#include "htslib/hts.h"
#include "htslib/tbx.h"
#include "htslib/faidx.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include "sort_list.h"

#define BLOCK_START 0
#define BLOCK_END   1

#define read_start(a, l) (a[BLOCK_START][l])
#define read_end(a, l) (a[BLOCK_END][l])
#define get_length_exon(a, l) (a[BLOCK_END][l] - a[BLOCK_START][l] + 1)

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

    // realign transcript sequence with reference sequence, record in CIGAR format.
    // M for match, I for insertion, D for deletion
    // for an example, 256M1D43M15I stand for 256 base matched and 1 deletion in transcript, and continued with 43 matched bases,
    // and 15 insertion in the tail (usually the tailAs).
    char *realn; // orignal string, should be free and point to NULL after parse.
    int n_cigar;
    int *cigars;
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
    int realn;
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
extern int genepred_read_line(struct genepred_spec *spec, struct genepred_line *line);
struct genepred_spec *genepred_spec_init();
void genepred_spec_destroy(struct genepred_spec *spec);
extern int parse_line(kstring_t *string, struct genepred_line *line);
extern int parse_line_locs(struct genepred_line *line);
extern void genepred2line(struct genepred_line *line, kstring_t *str);

extern int genepred_load_data(struct genepred_spec *spec, const char *fname);
extern int genepred_load_fasta(struct genepred_spec *spec, const char *fname);
extern int genepred_load_genes(struct genepred_spec *spec, const char *fname);
extern int genepred_load_trans(struct genepred_spec *spec, const char *fname);
extern struct genepred_line *genepred_retrieve_gene(struct genepred_spec *spec, const char *name);
extern struct genepred_line *genepred_retrieve_trans(struct genepred_spec *spec, const char *name);
extern struct genepred_line *genepred_retrieve_region(struct genepred_spec *spec, char *name, int start, int end);
extern void generate_dbref_database(struct genepred_line *line);
#endif
