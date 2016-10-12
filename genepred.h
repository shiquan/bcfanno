#ifndef GENEPRED_HEADER
#define GENEPRED_HEADER
#include <stdio.h>
#include <stdlib.h>
#include "htslib/kstring.h"
#include "htslib/tbx.h"
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

#define REG_NONCODING  1
#define REG_CODING     2
#define REG_UTR5       4
#define REG_UTR3       8
#define REG_MASK       0xF

#define TYPEBITS   4

#define read_type(a) ((a) & REG_MASK)
#define read_loc(a)  ((a) >> TYPEBITS)
#define compact_loc(a, t) ((a)<<TYPEBITS | ((t) & REG_MASK))

#define BLOCK_START 0
#define BLOCK_END   1

#define read_start(a, l) (a[BLOCK_START][l])
#define read_end(a, l) (a[BLOCK_END][l])

struct genepred {
    // mark the memory is allocated or not inside this struct,
    // clear == 0 for empty, clear == 1 recall clear_genepred to free it.
    int clear;
    // chromosome names, usually names is start with "chr", like chr1, chr2, ..., chrM, which is UCSC style.
    // however, for other insitutions they have different preference, Ensemble like to use ensemble id, NCBI like
    // use RefSeq accession to represent the chromosome and updated version. Please make sure all the databases have
    // same nomenclature before annotation.
    char *chrom;    
    uint32_t txstart;
    uint32_t txend;
    // '-' for minus, '+' for plus strand
    char strand;
    // usually transcript name
    char *name1;
    // gene name
    char *name2; 
    uint32_t cdsstart;
    uint32_t cdsend;
    // forward length and backward length are the UTRs in the transcript, for differnt strand forward or backward
    // could be UTR5 or UTR3; for noncoding transcript, no UTRs, forward and backward length should always be 0
    uint32_t forward_length;
    uint32_t backward_length;
    // length of this transcript
    uint32_t reference_length;
    uint32_t exoncount;
    // [start, end]
    uint32_t *exons[2];
    uint32_t *dna_ref_offsets[2];
    // transcript locations of each block, coding transcript consist of UTRs and CDS, so loc[] is not the
    // edge of coding sequences
    uint32_t *loc[2];
};

// Sometime transcript may span a huge region on genome, and parse genepred record may cost a lot of CPU,
// so for better performance, put the genepred records in the memory pool.
struct gp_mempool {
    int rid;
    uint32_t start;
    uint32_t end;
    // l for used length, m for max length, i for inited length
    int l, m, i;
    struct genepred *a;
};

#define GENEPRED_MEMORY_INIT {.rid = -1, .start = 0, .end = 0, .l = 0, .m = 0, .i = 0, .a = NULL }

// map each column right for each format.
// in default refgene have one more `bin` column in the first column than genepred format
struct genepred_format {
    int chrom;
    // usually name1 is transcript, name2 is gene
    int name1; 
    int name2; 
    int strand;
    int txstart;
    int txend;
    int cdsstart;
    int cdsend;
    int exoncount;
    int exonstarts;
    int exonends;
};

// set the format to parse genepred file
extern void set_format_refgene();
extern void set_format_genepred();
extern void set_format_refflat();

extern void genepred_parser(kstring_t *string, struct genepred *line);
extern char *gp_describe(struct genepred *line);
extern void gp_push_mempool(struct gp_mempool *pool, kstring_t *str);
extern void gp_update_mempool(struct gp_mempool *pool);
extern void gp_clear_mempool(struct gp_mempool *pool);
extern tbx_t *load_genepred_file(const char *fn);
extern void generate_dbref_database(struct genepred *line);
#endif
