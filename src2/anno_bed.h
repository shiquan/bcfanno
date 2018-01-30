#ifndef ANNO_BED_H
#define ANNO_BED_H
#include "utils.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/kstring.h"
#include "htslib/tbx.h"
#include "anno_pool.h"

struct anno_bed_tsv {
    int  n_field;
    int *fields;
    int  start;
    int  end;
    kstring_t string;
};
struct anno_bed_buffer {
    int no_such_chrom;
    int last_rid;
    int last_start;
    int last_end;
    int cached;
    int i;
    int end_pos_for_skip;
    int max;
    struct anno_bed_tsv **buffer;
    int mtmps;
    char *tmps;
};
struct anno_bed_file {
    //int id;
    const char *fname;
    htsFile *fp;
    tbx_t *idx;
    // set to 0 if records are NOT overlapped, records will be refreshed only if out of range
    int overlapped;
    int n_col;
    struct anno_col *cols;
    struct anno_bed_buffer *buffer;
};

extern int anno_bed_core(struct anno_bed_file *file, bcf_hdr_t *hdr, bcf1_t *line);
extern struct anno_bed_file *anno_bed_file_init(bcf_hdr_t *hdr, const char *fname, char *column);
extern struct anno_bed_file *anno_bed_file_duplicate(struct anno_bed_file *f);
extern void anno_bed_file_destroy(struct anno_bed_file *f);
extern int anno_bed_chunk(struct anno_bed_file *file, bcf_hdr_t *hdr, struct anno_pool *pool );

#endif
