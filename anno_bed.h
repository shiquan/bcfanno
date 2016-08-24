#ifndef ANNO_BED_HEADER
#define ANNO_BED_HEADER
#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>
#include <htslib/tbx.h>
#include "anno.h"

struct anno_bed_col {
    int hdr_id;
    int icol; // column number in bed file
    char *a;
};

struct anno_bed_line {
    int rid;
    int start; // 0based start
    int end; // 1based end
    struct anno_bed_col *cols;
};

struct anno_bed_file {
    htsFile *fp;
    tbx_t *idx;
    int n_cols;
    struct anno_col *cols;
    // memory pool
    int cached, max;
    struct anno_bed_line **buffer;
};

struct beds_options {
    
    int n_files;
    struct anno_bed_file *files;
    
};
// bed format function annotation
extern int setter_func_region(struct beds_options *opts, bcf1_t *line);

#endif
