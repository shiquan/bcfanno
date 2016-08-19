#ifndef ANNO_BED_HEADER
#define ANNO_BED_HEADER
#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>
#include <htslib/tbx.h>
#include "anno.h"

struct anno_bed_line {
    int rid;
    int start; // 0based start
    int end; // 1based end
    char *a;
};

struct anno_bed_file {
    htsFile *fp;
    tbx_t *idx;
    int rid;
    uint32_t start;
    uint32_t end;
    struct anno_col *col;
    // memory pool
    int cached, max;
    struct anno_bed_line *a;
};

struct beds_options {
    char **fnames;
    int n_files;
    struct anno_bed_files *files;
    
};
// bed format function annotation
extern int setter_func_region(struct beds_options *opts, bcf1_t *line);

#endif
