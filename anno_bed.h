#ifndef ANNO_BED_HEADER
#define ANNO_BED_HEADER
#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>
#include <htslib/tbx.h>
#include "anno.h"

/* struct bed_anno_file_col { */
/*     int hdr_id; */
/*     int icol; // column number in bed file */
/*     char *a; */
/* }; */

/* struct bed_anno_line { */
/*     int rid; */
/*     int start; // 0based start */
/*     int end; // 1based end */
/*     struct bed_anno_col *cols; */
/* }; */
struct beds_anno_tsv {
    int nfields;
    int *fields;
    kstring_t string;
};

struct beds_anno_file {
    int id; // file idx
    htsFile *fp;
    tbx_t *idx;
    char *fname;
    int last_id;
    int last_start;
    int last_end;
    int n_cols;
    struct anno_col *cols;
    // memory pool    
    int cached, max;
    struct beds_anno_tsv **buffer;
    // struct bed_anno_line **buffer;
};

struct beds_options {
    int beds_is_inited;
    bcf_hdr_t *hdr_out;
    int n_files;
    int m_files;    
    struct beds_anno_file *files;
    
};
// bed format function annotation
// extern int setter_func_region(struct beds_options *opts, bcf1_t *line);

#endif
