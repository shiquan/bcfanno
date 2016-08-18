#ifndef ANNO_VCF_HEADER
#define ANNO_VCF_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>
#include <htslib/tbx.h>
#include "anno.h"
#include "vcmp.h"

struct anno_vcf_file {
    htsFile *fp;
    hts_idx_t *bcf_idx;
    hts_itr_t *itr;
    char *fname;
    int cached, max;
    bcf1_t **buffer;
    char *columns;
    int ncols;
    struct anno_col *cols;
};

struct vcfs_options {
    // hdr is header struct of input vcf file, all the annotated tags should be inited in the hdr_out before export bcf lines
    bcf_hdr_t *hdr, *hdr_out;
    // pre-indexed vcf/bcf databases
    // columns structure of vcf databases
    // struct anno_cols *cols;    
    // int i_data;
    int n_files;
    struct anno_vcf_file *files;
    // for matching annotation and VCF lines by allele
    vcmp_t *vcmp;  
    int *sample_map, nsample_map, sample_is_file;   
    int mtmpi, mtmpf, mtmps;
    int mtmpi2, mtmpf2, mtmps2;
    int mtmpi3, mtmpf3, mtmps3;
    int32_t *tmpi, *tmpi2, *tmpi3;
    float *tmpf, *tmpf2, *tmpf3;
    char *tmps, *tmps2, **tmpp, **tmpp2;
    kstring_t tmpks;
};

extern int setter_filter(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_filter(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int setter_id(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_id(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int setter_qual(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_qual(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int setter_info_flag(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_info_flag(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int setter_ARinfo_int32(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, int nals, char **als, int ntmpi);
extern int setter_info_int(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_info_int(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int setter_ARinfo_real(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, int nals, char **als, int ntmpf);
extern int setter_info_real(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_info_real(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst);
extern int setter_ARinfo_string(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, int nals, char **als);
extern int setter_info_str(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_info_str(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_format_gt(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int count_vals(struct anno_line *tab, int icol_beg, int icol_end);
extern int setter_format_int(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int setter_format_real(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int setter_format_str(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_format_int(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_format_real(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_format_str(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data);

#endif
