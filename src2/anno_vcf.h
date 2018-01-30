#ifndef ANNO_VCF_H
#define ANNO_VCF_H

#include "utils.h"
#include "vcmp.h"
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "anno_pool.h"

struct anno_vcf_buffer {
    int no_such_chrom;
    int last_rid;
    int cached, max;
    int i_chunk;
    bcf1_t **buffer;
    vcmp_t *vcmp;
    int mtmpi, mtmpf, mtmps;
    int mtmpi2, mtmpf2, mtmps2;
    int mtmpi3, mtmpf3, mtmps3;
    int32_t *tmpi, *tmpi2, *tmpi3;
    float *tmpf, *tmpf2, *tmpf3;
    char *tmps, *tmps2, **tmpp, **tmpp2;
    kstring_t tmpks;
};

struct anno_vcf_file {
    const char *fname;
    htsFile *fp;
    // header of vcf
    bcf_hdr_t *hdr;
    hts_idx_t *bcf_idx;
    tbx_t     *tbx_idx;
    // iterator
    hts_itr_t *itr;

    int n_col;
    struct anno_col *cols;
    struct anno_vcf_buffer *buffer;
};

extern struct anno_vcf_file *anno_vcf_file_init(bcf_hdr_t *hdr, const char *fname, char *column);
extern void anno_vcf_file_destroy(struct anno_vcf_file *f);
extern struct anno_vcf_file *anno_vcf_file_duplicate(struct anno_vcf_file *f);
extern void anno_vcf_file_destroy(struct anno_vcf_file *f);
extern int anno_vcf_core(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line);
extern int anno_vcf_chunk(struct anno_vcf_file *f, bcf_hdr_t *hdr, struct anno_pool *pool);

// APIs from vcf_annos.c
extern int vcf_setter_filter(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_id(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_qual(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_info_flag(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_info_int(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_info_real(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, void *data);
extern int vcf_setter_info_str(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, void *data);

    
#endif
