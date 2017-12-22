#ifndef ANNO_HGVS_H
#define ANNO_HGVS_H
#include "utils.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/tbx.h"

struct anno_hgvs {
    htsFile *fp;
    faidx_t *rna_fai;
    tbx_t   *idx;
    void    *gene_hash;
    void    *trans_hash;
    int      n_col;
    struct anno_col *cols;
};

struct anno_hgvs_file {
    struct hgvs_handler *h;
    int n_allele; // assume n_allele == 1 
    struct hgvs **files;
    int n_col;
    struct anno_col *cols;
    char *tmps;
    int mtmps;
};

extern struct anno_hgvs_file *anno_hgvs_file_init(bcf_hdr_t *hdr, const char *column, const char *data, const char *rna, const char *reference);
extern struct anno_hgvs_file *anno_hgvs_file_duplicate(struct anno_hgvs_file *f);
extern void anno_hgvs_file_destroy(struct anno_hgvs_file *f);
extern void anno_hgvs_core(struct anno_hgvs_file *f, bcf_hdr_t *hdr, bcf1_t *line);
#endif
