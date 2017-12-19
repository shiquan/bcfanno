#ifndef ANNO_POOL_H
#define ANNO_POOL_H

#include <stdlib.h>
#include "htslib/vcf.h"
#include "htslib/hts.h"

struct anno_pool {
    // maximum number of records
    int m;
    int n_reader;
    bcf1_t **readers;
    void *arg;
};

extern struct anno_pool *anno_reader(htsFile *fp, bcf_hdr_t *hdr, int n_record);

#endif
