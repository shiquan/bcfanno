#ifndef SEQIDX_H
#define SEQIDX_H

#include "utils.h"
#include "htslib/faidx.h"

struct seqidx {
    const char *file;
    faidx_t *idx;
};

#endif 
