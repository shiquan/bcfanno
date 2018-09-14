#ifndef WRAP_PILEUP
#define WRAP_PILEUP

#include <stdlib.h>
//#include "htslib/sam.h"
#include "htslib/faidx.h"

struct plp_ref {
    char *ref[2];
    int ref_id[2];
    int ref_len[2];
};

extern struct plp_ref *plp_ref_init();
extern int plp_get_ref(struct plp_ref *r, char *seqname, faidx_t *fai, int tid, char **ref, int *ref_len);

#endif
