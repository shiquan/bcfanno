#include "utils.h"
#include "wrap_pileup.h"
#include <limits.h>
struct plp_ref *plp_ref_init()
{
    struct plp_ref *r = malloc(sizeof(struct plp_ref));
    r->ref_id[0] = -1;
    r->ref_id[1] = -1;
    r->ref[0] = NULL;
    r->ref[1] = NULL;
    r->ref_len[0] = 0;
    r->ref_len[1] = 1;
    return r;
}

//
// Assume input bam or VCF records are sorted by chromosome and genomic positions,
// so we just need load the reference sequence one time, and cached in the memory.
// This function was designed to check the memory and retrieve the cached sequence
// or update the memory.
//
// return 1 on success, 0 on failure.
//
// This function is copied and edited from bam_plcmd.c, all credit of this code
// goes to Li Heng's samtools.
//
int plp_get_ref(struct plp_ref *r, char *seqname, faidx_t *fai, int tid, char **ref, int *ref_len)
{
    if ( !fai || !r ) {
        ref = NULL;
        *ref_len = 0;
        return 0;
    }

    if ( tid == r->ref_id[0] ) {
        *ref = r->ref[0];
        *ref_len = r->ref_len[0];
        return 1;
    }

    if ( tid == r->ref_id[1] ) {
        // Last, swap over
        int tmp;
        tmp = r->ref_id[0];  r->ref_id[0]  = r->ref_id[1];  r->ref_id[1]  = tmp;
        tmp = r->ref_len[0]; r->ref_len[0] = r->ref_len[1]; r->ref_len[1] = tmp;

        char *tc;
        tc = r->ref[0]; r->ref[0] = r->ref[1]; r->ref[1] = tc;
        *ref = r->ref[0];
        *ref_len = r->ref_len[0];
        return 1;
    }

    // New, so migrate to old and load new
    if ( r->ref[1] != NULL ) free(r->ref[1]);
    
    r->ref[1]     = r->ref[0];
    r->ref_id[1]  = r->ref_id[0];
    r->ref_len[1] = r->ref_len[0];

    r->ref_id[0] = tid;
    r->ref[0] = faidx_fetch_seq(fai,
                                seqname,
                                0,
                                INT_MAX,
                                &r->ref_len[0]);
    if (!r->ref[0]) {
        r->ref[0] = NULL;
        r->ref_id[0] = -1;
        r->ref_len[0] = 0;
        *ref = NULL;
        return 0;
    }

    *ref = r->ref[0];
    *ref_len = r->ref_len[0];
    return 1;
}

