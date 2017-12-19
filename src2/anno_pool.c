// anno_pool.c - A pool of bcf structure
#include "anno_pool.h"
#include "utils.h"

static struct anno_pool *anno_pool_init(int m)
{    
    struct anno_pool *p = malloc(sizeof(*p));
    p->n_reader = 0;
    p->m = m;
    p->readers = malloc(m * sizeof(bcf1_t*));
    return p;
}

struct anno_pool *anno_reader(htsFile *fp, bcf_hdr_t *hdr, int n_record) {
    
    struct anno_pool *p = anno_pool_init(n_record);
    for ( ;; ) {
        p->readers[p->n_reader] = bcf_init();
        if ( bcf_read(fp, hdr, p->readers[p->n_reader]) )
            break;
        p->n_reader++;
        if ( p->n_reader == p->m )
            break;        
    }
    return p;
}

