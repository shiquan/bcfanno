// anno_pool.c - A pool of bcf structure
#include "anno_pool.h"
#include "utils.h"

static struct anno_pool *anno_pool_init(int m)
{    
    struct anno_pool *p = malloc(sizeof(*p));
    memset(p, 0, sizeof(*p));
    p->m = m;
    p->readers = malloc(m * sizeof(bcf1_t*));
    return p;
}

void update_chunk_region(struct anno_pool *pool)
{
    int i_chunk;
    int start = -1, end = -1, last_rid = -1, last_pos = -1;
    pool->curr_line = pool->readers[pool->n_chunk];
    start = pool->curr_line->pos;
    last_rid = pool->curr_line->rid;
    
    for ( i_chunk = pool->n_chunk; i_chunk < pool->n_reader; ++i_chunk ) {
        int pos = pool->readers[i_chunk]->pos;
        int rid = pool->readers[i_chunk]->rid;
        if ( last_rid == -1 ) last_rid = rid;
        if ( last_rid != rid ) break;
        if ( start == -1 ) start = pos;        
        if ( last_pos != -1 && pos - last_pos > CHUNK_MAX_GAP ) break;
        last_pos = pos;
        end = pos;
    }
    // 
    pool->curr_start = start;
    pool->curr_end   = end;

    pool->i_chunk = pool->n_chunk;
    pool->n_chunk = i_chunk;
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

