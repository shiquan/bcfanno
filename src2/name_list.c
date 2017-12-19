#include "utils.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/hts.h"

KHASH_MAP_INIT_STR(name_hash, char*)

//typedef khash_t(name) hash;

void *name_hash_init(const char *fname)
{
    int i, n = 0;
    khiter_t k;
    khash_t(name_hash) *hash;
    int ret;
    char **names = hts_readlist(fname, 1, &n);
    if ( n == 0 )
        return NULL;

    hash = kh_init(name_hash);
    for ( i = 0; i < n; ++i ) k = kh_put(name_hash, hash, names[i], &ret);
    free(names);
    return hash;
}
// 1 on exist
// 0 on empty
int name_hash_key_exists(void *hash, char *key)
{
    khash_t(name_hash) *h = (khash_t(name_hash)*)hash;
    khiter_t k;
    k = kh_get(name_hash, h, key);
    if ( k == kh_end(h))
        return 0;
    return 1;
}
int name_hash_key_delete(void *hash, char *key)
{
    khash_t(name_hash) *h = (khash_t(name_hash)*)hash;
    khiter_t k;
    k = kh_get(name_hash, h, key);
    if ( k == kh_end(h))
        return 0;
    kh_del(name_hash, h, k);        
    return 1;
}
// 1 on success add
// 0 on dup
int name_hash_key_add(void *hash, char *key)
{
    khash_t(name_hash) *h = (khash_t(name_hash)*)hash;
    khiter_t k;
    int ret;
    k = kh_get(name_hash, h, key);
    if ( k == kh_end(h)) {
        kh_put(name_hash, h, key, &ret);        
        return 1;
    }
    return 0;
}
void name_hash_destroy(void *hash)
{
    khash_t(name_hash) *h = (khash_t(name_hash)*)hash;
    khiter_t k;
    for ( k = 0; k < kh_end(h); ++k )
        if ( kh_exist(h, k) )
            free((char*)kh_key(h, k));
    kh_destroy(name_hash, h);    
}
