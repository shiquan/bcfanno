#ifndef NAME_LIST_H
#define NAME_LIST_H
extern void *name_hash_init(const char *fname);
extern int name_hash_key_exists(void *hash, char *key);
extern int name_hash_key_add(void *hash, char *key);
extern int name_hash_key_delete(void *hash, char *key);
extern void name_hash_destroy(void *hash);
#endif
