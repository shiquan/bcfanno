#ifndef TABLE_HASH_HEADER
#define TABLE_HASH_HEADER
#include <stdio.h>
#include <stdlib.h>
extern void table_set_flag();
extern void table_clear_flag();
extern int table_check_flag();
extern int table_read(const char *fname);
extern char *table_convert_name(char *name);
extern int table_release();

#endif
