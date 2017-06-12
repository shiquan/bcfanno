#ifndef NUMBER_HEADER
#define NUMBER_HEADER
#include <stdio.h>
#include <stdlib.h>

extern int get_numbase(const char *s);
extern int get_numbase_l(const char *s, int l);
extern int is_ieee_magic_val(const char *val);
extern double nondec2num(char *str, int length);
extern int check_num_likely(const char *str);
extern int check_num_likely_l(const char *str, int length);
extern double force2num(char *str);
extern double force2num_l(char *str, int l);
extern int str2int(char *str);
extern int str2int_l(char *str, int l);





#endif
