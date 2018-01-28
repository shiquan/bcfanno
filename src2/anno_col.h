#ifndef ANNO_COL_H
#define ANNO_COL_H

#include "utils.h"
#include "htslib/vcf.h"

#define CHUNK_MAX_GAP 10000

struct anno_col;
struct anno_vcf_file;
struct anno_bed_file;
struct anno_hgvs;

#define REPLACE_MISSING  0 // replace only missing values
#define REPLACE_ALL      1 // replace both missing and existing values
#define REPLACE_EXISTING 2 // replace only if tgt is not missing
#define SET_OR_APPEND    3 // set new value if missing or non-existent, append otherwise

typedef int (*setter_vcf)(struct anno_vcf_file*, bcf_hdr_t *hdr, bcf1_t *, struct anno_col *, void *);
//typedef int (*setter_hgvs)(struct anno_hgvs_file*, bcf_hdr_t *hdr, bcf1_t *, struct anno_col *);
typedef int (*setter_bed)(struct anno_bed_file*, bcf_hdr_t *hdr, bcf1_t *, struct anno_col *);

//typedef int (*setter_pl)(bcf1_t *, struct anno_col*, void*);
typedef union {
    setter_vcf vcf;
    setter_bed bed;
    //setter_hgvs hgvs;
    //setter_pl pl;
} setter_func;


struct anno_col {
    // column index
    int icol;
    // default is REPLACE_MISSING, todo: REPLACE_EXISTSING
    int replace;
    // number : BCF_VL_*
    int number;
    // tag name
    char *hdr_key;
    // setter function
    setter_func func;
    // point to hdr names, do not free it
    const char *curr_name;
    int curr_line;
};

extern void anno_col_copy(struct anno_col *src, struct anno_col *dest);
extern int bcf_update_info_fixed(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type);

#define bcf_update_info_int32_fixed(hdr,line,key,values,n)   bcf_update_info_fixed((hdr),(line),(key),(values),(n),BCF_HT_INT)
#define bcf_update_info_float_fixed(hdr,line,key,values,n)   bcf_update_info_fixed((hdr),(line),(key),(values),(n),BCF_HT_REAL)
#define bcf_update_info_flag_fixed(hdr,line,key,string,n)    bcf_update_info_fixed((hdr),(line),(key),(string),(n),BCF_HT_FLAG)
#define bcf_update_info_string_fixed(hdr,line,key,string)    bcf_update_info_fixed((hdr),(line),(key),(string),1,BCF_HT_STR)


#endif
