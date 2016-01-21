#ifndef VCFANNO_CONFIG_H
#define VCFANNO_CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include "kson.h"
//#include "anno_setter.h"

#define DEFAULT_INTRON_EDGE 30

enum api_type { api_is_vcf, api_is_sql, api_is_unknown };

//enum api_error_num { api_is_unreach, api_fail_open, api_fail_index };

struct summary {
    const char *  name;
    const char *  version;
    const char *  author;
    const char *  ref_version;
    char *path; // ip address or environment path
};

struct vcf_sql_api {
    enum api_type type;
    const char *vfile;
    const char *dynlib;
    const char *columns;
    int has_summary:1;
    struct summary *summary;
    //    enum api_error_num;
};

extern struct configs anno_config_file;
extern int has_hgvs;

extern void config_release();

#define safe_free(x) do                                                 \
    {                                                                   \
        void **_px = (void**)&(x);                                      \
        if (_px==NULL || *_px == NULL) error("[error] double free %s, %d", __FUNCTION__, __LINE__); \
        free(*_px);                                                     \
        *_px = NULL;                                                    \
    } while(0)

extern int load_config(const char *json);



#endif
