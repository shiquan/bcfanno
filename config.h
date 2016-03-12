#ifndef VCFANNO_CONFIG_H
#define VCFANNO_CONFIG_H

#include <stdio.h>
#include <stdlib.h>
//#include "anno_setter.h"

#define DEFAULT_INTRON_EDGE 30

#define ANNOCONFIG_INIT { NULL, NULL, NULL, 0, 0, NULL }

enum api_type { api_is_vcf, api_is_sql, api_is_unknown };

//enum api_error_num { api_is_unreach, api_fail_open, api_fail_index };

struct summary {
    const char *name;
    const char *version;
    const char *author;
    const char *ref_version;
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

struct anno_data_file {
    char *refgene_file_path;
    char *transcripts_list;
    char *genes_list;
    char *columns;
    int intron_edge;
};

struct configs {
    const char *path_string;
    struct summary *summary;
    struct anno_data_file *anno;
    int n_theads;
    int n_apis;
    struct vcf_sql_api * apis;
};

extern struct configs anno_config_file;

extern int has_hgvs;

extern void config_release();

extern int load_config(const char *json);



#endif
