#ifndef VCFANNO_CONFIG_H
#define VCFANNO_CONFIG_H

#include <stdio.h>
#include <stdlib.h>

struct refgene_config {
    // if refgene_is_set == 1, genepred_fname and columns are mandatory.
    int refgene_is_set;
    const char *genepred_fname;
    const char *refseq_fname;
    const char *trans_list_fname;
    const char *gene_list_fname;
    const char *columns;
};

struct file_config {
    // file path
    const char *fname;
    // columns string
    const char *columns;
};
struct vcfs_config {
    // vcf files number
    int n_vcfs;
    struct file_config *files;
};
struct beds_config {
    // bed files number
    int n_beds;
    struct file_config *files;
};

// skip other keys except author, config_id and reference_version
struct vcfanno_config {
    const char *author;
    const char *config_id;
    const char *reference_version;
    struct vcfs_config vcfs;
    struct beds_config beds;
    struct refgene_config refgene;
};

extern struct vcfanno_config *vcfanno_config_init(void);
extern struct vcfanno_config_destroy(struct vcfanno_config *);
extern int vcfanno_load_config(struct vcfanno_config *, const char *);


#define ANNOCONFIG_INIT { NULL, NULL, NULL, 0, 0, NULL }

enum anno_type { anno_is_vcf, anno_is_sql, anno_is_unknown };

//enum api_error_num { api_is_unreach, api_fail_open, api_fail_index };

struct summary {
    const char *name;
    const char *version;
    const char *author;
    const char *ref_version;
    char *path; // ip address or environment path
};

struct vcf_sql_api {
    enum anno_type type;
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
    //int intron_edge;
};

struct configs {
    const char *path_string;
    struct summary *summary;
    struct anno_data_file *anno;
    int n_theads;
    int n_apis;
    struct vcf_sql_api *apis;
};

extern struct configs anno_config_file;

extern int has_hgvs;

extern void config_release();

extern int load_config(const char *json);



#endif
