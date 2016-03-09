#ifndef VCFANNO_PLUGIN_HEADER
#define VCFANNO_PLUGIN_HEADER

#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/khash_str2int.h>
#include <dlfcn.h>

/*
 * try to convert any data from SQL into annot_line struct and annotated by core_annotate()
 */
typedef struct annot_line {
    char **cols;
    int ncols, mcols;
    char **als;
    kstring_t line;
    //int rid, start, end;
} annot_line_t;

/** 
 * Plugin API:
 * ----------
 *   const char *about(void)
 *       - short description used by vcfanno --list config.json
 *
 *   int test(void)
 *      - check if host is reachable
 *
 *   int init(int argc, char **argv, bcf_hdr_t *in_hdr, bcf_hdr_t *out_hdr)
 *      - called once at startup, allows to initialize local variables.
 *      Return 1 to suppress normal VCF/BCF header output, -1 on critical
 *      errors, 0 otherwise.
 *
 *   annot_line_t *process(const char *keys)
 *      - called for each online record, return NULL for no output
 *
 *   void destroy(void)
 *      - called after all lines have been processed to clean up 
 */

typedef int (*dl_init_func)(int, char**, bcf_hdr_t *, bcf_hdr_t *);
typedef char * (*dl_about_func)(void);
typedef int (*dl_test_func)(void);
typedef annot_line_t * (*dl_process_func)(annot_line_t *anno);
typedef void (*dl_destroy_func)(void);

struct plugin_functions {
    dl_init_func init;
    dl_about_func about;
    dl_process_func process;
    dl_destroy_func destroy;
};

/* struct plugin_handler { */
/*     struct plugin_functions plugin; */
/*     int nplugin_paths; */
/*     char **plugin_paths; */
/* }; */

extern void *dlopen_plugin(const char *fname);



#endif
