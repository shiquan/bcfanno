#ifndef VCFANNO_HEADER
#define VCFANNO_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include <dlfcn.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>
#include <htslib/tbx.h>
#include <htslib/faidx.h>
#include "plugin.h"
#include "config.h"

#define KSTRING_INIT {0, 0, 0}

// vcfanno is designed to annotate vcf files by retrieving tags from different databases. For now, vcfanno only
// accept three different format of databases.
// for genotype databases, BCF/VCF format is required and should be preindexed by bcftools,
// for function region databases, bed format is required and should be packed by bgzip and indexed by tabix,
// for hgvs databases, only used for generating HGVS names, here require genepred format, see our manual for more details.

// handler of function regions data in bed format, goto anno_bed.h for details
struct beds_options;
// handler of variantions data in vcf/bcf format, goto anno_vcf.h for details
struct vcfs_options;
// handler of refgene data in genepred format, goto anno_hgvs.h for details
struct refgene_options;

// anno_line also defined in plugin.h for compile dynamic library
struct anno_line {
    char **cols;
    int ncols, mcols;
    int nals, mals;
    char **als;
    kstring_t line;
    int rid, start, end;
};

struct anno_col;

typedef int (*setter_vcf)(struct vcfs_options *, bcf1_t *, struct anno_col *, void *);
typedef int (*setter_hgvs)(struct refgene_options *, bcf1_t *, struct anno_col *);
typedef int (*setter_bed)(struct beds_options *, bcf1_t *, struct anno_col *);

typedef union {
    setter_vcf vcf;
    setter_bed bed;
    setter_hgvs hgvs;
} setter_func;

struct anno_col {
    int ifile; // file idx
    int icol; // col idx for tab seperated file, for BCF and VCF file icol is unused
    int replace; // default is REPLACE_MISSING, for REPLACE_EXISTSING there are some bugs in htslib to support
    int number;  // number: one of BCF_VL_* types
    char *hdr_key;
    setter_func setter;
    //int (*setter)(struct vcfs_options *, bcf1_t *, struct anno_col *, void*);
};

typedef void (* rel_func)(void*);

extern void safe_release(void *p, rel_func func);

// regions struct for update pre-index cache
struct region1 {
    // 0 based.
    uint32_t start;
    // 1 based. 0 for init
    uint32_t end;
};

struct region {
    int n, m, c;
    struct region1 *regs;
};

// extern struct anno_col *init_columns(const char *rules, bcf_hdr_t *in, bcf_hdr_t *out, int *n, enum anno_type type);

struct sql_connect {
    void *connect;
    // index id
    int idx; 
    char *column;
    int ncols;
    struct anno_col *cols;
    int nbuffers, mbuffer;
    struct anno_line *buffers;
    // in memory cached regions
    struct region *regs; 
    // current position: chr name index to names[]
    int iseq; 
    int start, end;
    int prev_seq, prev_start;
};

#define REPLACE_MISSING  0 // replace only missing values
#define REPLACE_ALL      1 // replace both missing and existing values
#define REPLACE_EXISTING 2 // replace only if tgt is not missing
#define SET_OR_APPEND    3 // set new value if missing or non-existent, append otherwise


#endif
