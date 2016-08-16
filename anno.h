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
#include <htslib/khash_str2int.h>
#include <htslib/tbx.h>
#include <htslib/faidx.h>
#include "plugin.h"
#include "config.h"
#include "vcmp.h"

#define ANNOTCOLS_PACK_INIT { NULL, 0, 0, NULL }
#define KSTRING_INIT { 0, 0, 0 }

/* annot line also defined in plugin.h for compile dynamic library */
#ifndef ANNOT_LINE_FLAG
#define ANNOT_LINE_FLAG

struct annot_line {
    char **cols;
    int ncols, mcols;
    int nals, mals;
    char **als;
    kstring_t line;
    int rid, start, end;
};

#endif

struct annot_bed_line {
    int rid;
    int start; // 0based start
    int end; // 1based end
    char *a;
};

typedef void (* rel_func)(void*);

extern void safe_release(void * p, rel_func func);

/* struct annot_cols_pack { */
/*     struct annot_col *cols; */
/*     int ncols; */
/*     enum anno_type type;  */
/*     char *columns; */
/* }; */
// cached regions
struct region1 {
    uint32_t start;
    uint32_t end; // 0 for init
};

struct region {
    int n, m, c;
    struct region1 *regs;
};

/*
 * local HGVS nomenclature convert,
 * see refgene.h for online convert method..
 */


enum strand { strand_plus, strand_minus, strand_unknonw, };

extern struct annot_col *init_columns(const char *rules, bcf_hdr_t *in, bcf_hdr_t *out, int *n, enum anno_type type);

struct sql_connect {
    void *connect;
    // index id
    int idx; 
    char *column;
    int ncols;
    struct annot_col *cols;
    int nbuffers, mbuffer;
    struct annot_line *buffers;
    // in memory cached regions
    struct region *regs; 
    // current position: chr name index to names[]
    int iseq; 
    int start, end;
    int prev_seq, prev_start;
};

struct annot_cols {
    // for vcf files, data_id should be the index of bcf_srs_files
    int data_id;
    // number of tags
    int ncols;
    // col structure for each tag
    struct annot_col *cols;
    // columns string, retrieve from configure file
    const char *columns;
};


// for tabix-indexed database, usually four columns, three bed format with an extra tag column.
//
struct annot_local_bed {
    // file name
    const char *fname;
    // tabix indexed cache
    tbx_t *tbx;
    // file handler
    htsFile *fp;
    // only one col is accept
    struct annot_col cols;
    
};

struct anno_aux {
    // hdr is header struct of input vcf file, all the annotated tags should be inited in the hdr_out before export bcf lines
    bcf_hdr_t *hdr, *hdr_out;
    // file handler of input vcf
    htsFile *fp; 
    // file handler of output vcf
    htsFile *out_fh;
    // output directory, default is "-" for stdout
    const char *out;
    // output format, default is vcf
    int output_type;
    // pre-indexed vcf/bcf databases
    bcf_srs_t *files;
    // columns structure of vcf databases
    struct annot_cols *cols;
    
    int i_data; 
    // for matching annotation and VCF lines by allele
    vcmp_t *vcmp;  

    int *sample_map, nsample_map, sample_is_file;
    
    int mtmpi, mtmpf, mtmps;
    int mtmpi2, mtmpf2, mtmps2;
    int mtmpi3, mtmpf3, mtmps3;
    int32_t *tmpi, *tmpi2, *tmpi3;
    float *tmpf, *tmpf2, *tmpf3;
    char *tmps, *tmps2, **tmpp, **tmpp2;
    kstring_t tmpks;
};

#define REPLACE_MISSING  0 // replace only missing values
#define REPLACE_ALL      1 // replace both missing and existing values
#define REPLACE_EXISTING 2 // replace only if tgt is not missing
#define SET_OR_APPEND    3 // set new value if missing or non-existent, append otherwise

struct annot_col {
    int icol, replace, number;  // number: one of BCF_VL_* types
    char *hdr_key;
    int (*setter)(struct anno_aux *, bcf1_t *, struct annot_col *, void*);
};

extern int setter_filter(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int vcf_setter_filter(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int setter_id(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int vcf_setter_id(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int setter_qual(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int vcf_setter_qual(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int setter_info_flag(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int vcf_setter_info_flag(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int setter_ARinfo_int32(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, int nals, char **als, int ntmpi);
extern int setter_info_int(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int vcf_setter_info_int(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int setter_ARinfo_real(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, int nals, char **als, int ntmpf);
extern int setter_info_real(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int vcf_setter_info_real(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst);
extern int setter_ARinfo_string(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, int nals, char **als);
extern int setter_info_str(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int vcf_setter_info_str(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int vcf_setter_format_gt(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int count_vals(struct annot_line *tab, int icol_beg, int icol_end);
extern int setter_format_int(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int setter_format_real(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int setter_format_str(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int vcf_setter_format_int(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int vcf_setter_format_real(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);
extern int vcf_setter_format_str(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data);

/**/
extern void rebuild_annot_lines(struct anno_aux *hand);


#endif
