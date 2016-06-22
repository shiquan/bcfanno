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

typedef annot_line annot_line_t;

#endif

/* memory handle macros */
#define LINE_CACHE 1024
#define check_double_free(p) do						\
    {                                                                   \
        void **_pp = (void**)&(p);                                      \
        if (_pp==NULL || *_pp == NULL) {				\
	    fprintf(stderr,"[error] double free %s, %d", __FUNCTION__, __LINE__); \
	    exit(EXIT_FAILURE);						\
	}								\
    } while(0)

#define ignore_free(p) do						\
    {									\
	void **_pp = (void**)&(p);					\
	if (*_pp!=NULL && _pp != NULL) {				\
	    free(*_pp);							\
	    *_pp = NULL;						\
	}								\
    } while(0)

#define safe_free(p) do				\
    {						\
	void **_pp = (void**)&(p);					\
        if (_pp==NULL || *_pp == NULL) {				\
	    fprintf(stderr,"[error] double free %s, %d", __FUNCTION__, __LINE__); \
	    exit(EXIT_FAILURE);						\
	}								\
	free(*_pp);							\
        *_pp = NULL;                                                    \
    } while(0)

#define check_mem(p) do				\
    {						\
	void **_pp = (void**)&p;		\
	if (_pp == NULL || *_pp == NULL) {				\
	    fprintf(stderr, "[memory out] func: %s, line: %d\n", __FUNCTION__, __LINE__);\
	    exit(EXIT_FAILURE);						\
	}								\
    }while(0)

typedef void (* rel_func)(void*);

extern void safe_release(void * p, rel_func func);

/* ref: http://c.learncodethehardway.org/book/ex20.html */
#define str_errno() (errno == 0 ? "None" : strerror(errno))

#define clear_errno() do \
    {\
	fprintf(stderr, "%s\n", str_errno());\
	errno = 0;\
    }while(0)

#define error(line, ...) do						\
    {									\
	fprintf(stderr, "[error] func : %s, line : %d, errno : %s. " line "\n", __FUNCTION__, __LINE__, str_errno(), ##__VA_ARGS__); \
	errno = 0;							\
	exit(EXIT_FAILURE);						\
    }while(0)

#define warnings(line, ...) do						\
    {									\
	if (errno == 0) {						\
	    fprintf(stderr, "[warnings] " line "\n", ##__VA_ARGS__);	\
	} else {							\
	    fprintf(stderr, "[warnings] Errno: %s. " line "\n", str_errno(), ##__VA_ARGS__); \
	}								\
    }while(0)

#define debug_print(line, ...) do						\
    {									\
	fprintf(stderr, "[debug] func : %s, line : %d, errno : %s. " line "\n", __FUNCTION__, __LINE__, str_errno(), ##__VA_ARGS__); \
    } while(0)


typedef struct annot_col annot_col_t;

struct annot_cols_pack
{
    annot_col_t *cols;
    int ncols;
    enum anno_type type; 
    char *columns;
};

typedef struct
{
    uint32_t start, end; // 0 for init
}
region1_t;

struct region
{
    region1_t *regs;
    int n, m, c;
};

/*
 * local HGVS nomenclature convert,
 * see refgene.h for online convert method..
 */

/* See UCSC website for more details about refgene format, there is no 
 * version information in refgene file in default, however, we can retrieve
 * this version number from other file in the same directrary. But remeber 
 * to check the names accordly in the refrna file.
 */
struct refgene_file_option {
    const char *file_path; 
    const tbx_t *idx; // refgene should be sorted and indexed by tabix
};

/* refrna should be indexed by samtools faidx, and the names should be
 * keep consistency with refgene entries 
 */
struct refrna_file_option {
    const char *file_path;
    const faidx_t *idx;
};

enum strand { strand_plus, strand_minus, strand_unknonw, };

extern annot_col_t *init_columns(const char *rules, bcf_hdr_t *in, bcf_hdr_t *out, int *n, enum anno_type type);

struct sql_connect
{
    void *connect;
    int idx; // index id
    char *column;
    int ncols;
    annot_col_t *cols;
    int nbuffers, mbuffer;
    annot_line_t *buffers;
    /* kstring_t line; */
    /* char *fname; */
    /* char **als; // parsed alleles  */
    /* kstring_t als_str; */
    /* int nals, mals;   */
    /* int als_type;  // alleles type, VCF_SNP or VCF_INDEL */
    
    struct region *regs; // in memory cached regions
    /* void *seq_hash; */
    /* char **seq_names; */
    /* int nseqs; */
    int iseq; // current position: chr name index to names[]
    int start, end;
    int prev_seq, prev_start;
};

/* struct sql_connect_readers */
/* { */
/*     struct sql_connect *connects; */
/*     int nreaders; */
/* }; */
/* anno_handler is adapt form args_t struct in vcfannotation.c */
/* struct filter_pack */
/* { */
/*     filter_t *filter; */
/*     char *string; */
/*     annot_col_t *cols; */
/*     int ncols; */
/* }; */

struct annot_cols {
    int ncols;
    annot_col_t *cols;
    //char *columns;
};

struct annot_cols_vector {
    int n, m;
    struct annot_cols *vcols;
};

struct anno_handler
{
    //int vcf_anno_count;
    //int sql_anno_count;
    bcf_hdr_t *hdr, *hdr_out;
    htsFile *out_fn;
    char *out;
    int output_type, n_threads;
    //struct anno_hgvs_option *hgvs_opts;
    bcf_srs_t *files;
    struct annot_cols_vector *vcf_cols;
    struct annot_cols_vector *sql_cols;

    int ti; 
    //struct annot_cols_pack *cols;
    //struct sql_connect *connects;

    /* filter tag, usually in FORMAT field */
    //struct filter_pack *filter;
    
    vcmp_t *vcmp;  // for matching annotation and VCF lines by allele
    //annot_line_t *alines; // buffered annotation lines
    //int nalines, malines;
    //int ref_idx, alt_idx, from_idx, to_idx; // -1 for not present
    //struct annot_col *cols;
    //int ncols;

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

struct annot_col
{
    int icol, replace, number;  // number: one of BCF_VL_* types
    char *hdr_key;
    int (*setter)(struct anno_handler *, bcf1_t *, annot_col_t *, void*);
};

extern int setter_filter(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int vcf_setter_filter(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int setter_id(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int vcf_setter_id(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int setter_qual(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int vcf_setter_qual(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int setter_info_flag(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int vcf_setter_info_flag(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int setter_ARinfo_int32(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, int nals, char **als, int ntmpi);
extern int setter_info_int(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int vcf_setter_info_int(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int setter_ARinfo_real(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, int nals, char **als, int ntmpf);
extern int setter_info_real(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int vcf_setter_info_real(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst);
extern int setter_ARinfo_string(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, int nals, char **als);
extern int setter_info_str(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int vcf_setter_info_str(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int vcf_setter_format_gt(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int count_vals(annot_line_t *tab, int icol_beg, int icol_end);
extern int setter_format_int(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int setter_format_real(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int setter_format_str(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int vcf_setter_format_int(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int vcf_setter_format_real(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);
extern int vcf_setter_format_str(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data);

/**/
extern void rebuild_annot_lines(struct anno_handler *hand);


#endif
