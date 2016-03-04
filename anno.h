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

/* memory handle macros */

#define LINE_CACHE 1024

#define check_double_free(p) do						\
    {                                                                   \
        void **_pp = (void**)&(p);                                      \
        if (_pp==NULL || *_pp == NULL) { \
	    fprintf(stderr,"[error] double free %s, %d", __FUNCTION__, __LINE__); \
	    exit(EXIT_FAILURE);						\
	}								\
    } while(0)

#define ignore_free(p) do \
    {\
	void **_pp = (void**)&(p);					\
	if (*_pp!=NULL && _pp != NULL) {				\
	    free(*_pp);							\
	    *_pp = NULL;						\
	}								\
    } while(0)

#define safe_free(p) do \
    {\
	void **_pp = (void**)&(p);	 \
        if (_pp==NULL || *_pp == NULL) {				\
	    fprintf(stderr,"[error] double free %s, %d", __FUNCTION__, __LINE__); \
	    exit(EXIT_FAILURE);						\
	}								\
	free(*_pp);							\
        *_pp = NULL;                                                    \
    } while(0)

#define check_mem(p) do \
    {\
	void **_pp = (void**)&p; \
	if (_pp == NULL || *_pp == NULL) {\
	    fprintf(stderr, "[memory out] func: %s, line: %d\n", __FUNCTION__, __LINE__);\
	    exit(EXIT_FAILURE);						\
	}\
    }while(0)

typedef void (* rel_func)(void*);

extern void safe_release(void const * p, rel_func func);

/* ref: http://c.learncodethehardway.org/book/ex20.html */
#define str_errno() (errno == 0 ? "None" : strerror(errno))

#define clear_errno() do \
    {\
	fprintf(stderr, "%s\n", str_errno());\
	errno = 0;\
    }while(0)

#define error(line, ...) do						\
    {									\
	fprintf(stderr, "[error] func:%s, line:%d, errno:%s. " line "\n", __FUNCTION__, __LINE__, str_errno(), ##__VA_ARGS__); \
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
	fprintf(stderr, "[debug] func:%s, line:%d, errno:%s. " line "\n", __FUNCTION__, __LINE__, str_errno(), ##__VA_ARGS__); \
    } while(0)


/*
 * local HGVS nomenclature convert,
 * see refgene.h for online convert method..
 */

/* See UCSC website for more details about refgene format, there is no 
 * version information in refgene file in default, however, we can retrieve
 * this version number from other file in the same directrary. But remeber 
 * to check the names accordly in the refrna file.
 */
struct refgene_file {
    const char *file_path; 
    const tbx_t *idx; // refgene should be sorted and indexed by tabix
};

/* refrna should be indexed by samtools faidx, and the names should be
 * keep consistency with refgene entries 
 */
struct refrna_file {
    const char *file_path;
    const faidx_t *idx;
};

enum strand { strand_plus, strand_minus, strand_unknonw, };

struct region {
    uint32_t start;
    uint32_t stop;
};

#define clear_all(x) (x ＝ 0)
#define clear_flag(x, n) (x &= ~(n))

#define REFGENE_PRASE_BIN    1
#define REFGENE_PRASE_NAME1  2
#define REFGENE_PRASE_REG    4
#define REFGENE_PRASE_STRAND 8
#define REFGENE_PRASE_EXONS  (1<<4 | 8)
#define REFGENE_PRASE_ALL (REFGENE_PRASE_BIN | REFGENE_PRASE_NAME1 | REFGENE_PRASE_REG| REFGENE_PRASE_EXONS)

struct refgene_entry {
    char *buffer;
    int n_buffer;
    int m_buffer;
    int prase_flag;
    
    int bin; // see USCS bin scheme for details, be used to check regions
    
    const char *name1; // gene name
    const char *name2; // rna name usually

    int rid; // chromosome id, should be one of contigs in the VCF header, -1
    
    int strand; // strand_plus, strand_minus, strand_unknown
    uint32_t exon_begin; // exon region begin from xxx
    uint32_t exon_end; // exon end to xxx, remember exon contains cds and utr
    uint32_t cds_start; // cds region start position, 5'UTR is exon_begin to cds_start-1.
    uint32_t cds_stop; // 3'UTR is cds_end+1 to exon_end in plus strand
    int total_exons; // exon number
    struct region *exons;
    struct region *cds;
};

struct hgvs_annos {
    const char *columns;
    struct refgene_file refgene;
    struct refrna_file refrna;
};

extern void extract_refgene(struct refgene_entry *entry, int type);
extern bcf1_t *hgvs_anno_local(struct refgene_entry *entry, bcf1_t *line);


#endif
