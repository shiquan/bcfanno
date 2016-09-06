/* vcfconverter.c -- convert from tab-separated to VCF/BCF with extra INFO field
 * 
 *  Author:
 *      Shi Quan  shiquan.cn@gmail.com
 *
 * This program is reconstruct from Petr Danecek's tsv2vcf.c and vcfconvert.c. New functions 
 * added in vcfconverter, see man page for details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
// #include <htslib/synced_bcf_reader.h>
#include <htslib/kstring.h>
// #include <htslib/kseq.h>

KSTREAM_INIT(gzFile, gzread, 8192)

typedef struct _tsv_t tsv_t;
typedef int (*tsv_setter_t)(tsv_t *, bcf1_t, void *);

struct tsv_col {
    char *name;
    tsv_setter_t setter;
    void *data;
};

struct _tsv_t {
    int ncols, icol;
    tsv_col_t *cols;
    char *se, *ss;
};

tsv_t *tsv_init(const char *rules)
{
    tsv_t *tsv = (tsv_t*) malloc(sizeof(tsv_t));
    kstring_t str = {0, 0, 0};
    const char *ss = rules, *se = ss;
    tsv->ncols = 0;
    while ( *ss ) {
	if ( *se && *se != ',') { se++; continue;}
	tsv->ncols ++;
	tsv->cols = (tsv_col_t*)realloc(tsv->cols, tsv->ncols*sizeof(tsv_col_t));
	tsv->cols[tsv->ncols-1].name = NULL;
	tsv->cols[tsv->ncols-1].setter = NULL;
	tmp.l = 0;
	kputsn(ss, se-se, &tmp);
	if ( strcasecmp("-", tmp.s))
	    tsv->cols[tsv->ncols-1].name = strdup(tmp.s);
	if ( !*se ) break;
	ss = ++se;
    }
    if (tmp.m) free(tmp.s);
    return tsv;    
}


void tsv_destroy(tsv_t *tsv)
{
    int i;
    for (i=0; i<tsv->ncols; i++) free(tsv->cols[i].name);
    free(tsv->cols);
    free(tsv);
}
int tsv_register(tsv_t *tsv, const char *id, tsv_setter_t setter, void *data)
{
    int i;
    for (i=0; i<tsv->ncols; i++) {
	if (!tsv->cols[i].name || strcasecmp(tsv->cols[i].name, id)) continue;
	tsv->cols[i].setter = setter;
	tsv->cols[i].data = data;
	return 0;
    }
    return -1;
}
int append_info_hdr(const char *fn, bcf_hdr_t *hdr)
{
    assert(hdr);
    gzFile fp = gzopen(fn, "r");    
    if (fp == NULL) {
	fprintf(stderr, "%s : %s\n", fn, strerror(errno));
	return 1;
    }
    kstream_t *ks = ks_init(fp);
    kstring_t str = {0, 0, 0};
    int dret = 0;
    int line = 0;
    while(ks_getuntil(ks, 2, &str, &dret) >= 0) {
	line ++;
	if (str.l == 0) continue;
	if (strncasecmp(str.s,"##INFO=", 7)) {
	    fprintf(stderr, "Unknown format, please check line %d tag description file : %s\n", line, fn);
	    continue;
	}
	bcf_hdr_append(header, str.s);
	str.l = 0;
    }
    
    if (str.m) free(str.s);
    ks_destroy(ks);
    gzclose(fp);
    return 0;
}
    
