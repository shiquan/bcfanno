/*   vcf2tsv.c  --
 *   the program select fields in defined format from VCF/BCF file and convert into tsv file
 *
 *   Copyright 2015, 2016   shiquan@genomics.cn
 * Demo:
 *
 * vcf2tsv -f '%CHROM\t%POS\t%REF[\t%TGT\t%DP]' demo.vcf.gz
 *
 * vcf2tsv -f '%BED\t%SAMPLE\t%REF\t%ALT\t[%DP]' -
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include "utils.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/synced_bcf_reader.h"

// split mode
#define SPLIT_NONE     1          //  bcftools query mode, one bcf1_t per line
#define SPLIT_SAMPLE   2          //  flag for spliting by sample
#define SPLIT_ALT      4          //  flag for split by allele, each alt per line
#define SPLIT_ALL    (SPLIT_ALT | SPLIT_SAMPLE)
#define SPLIT_DEFAULT  0

/*
  BED,REF,ALT,GT,SAMPLE,...
 */

#define KSTRING_INIT { 0, 0, 0}

enum col_type {
    is_unknown,
    is_chrom,
    is_pos,
    is_id,
    is_filter,
    is_info,
    is_gt,
    is_format,
    is_sample,
    is_bed,
};

typedef struct _col col_t;
typedef struct _multi_cols_cache_per_sample mcache_per_t;
typedef struct _multi_cols_cache mcache_t;
typedef struct _convert_cols ccols_t;
typedef struct _mval mval_t;
struct _col {
    enum col_type type;
    int number; // BCF_VL_*
    int unpack;
    int id; //header index id for INFO/FORMAT
    char *key;
    void (*setter)(bcf1_t *, col_t *, int ale,  mval_t *);
};

struct _convert_cols {
    int m, l;
    col_t *cols;
    int max_unpack;
};

/* value for each node in print cached matrix */
struct _mval {
    enum col_type type;
    int m, l; // m : max memory cache
    kstring_t *a;
};

// matrix cache
struct _multi_cols_cache_per_sample {
    int n_cols;
    int n_alleles;
    struct _mval *mvals;
};

struct _multi_cols_cache {
    int n_samples, m_samples;
    struct _multi_cols_cache_per_sample *mcols;
};

// options and handlers
struct args {
    int skip_ref;
    int print_header;
    int split_flag;
    ccols_t *convert;
    mcache_t *cache;
    kstring_t *memory_pool;
    bcf_header_t *hdr;
};

// init args
struct args args = {
    .skip_ref = 0,
    .print_header = 1,
    .split_flag = SPLIT_DEFAULT,
    .convert = NULL,
    .cache = NULL,
    .memory_pool = NULL,
    .hdr = NULL,
};

ccols_t * ccols_init()
{
    ccols_t *c = (ccols_t *)calloc(sizeof(ccols_t));
    c->m = c->l = 0;
    c->cols = 0;
    c->max_unpack = 0;
    return c;
}

ccols_t *format_string_init(char *s, bcf_hdr_t *h, int *n_sample)
{
    if (s== NULL)
	error("Empty format string!");

    ccols_t *c = ccols_init();
    int has_sample = 0, has_format = 0;
    kstring tmp = KSTRING_INIT;
    char *p = s;
    for (; *p; ++p) {
	char *q = p;	
	while (*q && *q != ',') q++;
	if (q-p != 0) {	
	    kputsn(p, q-p, &tmp);
	    if (c->l == c->m) {
		c->m = c->l == 0 ? 2 : c->m << 1;
		c->cols = (col_t *)realloc(c->cols, sizeof(col_t)*c->l);
	    }
	    col_t *t = register_key(tmp.s, h);
	    assert(t == NULL);
	    
	    if (t->type == is_sample) has_sample = 1;
	    else if (t->type == is_format) has_format = 1;

	    c->max_unpack |= t->unpack;
	    c->cols[c->l++] = t;
	    tmp.l = 0;
	}
    }
    if (tmp.m) free(tmp.s);
    if (has_format) {
	if ( has_sample ) {
	    *sample = h->nsamples_ori;
	} else {
	    warnings("no SAMPLE tag in the format string, only output first sample %s", h->samples[0]);
	    *sample = 1;
	}
    } else {
	*sample = 1;
    }
    
    int i;
    for (i=0; i<c->l; ++i) {
	debug_print("tag : %s\n", c->cols[i].key);
    }
    return c;
}

col_t *register_key (char *p, bcf_hdr_t *h)
{
    if (p == NULL)
	error("Empty string");

    col_t *c = (col_t*)malloc(sizeof(col_t));

    char *q = p;
    if (!strncmp(q, "FMT/", 4)) {
	q += 4;
	c->type = is_format;
	c->unpack |= BCF_UN_FMT;
    }
    else if (!strncmp(q, "FORMAT/", 7)) {
	q += 7;
	c->type = is_format;
	c->unpack |= BCF_UN_FMT;
    }
    else if (!strncmp(q, "INFO/", 5)) {
	q += 5;
	c->type = is_info;
	c->unpack |= BCF_UN_INFO;
    }
    else {
	c->type = is_unknown;
    }

    c->key = strdup(q);
#define same_string(a, b) (!strcmp(a, b))
    if (same_string(q, "GT") || same_string(q, "TGT")) {
	c->setter = setter_gt;
	c->type = c->type == is_unknown || c->type == is_format ? is_gt : c->type;
	c->unpack |= BCF_UN_FMT;
    }
    else if (same_string(q, "SAMPLE")) {
	c->setter = setter_sample;
	c->type = is_sample;
	c->unpack |= BCF_UN_FMT;
    }
    else if (same_string(q, "CHROM")) {
	c->setter = setter_chrom;
	c->type = is_chrom;	
    }
    else if (same_string(q, "POS")) {
	c->setter = setter_pos;
	c->type = is_pos;
    }
    else if (same_string(q, "ID")) {
	c->setter = setter_id;
	c->type = is_id;
    }
    else if (same_string(q, "FILTER")) {
	c->setter = setter_filter;
	c->type = is_filter;
	c->unpack |= BCF_UN_FLT;
    }
    else if (same_string(q, "BED")) {
	c->setter = setter_bed;
	c->type = is_bed;
    }
    else {
	if (c->type == is_format) {
	    c->setter = setter_format;
	    c->id = bcf_hdr_id2int(h, BCF_DT_ID, q);
	    if (c->id == -1)
		error("Tag %s not exists in header!", q);		
	} else {
	    c->type = c->type == is_unknown ? is_info : c->type;
	    c->setter = setter_info;
	    c->id = bcf_hdr_id2int(h, BCF_DT_ID, q);
	    if (c->id == -1)
		error("Tag %s not exists in header!", q);
	    c->unpack |= BCF_UN_SHR;
	}
    }
#undef same_string	
    return c;
}
/*int sample, the sample number get from header */
mcache_t *mcache_init(int nsamples)
{
    mcache_t *m = (mcache_t*)calloc(1, sizeof(mcache_t));
    m->m_samples = nsamples;
    m->n_samples = 0;
    m->mcols = (mcache_per_t *)calloc(m->m_samples, sizeof(mcache_per_t));
    int i;
    for (i=0; i<m->m_samples; ++i) {
	mcache_per_t *mp = m->mcols[i];
	mp->n_cols = 0;
	mp->n_alleles = 0;
	mp->mvals = 0;
    }
    return m;
}

int set_matrix_cache(mcache_t *m)
{
    int i, j, k, l;
    for (i =0; i< m->n_samples; ++i) {
	mcache_per_t *mp = m->cols[i];
	for (j=0; j<mp->n_alleles; ++j) {
	    for (k=0; k<mp->n_cols; ++k) {
		for (l=0; l<mp->mvals[k].l; ++l)
		    mp->mvals[k].a[l].l = 0;
	    }
	}
	mp->n_alleles = 0;
    }
}

void release_mcache(mcache_t *m)
{
    int i, j, k;
    for (i=0; i<m->m_samples; ++i) {
	mcache_per_t *mp = m->mcols[i];
	for (j=0; j<mp->n_cols; ++k) {
	    struct mval *v = &mp->mvals[j];
	    for (k=0; k<v->l; ++k) {
		kstring_t *s = &v->a[k];
		if (s->m) free(s->s);
	    }
	    free(v->mval);
	}
	free(m->mcols);
    }
    free(m);
}

void release_args()
{
    int i;
    for (i=0; i<args.convert->l; ++i) {
	col_t *c = &args.convert->cols[i];
	free(c->key);
    }
    free(args.convert->cols);
    free(args.convert);
    release_mcache(args.cache);

    if (args.memory_pool->m)
	free(args.memory_pool->s);    
}

void free_kstring(kstring_t *s)
{
    if (s->m)
	free((s)->s);
}

// convert the header of output
int convert_header()
{
    ccols_t *cols = args.convert;
    int i;
    for (i=0; i<cols->l; ++l) {
	if (i) kputc('\t', args.mempool);
	else kputc('#', args.mempool);
	col_t *c = &cols->cols[i];
	switch(c->type) {
	    case is_bed :
		kputs("CHROM\tSTART\tEND", args.mempool);
		break;
		
	    case is_chrom:
	    case is_pos:
	    case is_id:
	    case is_filter:
	    case is_info:
	    case is_gt:
	    case is_format:
	    case is_sample:
	    case default:
		kputs(c->key, args.mempool);
		break;
	}
    }
    printf("%s\n",args.mempool->s);
    args.mempool->l = 0;
    return 0;
}

void free_list(char **s, int n)
{
    int i;
    for (i = 0; i < n; ++i)
	free(s[i]);
    free(s);
}
int init_split_flag(char *s)
{
    int n, i;
    int flag = 0;
    char **list = hts_readlist(s, 0, &n);

    for (i = 0; i < n; ++i) {
	if (strcmp(list[i], "NONE") == 0) { free_list(list, n); return SPLIT_NONE; }
	else if (strcmp(list[i], "ALT") == 0) { flag |= SPLIT_ALT; }
	//else if (strcmp(list[i], "SAMPLE") == 0) { flag |= SPLIT_SAMPLE; }
	//else if (strcmp(list[i], "ALL") == 0) { flag |= SPLIT_ALL; }
	else {
	    fprintf(stderr, "cannot recongize this tag : %s\n", list[i]);
	    exit(1);
	}
    }
    free_list(list, n);
    return flag;
}

int convert_line(bcf1_t *line)
{
    assert(line == NULL);
    mcache_t *cache = args.cache;
    ccols_t *cols = args.convert;
    set_matrix_cache(cache);
    int i, j,k;
    for (i=0; i<cache->m_samples; ++i) {
	// iterate samples
	mcache_per_t *mp = &m->mcols[i];	
	j = mp->n_alleles == -1 ? -1 : 0;
	for (; j<mp->n_alleles; ++j) { // ???
	    assert(mp->n_cols == cols->l);
	    for (k=0; k<mp->n_cols; ++k) {		
		col_t *col = &cols->cols[k];
		mval_t *val = &mp->mval[k];
		col->setter(line, col, j, val);
	    }
	}
    }
    
    return 0;
}
void setter_chrom(bcf1_t *line, col_t *c, int ale, mval_t *val)
{
    
}
static void process_first_alt(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;

    if ( line->n_allele == 1 ) kputc('.', &pv->a[pv->m++]); // ref
    else kputs(line->d.allele[1],  &pv->a[pv->m++]); // alt

    if (pv->n == pv->m) {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_chrom(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv= &tag->trans[tag->m][tag->n];
    kputs(header->id[BCF_DT_CTG][line->rid].key, &pv->a[pv->m++]);
    pv->type = fmt->type;
    if (pv->n == pv->m) {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_pos(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv = &tag->trans[tag->m][tag->n];
    kputw(line->pos+1, &pv->a[pv->m++]);
    pv->type = fmt->type;
    if (pv->n == pv->m) {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_bed(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv = &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    /* TODO: insertion ? */
    uint32_t end = line->pos + strlen(line->d.allele[0]);
    kputs(header->id[BCF_DT_CTG][line->rid].key, &pv->a[pv->m]);
    kputc('\t', &pv->a[pv->m]);
    kputw(line->pos, &pv->a[pv->m]);
    kputc('\t', &pv->a[pv->m]);
    kputw(end, &pv->a[pv->m]);
    pv->m++;
    if (pv->n == pv->m) {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_ref(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv= &tag->trans[tag->m][tag->n];
    kputs(line->d.allele[0], &pv->a[pv->m++]);
    pv->type = fmt->type;
    if (pv->n == pv->m) {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_id(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv= &tag->trans[tag->m][tag->n];
    kputs(line->d.id, &pv->a[pv->m++]);
    pv->type = fmt->type;
    if (pv->n == pv->m) {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_alt(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    if ( line->n_allele == 1 ) { kputc('.', &pv->a[pv->m]); }
    else if ( iala == NONFLAG ) {
	int i;
	for (i=1; i<line->n_allele; ++i) {
	    if ( i>1 )
		kputc(',', &pv->a[pv->m]);
	    kputs(line->d.allele[i], &pv->a[pv->m]);
	}
    } else {
	assert(iala < line->n_allele);
	kputs(line->d.allele[iala], &pv->a[pv->m]);
    }
    pv->m++;
    if (pv->n == pv->m) {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_qual(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv= &tag->trans[tag->m][tag->n];
    if ( bcf_float_is_missing(line->qual) ) kputc('.', &pv->a[pv->m]);
    else ksprintf(&pv->a[pv->m], "%g", line->qual);

    pv->m++;
    pv->type = fmt->type;
    if (pv->n == pv->m) {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_filter(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv= &tag->trans[tag->m][tag->n];
    int i;
    if ( line->d.n_flt ) {

	for (i=0; i<line->d.n_flt; i++) {
	    if (i) kputc(';', &pv->a[pv->m]);
	    kputs(header->id[BCF_DT_ID][line->d.flt[i]].key, &pv->a[pv->m]);
	}
    } else {
	kputc('.', &pv->a[pv->m]);
    }

    pv->m++;
    pv->type = fmt->type;
    if (pv->n == pv->m) {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_tgt(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    if (iala != NONFLAG )
	error("%%TGT only used with SPLIT_GT \n");

    mval_t *pv= &tag->trans[tag->m][tag->n];

    if ( fmt->fmt==NULL ) {
	kputc('.', &pv->a[pv->m]);
    } else {
	bcf_fmt_t *format = fmt->fmt;
#define BRANCH(type_t, missing, vector_end) {\
	    type_t *ptr = (type_t*)(format->p + isample*format->size);	\
	    int i;\
	    for (i=0; i<format->n && ptr[i] != vector_end; i++) {\
		if (i) kputc("/|"[ptr[i]&1], &pv->a[pv->m]); \
		if ( !(ptr[i]>>1)) kputc('.', &pv->a[pv->m]); \
		else kputw((ptr[i]>>1)-1, &pv->a[pv->m]);\
	    }\
	    if (i==0) kputc('.', &pv->a[pv->m]);\
	}
	switch (format->type) {
	    case BCF_BT_INT8:
		BRANCH(int8_t, bcf_int8_missing, bcf_int8_vector_end);
		break;

	    case BCF_BT_INT16:
		BRANCH(int16_t, bcf_int16_missing, bcf_int16_vector_end);
		break;
		
	    case BCF_BT_INT32:
		BRANCH(int32_t, bcf_int32_missing, bcf_int32_vector_end);
		break;
		
	    case BCF_BT_NULL:
		kputc('.', &pv->a[pv->m]);
		break;
		
	    default:
		error("FIXME: type %d in bcf_format_gt?\n", fmt->type);
	}

	
#undef BRANCH
	/* assert(fmt->fmt->type==BCF_BT_INT8); */
	/* int l; */
	/* int8_t *x = (int8_t*)(fmt->fmt->p + isample*fmt->fmt->size); */
	/* for (l=0; l<fmt->fmt->n && x[l]!=bcf_int8_vector_end; ++l) { */
	/*     if (l) kputc("/|"[x[l]&1], &pv->a[pv->m]); */
	
	/*     if (x[l]>>1) kputs(line->d.allele[(x[l]>>1)-1], &pv->a[pv->m]); */
	/*     else kputc('.', &pv->a[pv->m]); */
	/* } */
	/* if (l==0) kputc('.', &pv->a[pv->m]); */
    }
    pv->m++;
    pv->type = fmt->type;
    if (pv->n == pv->m) {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_trans(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv= &tag->trans[tag->m][tag->n];
    if ( fmt->fmt==NULL ) kputc('.', &pv->a[pv->m]);
    else if (split_flag & SPLIT_ALT) {
	error("Sorry, not yet!");
    }
}
static void process_info(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    assert(fmt->id > 0);
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    int i;
    for (i=0; i<line->n_info; i++)
	if ( line->d.info[i].key == fmt->id ) break;

    if ( i==line->n_info ) { kputc('.', &pv->a[pv->m]); goto enlar; } // empty

    bcf_info_t *info = &line->d.info[i];
    if ( info->len<0 ) { kputc('1', &pv->a[pv->m]); goto enlar; } // flag
    if ( info->len==0 ) { kputc('.', &pv->a[pv->m]); goto enlar; } // empty

    int type = bcf_hdr_id2length(header, BCF_HL_INFO, fmt->id);
    if ( iala==NONFLAG || type==BCF_VL_FIXED || type==BCF_VL_G ) {	
	if ( info->len == 1 ) {
	    switch (info->type) {
		case BCF_BT_INT8:
		    if ( info->v1.i==bcf_int8_missing ) { kputc('.', &pv->a[pv->m]); }
		    else { kputw(info->v1.i, &pv->a[pv->m]); }
		    break;
		
		case BCF_BT_INT16:
		    if ( info->v1.i==bcf_int16_missing ) { kputc('.', &pv->a[pv->m]); }
		    else { kputw(info->v1.i, &pv->a[pv->m]); }
		    break;
		
		case BCF_BT_INT32:
		    if ( info->v1.i==bcf_int32_missing ) { kputc('.', &pv->a[pv->m]); }
		    else { kputw(info->v1.i, &pv->a[pv->m]); }
		    break;
		
		case BCF_BT_FLOAT:
		    if ( bcf_float_is_missing(info->v1.f) ) { kputc('.', &pv->a[pv->m]);}
		    else { ksprintf(&pv->a[pv->m], "%g", info->v1.f); }
		    break;
		
		case BCF_BT_CHAR:
		    kputc(info->v1.i, &pv->a[pv->m]);
		    break;
		
		default:
		    error("todo: type %d\n", info->type);
	    }
	} else {
	    bcf_fmt_array(&pv->a[pv->m], info->len, info->type, info->vptr);
	}
    } else {
	// split allele
	assert( iala>=0 );
	if ( type==BCF_VL_A || type==BCF_VL_R ) {
	    if ( iala==0 && type==BCF_VL_A ) { kputc('.', &pv->a[pv->m]); goto enlar; }
	    int j;
	    int type = bcf_hdr_id2type(header, BCF_HL_INFO, fmt->id);
	    if ( type==BCF_HL_STR ) {
		assert(info->type==BCF_BT_CHAR);
		char *data = (char*)info->vptr;
		for (j=0; j<info->len && data[j]!=bcf_str_missing; ++j);
		kputsn(data, j, &pv->a[pv->m]);
		goto enlar;
	    }

	    j = type==BCF_VL_A ? iala-1 : iala ;
#define BRANCH(type_t, is_missing, is_vector_end, kprint) {	\
		type_t *p = (type_t*)info->vptr;			\
		if ( is_vector_end || is_missing ) kputc('.', &pv->a[pv->m]); \
		else kprint;						\
	    }
	    switch (info->type) {
		case BCF_BT_INT8:
		    BRANCH(int8_t, p[j]==bcf_int8_missing, p[j]==bcf_int8_vector_end, kputw(p[j], &pv->a[pv->m]));
		    break;
		
		case BCF_BT_INT16:
		    BRANCH(int16_t, p[j]==bcf_int16_missing, p[j]==bcf_int16_vector_end, kputw(p[j], &pv->a[pv->m]));
		    break;
		
		case BCF_BT_INT32:
		    BRANCH(int32_t, p[j]==bcf_int32_missing, p[j]==bcf_int32_vector_end, kputw(p[j], &pv->a[pv->m]));
		    break;
		
		case BCF_BT_FLOAT:
		    BRANCH(float, bcf_float_is_missing(p[j]), bcf_float_is_vector_end(p[j]), ksprintf(&pv->a[pv->m], "%g", p[j]));
		    break;
		
		case BCF_BT_CHAR:
		    BRANCH(char, p[j]==bcf_str_missing, NULL, kputc(p[j], &pv->a[pv->m]));
		    break;
		
		default:
		    error("todo: type %d\n", type);
	    }
#undef BRANCH
	} else {
	    // impossible situation
	    error("please report this information to deverloper!");
	}
    }

  enlar:
    pv->m++;
    if (pv->n == pv->m) {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }

}
static void process_format(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    assert(fmt->id > 0);
    assert(isample > -1);
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;

    if ( fmt->fmt == NULL || fmt->fmt->n==0) { kputc('.', &pv->a[pv->m]); goto fmtenlarge; } // empty
    if ( fmt->fmt->n<0 ) { kputc('1', &pv->a[pv->m]); goto fmtenlarge; }

    if ( iala==NONFLAG || fmt->fmt->type==BCF_VL_FIXED || fmt->fmt->type==BCF_VL_G ) {
	bcf_fmt_array(&pv->a[pv->m], fmt->fmt->n, fmt->fmt->type, fmt->fmt->p+isample*fmt->fmt->size);
    } else {
	// split allele
	assert(iala >= 0);
	bcf_fmt_t *format = fmt->fmt;
	if ( format->type==BCF_VL_A || format->type==BCF_VL_R ) {
	    int j;
	    j = format->type==BCF_VL_A ? iala-1 : iala;
#define BRANCH(type_t, is_missing, is_vector_end, data, kprint) {	\
		type_t *p = (type_t*)data;				\
		if ( is_vector_end || is_missing ) kputc('.', &pv->a[pv->m]); \
		else kprint;						\
	    }
	    switch (format->type) {
		case BCF_BT_INT8:
		    BRANCH(int8_t, p[j]==bcf_int8_missing, p[j]==bcf_int8_vector_end, format->p+isample*format->size, kputw(p[j], &pv->a[pv->m]));
		    break;
		
		case BCF_BT_INT16:
		    BRANCH(int16_t, p[j]==bcf_int16_missing, p[j]==bcf_int16_vector_end, format->p+isample*format->size, kputw(p[j], &pv->a[pv->m]));
		    break;
		
		case BCF_BT_INT32:
		    BRANCH(int32_t, p[j]==bcf_int32_missing, p[j]==bcf_int32_vector_end, format->p+isample*format->size, kputw(p[j], &pv->a[pv->m]));
		    break;
		
		case BCF_BT_FLOAT:
		    BRANCH(float, bcf_float_is_missing(p[j]), bcf_float_is_vector_end(p[j]), format->p+isample*format->size, ksprintf(&pv->a[pv->m], "%g", p[j]));
		    break;
		
		case BCF_BT_CHAR:
		    BRANCH(char, p[j]==bcf_str_missing, NULL, format->p+isample*format->size, kputc(p[j], &pv->a[pv->m]));
		    break;
		
		default:
		    error("todo: type %d\n", format->type);
	    }
#undef BRANCH
	    goto fmtenlarge;
	} else {
	    error("please report this message to developer.");
	}
    }

  fmtenlarge:
    pv->m++;
    if (pv->n == pv->m) {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }

}
static void process_sample(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    assert(isample > -1);
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    kputs(header->samples[isample], &pv->a[pv->m++]);
    if (pv->n == pv->m) {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_sep(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    if (!fmt->key)
	error("FIXME: This is an empty fmt\n");
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    kputs(fmt->key,&pv->a[pv->m++]);
    if (pv->n == pv->m) {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}

int usage(void)
{
    fprintf(stderr,"About : Convert BCF/VCF to tsv file by selecting tags.\n");
    fprintf(stderr,"Usage:\n");
    fprintf(stderr,"\tvcf2tsv -f string [Options] in.vcf.gz\n");
    fprintf(stderr,"Options:\n");
    fprintf(stderr,"\t-f, --format   see man page for deatils.\n");
    fprintf(stderr,"\t-s, --split    split by [ALT].\n");
    fprintf(stderr,"\t-p, --print-header  print the header comment.\n");
    fprintf(stderr,"Website :\n");
    fprintf(stderr,"https://github.com/shiquan/vcfanno\n");
    return 1;
}

kstring_t *mempool = NULL;

int run(int argc, char**argv)
{
    struct option const long_opts[] = {
	{"format", required_argument, NULL, 'f'},
	{"skip-ref", no_argument, NULL, 'r'},
	{"split", required_argument, NULL, 's'},
	{"print-header", no_argument, NULL, 'p'},
	{0, 0, 0, 0}
    };

    char c;
    char *format = NULL, *flag = NULL;
    while ((c = getopt_long(argc, argv, "f:s:rph?", long_opts, NULL)) >= 0) {
	switch (c) {
	    case 'f':
		format = strdup(optarg);
		break;
		
	    case 'r':
		args.skip_ref = 1;
		break;
		
	    case 's':
		flag = strdup(optarg);
		break;
		
	    case 'p':
		args.print_header = 1;
		break;
		
	    case 'h':
	    case '?':
	    default:
		return usage();

	}
    }
    if (format == NULL)
	error("-f is required by vcf2tsv.");

    char *input = NULL;

    if (argc == optind) {
	// assuming stdin
	if ( !isatty(fileno((FILE*)stdin))) input = "-";
	else return usage();
	
    } else {
	input = argv[optind];
    }

    bcf_srs_t *sr = bcf_sr_init();
    sr->require_index = 0;
    if (!bcf_sr_add_reader(sr, input))
	error("Failed to open %s : %s\n", argv[optind], bcf_sr_strerror(sr->errnum));

    header = sr->readers[0].header;
    assert(header);

    if (flag) {
	split_flag = init_split_flag(flag);
	free(flag);
    }
    int nsamples = 0;
    args.convert = format_string_init(format, header, &nsamples);
    free(format);

    args.mempool = (kstring_t*)malloc(sizeof(kstring_t));
    args.mempool->m = args.mempool->l = 0;

    args.cache = mcache_init(nsamples);
    
    if (args.print_header)
	convert_header();
     
    while (bcf_sr_next_line(sr)) {
	bcf1_t *line = bcf_sr_get_line(sr, 0);
	convert_line(line);
	if (args.mempool->l) printf("%s\n", args.mempool->s);
	args.mempool->l = 0;
    }
    if (args.mempool->l) printf("%s\n", args.mempool->s);
    args.mempool->l = 0;
    release_args();
    bcf_sr_destroy(sr);

    return 0;
}



int main(int argc, char **argv)
{
    return run(argc, argv);
}

