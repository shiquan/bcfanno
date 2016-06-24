/*   select.c  --
 *   This is a plugin of bcftools for selecting fields in defined format from VCF/BCF file
 *
 *   The extracting rules are most similar with `bcftools-query -f`, but still
 *   there are two different rules in this plugin. 
 *   1) each sample can be exported per line by specify %SAMPLE in the INFO region
 *   2) split entry by specified tag  '--split NONE,ALT,GT,SAMPLE,HGVS,ALL'
 *
 *   Copyright (C) 2015 BGI Research
 *
 *   Author : SHI Quan (shiquan@genomics.cn)
 *
 *   License : MIT
 *
 *   This program is inspired by pd3's vcfquery.c .
 *
 * 
 * Demo:
 *
 * bcfselect -f '%CHROM\t%POS\t%REF[\t%TGT\t%DP]' demo.vcf.gz
 *
 * bcfselect -f '%BED\t%SAMPLE\t%REF\t%ALT\t[%DP]' demo.vcf.gz
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
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

// split mode
#define SPLIT_NONE     1          //  bcftools query mode, one bcf1_t per line
//#define SPLIT_TAG      0
#define SPLIT_SAMPLE   2          //  flag for spliting by sample
#define SPLIT_ALT      4          //  flag for split by allele, each alt per line
//#define SPLIT_GT       2          //  flag for split by genotype (0/0, 0/1, 1/1), it's same with split by samples ^ SPLIT_ALT
//#define SPLIT_TRANS    (1<<4)     //  flag for split by TRANS name, if use this flag any tags related with transcripts should be split 
#define SPLIT_DEFAULT  0
//#define SPLIT_ALL      (SPLIT_SAMPLE | SPLIT_ALT | SPLIT_TRANS)

#define TAG_UNKNOWN   0
#define TAG_CHROM     1
#define TAG_POS       2
#define TAG_ID        4
#define TAG_REF       8
#define TAG_ALT       16
#define TAG_BED       32
#define TAG_QUAL      64
#define TAG_FILTER    128
#define TAG_INFO      256
#define TAG_FORMAT    513
#define TAG_SAMPLE    1024
#define TAG_SEP       2048
#define TAG_TYPE      4096
#define TAG_GT        8192
#define TAG_TGT       8192
//#define TAG_IUPAC_GT  16384
//#define TAG_FIRST_ALT 32768
//#define TAG_TRANS     65536

#define IS_SEP    (TAG_SEP)
#define IS_SAM    (TAG_SAMPLE)
//#define IS_ALT    (TAG_BED | TAG_ALT | TAG_INFO | TAG_FORMAT | TAG_TRANS)
#define IS_ALT    (TAG_BED | TAG_ALT | TAG_INFO | TAG_FORMAT)
#define IS_GT     (TAG_GT | TAG_TGT)
//#define IS_TRANS  (TAG_TRANS)

#define MEMPOOL     2621440

#define NONFLAG  -2

#define str_errno() (errno == 0 ? "None" : strerror(errno))

#define error(line, ...) do						\
    {									\
	fprintf(stderr, "[error] func : %s, line : %d, errno : %s. " line "\n", __FUNCTION__, __LINE__, str_errno(), ##__VA_ARGS__); \
	errno = 0;							\
	exit(EXIT_FAILURE);						\
    } while(0)

/* split flag */
static int split_flag = SPLIT_DEFAULT;

/* bcf_hdr_t * header; */

/* the fmt_t and convert_t are adapted from pd3's convert.c 
 * if we want export the sample info and allele info per line
 * we must rewrite the process_* functions in convert.c
 */
/* typedef struct _convert convert_t; */
/* typedef struct _tags     tags_t; */

/* struct _fmt { */
/*     int id, type, is_format; */
/*     char *key; */
/*     bcf_fmt_t *fmt; */
/*     void (*handler)(bcf1_t *, struct _fmt *, int l, int isample, tags_t *); */
/*     // int l, allele */
/* }; */

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
    is_sep, 
};

enum field_type {
    field_fix,
    field_info,
    field_format,
    field_error,
};

typedef struct _col col_t;
typedef struct _matrix_cache mcache_t;
typedef struct _convert_cols ccols_t;

struct _col {
    enum col_type type;
    char *key;
    void (*setter)(bcf1_t *, col_t *, int ale,  mcache_t*);
};

struct _convert_cols {
    int m, n;
    col_t *cols;
    int max_unpack;
    char *format_string;
};

struct matrix_value {
    enum col_type type;
    int m, n; // m: max memory cache
    kstring_t *a; 
};

/* typedef struct _fmt fmt_t; */

/* struct _convert { */
/*     int mfmt, nfmt; // mfmt: max memory cache, nfmt: used memory cache */
/*     int max_unpack; // unpack flag */
/*     fmt_t *fmt; */
/*     char *format_str; */
/* }; */


/* struct _mval_t { */
/*     int type; // tag types */
/*     int m, n;    // m: max memory cache, n: used */
/*     kstring_t *a;  */
/* }; */

/* typedef struct _mval_t mval_t; */

struct _matrix_cache {
    int allele_count;
    int sample_count;
    int mmtx, nmtx;
    struct matrix_value **matrix;
};

struct args {
    int skip_reference;
    int print_header;
    int split_flag;
    ccols_t *convert;
    mcache_t *cache;
    kstring_t *memory_pool;
    bcf_header_t *hdr;
};

struct args args = {
    .skip_reference = 0,
    .print_header = 1,
    .split_flag = SPLIT_DEFAULT,
    .convert = NULL,
    .cache = NULL,
    .memory_pool = NULL,
    .hdr = NULL,
};

/* struct _tags { */
/*     int n, m; // n : fields , m : genotypesn*samples */
    
/*     int allele_count; */
/*     int k, l; // the sizes of matrix */
/*     mval_t **trans; */
/* }; */

/* // tags format */
/* static tags_t *tag = NULL; */

/* static int skip_ref = 0; */
/* static int print_header = 0; */
void format_string_init(char *s)
{
    if (s== NULL)
	error("Empty format string!");
    
    args.convert = (ccols_t *)calloc(sizeof(ccols_t));
    ccols_t *convert = args.convert;
    convert->format_string = strdup(s);
    convert->
}
enum col_type parse_key (char *p, enum field_type type)
{
    if (p == NULL) return is_unknown;
    char *q = p;

    if (*q == '%') p++;
    else return is_sep;

    if (*p
#define same_string(a, b) (!strcmp(a, b))
    
    if (type == field_format) {
	if (same_string(p, "GT") || same_string(p, "TGT")) return is_gt;
	else if (same_string(p, "SAMPLE")) return is_sample;
	else 
    } else {

    }

#undef same_string
}
int init_args(const char* string)
{
}
int set_tag(tags_t *t)
{
    assert(t);
    
    int i, j, k;
    for (i=0; i<=t->m; ++i) {    
	for (j=0; j<=t->n; ++j) {
	    for (k=0; k < t->trans[i][j].m; ++k) {
		t->trans[i][j].a[k].l = 0;
	    }
	    t->trans[i][j].m = 0;
	}
    }
    t->m = 0;
    t->n = 0;
    return 0;
}

void free_kstring(kstring_t *s)
{
    if (s->m)
	free((s)->s);
}

void destroy_tag(tags_t * t)
{
    int i, j, k;
    for (i=0; i<=t->n; ++i) {
	for (j=0; j<=t->m; ++j) {
	    for (k=0; k < t->trans[i][j].m; ++k) 
		free_kstring(&t->trans[i][j].a[k]);	    	    
	    if (t->trans[i][j].n) free(t->trans[i][j].a);
	}
	free(t->trans[i]);
    }
    free(t->trans);
    free(t);
}
int tags2str(const tags_t* t, kstring_t *s)
{
    int i, j, k, d;
    int l_ori = s->l;
    // allele
    for ( i=0; i<=t->m; ++i ) {
	d=1; k=0;
	// trans
	for (; k<d;k++) {
	    //rows
	    for (j=0; j<=t->n; ++j) {
		
		mval_t *f= &t->trans[i][j];
		fprintf(stderr, "%s\n", f->a[0].s);
		if ( (f->type & IS_TRANS) && f->m > 1) {
		    if (d != 1 && d != f->m)
			error("TRANS tags have different transcripts!\n");
		    d=f->m;
		}
		if (f->type & IS_TRANS) { kputs(f->a[k].s, s); }
		else { kputs(f->a[0].s, s); }
	    } // end rows
	    //   kputc('\n', s);
	} // end trans
    } // end allele
    return s->l -l_ori;
}

static convert_t *convert = NULL;

// convert the header of output
int convert_header(kstring_t *str)
{
    int l_ori = str->l;
    int i;
    kputc('#', str);
    for ( i=0; i<convert->nfmt; ++i ) {
	// is FORMAT field
	if ( convert->fmt[i].is_format == 1) {

	    if ( !(split_flag & SPLIT_SAMPLE )) {
		int j = i, js, k;
		while ( j<convert->nfmt && convert->fmt[j].is_gtf ) { j++; }
		for ( js = 0; js < bcf_hdr_nsamples(header); ++js ) {
		    
		    for ( k=i; k<j; ++k ) {
			
			if ( convert->fmt[k].type & IS_SEP ) {
			    if ( convert->fmt[k].key )
				kputs( convert->fmt[k].key, str);
			} else {
			    ksprintf(str, "%s:%s", header->samples[js], convert->fmt[k].key);
			}
		    }
		}
		i = j -1; // skip sample fields
	    } else {
		if ( convert->fmt[i].type & IS_SAM ) {
		    if ( !(split_flag & SPLIT_SAMPLE) )
			error("split tag has inconsistent format\n");
		    kputs("SAMPLE", str);
		} else {
		    kputs(convert->fmt[i].key, str);
		}
	    }
	} else {
	    // Fixed fields or INFO fields
	    if ( convert->fmt[i].key ) {	    
		if (!strcmp(convert->fmt[i].key, "BED")) { kputs("CHROM\tSTART\tSTOP", str); }
		else { kputs(convert->fmt[i].key, str); }
	    } else {
		// empty fields
		assert(1);
	    }
	}
    }
    return str->l - l_ori;
}
void init_tags(tags_t *t, const bcf1_t* line)
{
    int i, j;
    
    int mrow = split_flag & SPLIT_ALT ? line->n_allele : 1;

    if ( split_flag & SPLIT_SAMPLE )
	mrow *= header->n[BCF_DT_SAMPLE];
    if (t->k == 0)
	t->k = convert->nfmt;
    
    if ( t->l < mrow ) {

	t->trans = (mval_t**)realloc(t->trans, mrow*sizeof(mval_t*)); 
	
	for (i=t->l; i <mrow; ++i) {

	    t->trans[i] = (mval_t*)calloc(t->k, sizeof(mval_t));
	    for (j=0; j<t->k; ++j) {
		t->trans[i][j].m = 0;
		t->trans[i][j].n = 1;
		t->trans[i][j].type = TAG_UNKNOWN;
		t->trans[i][j].a = (kstring_t*)calloc(t->trans[i][j].n, sizeof(kstring_t));
	    }
	}
	t->l = mrow;
    }

    for (i=0; i<convert->nfmt; ++i) {
    	//if (convert->fmt[i].id == -1) continue;
	if ( convert->fmt[i].type == TAG_FORMAT) {
	    for (j=0; j<(int)line->n_fmt; j++) {
		if ( line->d.fmt[j].id==convert->fmt[i].id ) { convert->fmt[i].fmt = &line->d.fmt[j]; break; }
	    }
	}
    }
    t->m = t->n = -1;
}
int convert_line(bcf1_t *line, kstring_t *str)
{
    int l_ori = str->l;
    bcf_unpack(line, convert->max_unpack);
    int i, k=0;

    init_tags(tag, line);
    int isample=-1, nsample=0;
    int row = 0;

    if ( split_flag & SPLIT_SAMPLE ) {    
	isample = 0;
	nsample = header->n[BCF_DT_SAMPLE];
    }
    
    for ( ; isample<nsample; ++isample )
    {
	// allele iter
	int i_ala = NONFLAG;
	int l_ala = i_ala;

	if ( split_flag & SPLIT_ALT ) {

	    bcf_fmt_t *fgt = bcf_get_fmt(header, line, "GT");
	    
#define BRANCH(type_t, vector_end) do {					\
		type_t *ptr = (type_t*) (fgt->p + isample*fgt->size);	\
		i_ala = (ptr[k]>>1)-1;					\
	    } while(0)
	    
	    for ( k=0; k<fgt->n; ++k) {

		switch (fgt->type) {
		    case BCF_BT_INT8:
			BRANCH(int8_t, bcf_int8_vector_end);
			break;
			
		    case BCF_BT_INT16:
			BRANCH(int16_t, bcf_int16_vector_end);
			break;
			
		    case BCF_BT_INT32:
			BRANCH(int32_t, bcf_int32_vector_end);
			break;
			
		    default:
			error("FIXME: type %d in bcf_format_gt?\n", fgt->type);
		}
		if (i_ala == -1) continue;
		if (i_ala != NONFLAG) {
		    if (l_ala == i_ala) continue;
		    l_ala = i_ala;
		}
		if ( skip_ref && i_ala==0) continue;
		tag->m = row++;
		for ( i=0; i<convert->nfmt; i++ ) {
		    tag->n = i;
		    if ( convert->fmt[i].handler )
			convert->fmt[i].handler(line, &convert->fmt[i], i_ala, isample, tag);
		}
	    }
#undef BRANCH
	} else {
	    tag->m = row++;	    
	    for ( i=0; i<convert->nfmt; i++ ) {
		tag->n = i;
		if ( convert->fmt[i].handler )
		    convert->fmt[i].handler(line, &convert->fmt[i], i_ala, isample, tag);
	    }
	}
    }
    if (tag->m != -1 && tag->n != -1)
	tags2str(tag, str);
    set_tag(tag);
    return str->l - l_ori;
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
	else if (strcmp(list[i], "SAMPLE") == 0) { flag |= SPLIT_SAMPLE; }
	else if (strcmp(list[i], "TRANS") == 0) { flag |= (SPLIT_TRANS|SPLIT_ALT); }
	else if (strcmp(list[i], "ALL") == 0) { flag |= SPLIT_ALL; }
	else {
	    fprintf(stderr, "cannot recongize this tag : %s\n", list[i]);
	    exit(1);
	}
    }
    free_list(list, n);
    return flag;
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

/* The register_tag and parse_tag functions are adapted from pd3's convert.c 
 * use TAGS-tags instead of T-tags 
 */
fmt_t *register_tag(int type, char *key, int is_gtf)
{
    convert->nfmt++;
    if (convert->nfmt == convert->mfmt) {
	convert->mfmt += 10;
	convert->fmt = (fmt_t*)realloc(convert->fmt, convert->mfmt*sizeof(fmt_t));
    }
    fmt_t *fmt = &convert->fmt[ convert->nfmt-1 ];
    fmt->type = type;
    fmt->is_gtf = is_gtf;
    fmt->key = key ? strdup(key) : NULL;
    // allow non-format tags to appear amongst the format fields, split-samples mode
    if ( key ) {
	fmt->id = bcf_hdr_id2int(header, BCF_DT_ID, fmt->key);
	if (fmt->type==TAG_FORMAT && !bcf_hdr_idinfo_exists(header, BCF_HL_FMT, fmt->id)) {

	    if ( !strcmp("CHROM", key) ) fmt->type = TAG_CHROM;
	    else if ( !strcmp("POS", key) ) fmt->type = TAG_POS;
	    else if ( !strcmp("BED", key) ) fmt->type = TAG_BED;
	    else if ( !strcmp("ID", key) ) fmt->type = TAG_ID;
	    else if ( !strcmp("REF", key) ) fmt->type = TAG_REF;
	    else if ( !strcmp("ALT", key) ) fmt->type = TAG_ALT;
	    //else if ( !strcmp("FIRST_ALT", key) ) fmt->type = TAG_FIRST_ALT;
	    else if ( !strcmp("QUAL", key) ) fmt->type = TAG_QUAL;
	    else if ( !strcmp("FILTER", key) ) fmt->type = TAG_FILTER;
	    else if ( !strcmp("SAMPLE", key) ) fmt->type = TAG_SAMPLE;
	    //else if ( !strcmp("HGVS", key) ) fmt->type = TAG_TRANS;
	    else if ( fmt->id>=0 && bcf_hdr_idinfo_exists(header, BCF_HL_INFO, fmt->id)) fmt->type = TAG_INFO;
	}
    }

    switch ( fmt->type ) {
    
	/* case TAG_FIRST_ALT: */
	/*     fmt->handler = &process_first_alt; */
	/*     break; */
	    
	case TAG_CHROM:
	    fmt->handler = &process_chrom;
	    break;
	    
	case TAG_POS:
	    fmt->handler = &process_pos;
	    break;
	    
	case TAG_BED:
	    fmt->handler = &process_bed;
	    break;

	case TAG_ID:
	    fmt->handler = &process_id;
	    break;

	case TAG_REF:
	    fmt->handler = &process_ref;
	    break;

	case TAG_ALT:
	    fmt->handler = &process_alt;
	    break;

	case TAG_QUAL:
	    fmt->handler = &process_qual;
	    break;
	    
	case TAG_FILTER:
	    fmt->handler = &process_filter;
	    convert->max_unpack |= BCF_UN_FLT;
	    break;
	    
	case TAG_INFO:
	    fmt->handler = &process_info;
	    convert->max_unpack |= BCF_UN_INFO;
	    break;
	    
	/* case TAG_TRANS: */
	/*     fmt->handler = &process_trans; */
	/*     convert->max_unpack |= BCF_UN_INFO; */
	/*     break; */
	    
	case TAG_FORMAT:
	    fmt->handler = &process_format;
	    convert->max_unpack |= BCF_UN_FMT;
	    break;
	    
	case TAG_SAMPLE:
	    fmt->handler = &process_sample;
	    break;
	    
	case TAG_SEP:
	    fmt->handler = &process_sep;
	    break;
	    
	case TAG_TGT:
	    fmt->handler = &process_tgt;
	    convert->max_unpack |= BCF_UN_FMT;
	    break;
	    
	default:
	    error("TODO: handler for type %d\n", fmt->type);
    }
    if ( key ) {
	if ( fmt->type==TAG_INFO ) {
	    fmt->id = bcf_hdr_id2int(header, BCF_DT_ID, key);
	    if ( fmt->id==-1 ) error("Error: no such tag defined in the VCF header: INFO/%s\n", key);
	}
    }
    return fmt;
}

static char *parse_tag(char *p, int is_gtf)
{
    char *q = ++p;
    while ( *q && (isalnum(*q) || *q=='_' || *q=='.')) q++;
    kstring_t str = {0,0,0};
    if ( q-p==0 ) error("Could not parse format string: %s\n", convert->format_str);
    kputsn(p, q-p, &str);

    if ( is_gtf ) {
	if ( !strcmp(str.s, "GT") || !strcmp(str.s, "TGT")) {
	    register_tag(TAG_TGT, "GT", is_gtf);
	}
	else if ( !strcmp(str.s, "SAMPLE") ) {
	    register_tag(TAG_SAMPLE, "SAMPLE", is_gtf); // is_gtf ==> 0
	}
	//else if ( !strcmp(str.s, "IUPACGT") ) register_tag(TAG_IUPAC_GT, "GT", is_gtf);
	else {
	    register_tag(TAG_FORMAT, str.s, is_gtf);
	}
    } else {
	if ( !strcmp(str.s, "CHROM") ) {
	    register_tag(TAG_CHROM, "CHROM", is_gtf);
	}
	else if ( !strcmp(str.s, "POS") ) {
	    register_tag(TAG_POS, "POS", is_gtf);
	}
	else if ( !strcmp(str.s, "BED") ) {
	    register_tag(TAG_BED, "BED", is_gtf);
	}
	else if ( !strcmp(str.s, "POS") ) {
	    register_tag(TAG_CHROM, "POS", is_gtf);
	}
	else if ( !strcmp(str.s, "REF") ) {
	    register_tag(TAG_REF, "REF", is_gtf);
	}
	else if ( !strcmp(str.s, "ALT") ) {
	    register_tag(TAG_ALT, "ALT", is_gtf);
	}
	//else if ( !strcmp(str.s, "HGVS") ) register_tag(TAG_TRANS, "HGVS", is_gtf);
	else if ( !strcmp(str.s, "GT") || !strcmp(str.s, "TGT")) {
	    split_flag |= SPLIT_SAMPLE;
	    register_tag(TAG_TGT, "GT", is_gtf);
	}
	//else if ( !strcmp(str.s, "IUPACGT") ) { split_flag |= SPLIT_SAMPLE; register_tag(TAG_IUPAC_GT, "GT", is_gtf); }
	//else if ( !strcmp(str.s, "FUNC") ) register_tag(TAG_TRANS, "FUNC", is_gtf);
	else if ( !strcmp(str.s, "SAMPLE") ) {
	    split_flag |= SPLIT_SAMPLE;
	    register_tag(TAG_SAMPLE, "SAMPLE", is_gtf);
	}
	//else if ( !strcmp(str.s, "FIRST_ALT") ) register_tag(TAG_FIRST_ALT, "FIRST_ALT", is_gtf);
	else if ( !strcmp(str.s, "QUAL") ) {
	    register_tag(TAG_QUAL, "QUAL", is_gtf);
	}
	else if ( !strcmp(str.s, "INFO") ) {
	    if ( *q=='/' ) error("Could not parse format string: %s\n", convert->format_str);
	    p = ++q;
	    str.l = 0;
	    while ( *q && (isalnum(*q) || *q=='_' || *q=='.') ) q++;
	    if ( q-p==0 ) error("Could not parse format string: %s\n", convert->format_str);
	    kputsn(p, q-p, &str);
	    register_tag(TAG_INFO, str.s, is_gtf);
	} else {
	    register_tag(TAG_INFO, str.s, is_gtf);
	}
    }
    free(str.s);
    return q;
}
static char *parse_sep(char *p, int is_gtf)
{
    char *q = p;
    kstring_t str = { 0, 0, 0};
    while ( *q && *q!='[' && *q!=']' && *q!='%' ) {
	if ( *q=='\\' ) {
	    q++;
	    if ( *q=='n' ) kputc('\n', &str);
	    else if ( *q == 't') kputc('\t', &str);
	    else kputc(*q, &str);
	} else {
	    kputc(*q, &str);
	}
	q++;
    }
    if ( !str.l )
	error("Could not parse format string: %s\n", convert->format_str);
    
    register_tag(TAG_SEP, str.s, is_gtf);
    free(str.s);
    return q;
}
void convert_init(char *s)
{
    convert = (convert_t*)malloc(sizeof(convert_t));
    convert->format_str = strdup(s);
    convert->nfmt = 0;
    convert->mfmt = 2;
    convert->fmt = (fmt_t*)calloc(convert->mfmt, sizeof(fmt_t)); 
    convert->max_unpack = 0;
    int is_gtf = 0;
    char *p = convert->format_str;
    while ( *p ) {
	switch (*p) {
	case '[':
	    is_gtf = 1;
	    p++;
	    break;
	    
	case ']':
	    is_gtf = 0;
	    p++;
	    break;
	    
	case '%':
	    p = parse_tag(p, is_gtf);
	    break;
	    
	default:
	    p = parse_sep(p, is_gtf);
	    break;
	}
    }
}
const char *about(void)
{
    return "Select tags in pre-defined format.\n";
}

int usage(void)
{
    fprintf(stderr,"About : Select tags from VCF/BCF file.\n");
    fprintf(stderr,"Usage:\n");
    fprintf(stderr,"Standalone mode:\n");
    fprintf(stderr,"\tbcfselect [Options] in.vcf.gz\n");
    fprintf(stderr,"Options:\n");
    fprintf(stderr,"\t-f, --format   see man page for deatils.\n");
    fprintf(stderr,"\t-s, --split    split by [ALT,SAMPLE,TRANS].\n");
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
		skip_ref = 1;
		break;
		
	    case 's':
		flag = strdup(optarg);
		break;
		
	    case 'p':
		print_header = 1;
		break;
		
	    case 'h':
	    case '?':
	    default:
		return usage();

	}
    }
    if (format == NULL)
	error("-f is required by bcfselect.");

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
    
    if (format) {
	convert_init(format);
	free(format);
    } else {
	error("Must set a format string with -f \n");
    }
    
    mempool = (kstring_t*)malloc(sizeof(kstring_t));
    mempool->m = mempool->l = 0;
    mempool->s = NULL;
    
    tag = (tags_t*)malloc(sizeof(tags_t));
    tag->l = tag->k = 0;

    if (print_header)
	convert_header(mempool);

    while (bcf_sr_next_line(sr)) {

	bcf1_t *line = bcf_sr_get_line(sr, 0);
	convert_line(line, mempool);
	if (mempool->l > MEMPOOL) {
	    fprintf(stdout, "%s", mempool->s);
	    mempool->l = 0;
	}
    }
    if (mempool->l) {
	fprintf(stdout, "%s", mempool->s);
	mempool->l = 0;
    }
    bcf_sr_destroy(sr);
    destroy_tag(tag);
    if (mempool->m) free(mempool->s);
    free(mempool);
    return 0;
}



int main(int argc, char **argv)
{
    return run(argc, argv);
}

