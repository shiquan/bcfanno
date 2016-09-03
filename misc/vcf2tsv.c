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
    _pormote_to_int = -1, // use a minus value to promote every enum type to a int type
    is_unknown = 0,
    is_chrom,
    is_pos,
    is_id,
    is_filter,
    is_ref,
    is_alt,
    is_info,
    is_gt,
    is_format,
    is_sample,
    is_bed,
};

/*
  SAMPLE1  ALLELE 1
  SAMPLE1  ALLELE 2
  SAMPLE2  ALLELE 1
  SAMPLE2  ALLELE 2
 */

typedef struct _col col_t;
typedef struct _multi_cols_cache_per_sample mcache_ps_t;
typedef struct _multi_cols_cache_per_allele mcache_pa_t;
typedef struct _multi_cols_cache mcache_t;
typedef struct _convert_cols ccols_t;
typedef struct _mval mval_t;

struct _col {
    enum col_type type;
    int number; // BCF_VL_*
    int unpack;
    int id; //header index id for INFO/FORMAT
    char *key;
    void (*setter)(bcf_hdr_t *hdr, bcf1_t *, col_t *, int ale,  mval_t *);
};

struct _convert_cols {
    int m, l;
    col_t **cols;
    int max_unpack;
};

/* value for each node in print cached matrix */
struct _mval {
    enum col_type type;
    int sample_id;
    kstring_t a;
};

struct _multi_cols_cache_per_allele {
    int n_cols;
    struct _mval *mvals;
};
// matrix cache
struct _multi_cols_cache_per_sample {
    int n_alleles;
    int m_alleles;
    struct _multi_cols_cache_per_allele * alvals;
};

struct _multi_cols_cache {
    int m_samples;
    struct _multi_cols_cache_per_sample *mcols;
};

/* declared setter functions */
void setter_gt(bcf_hdr_t *, bcf1_t *, col_t *, int , mval_t*);
void setter_sample(bcf_hdr_t *, bcf1_t *, col_t *, int , mval_t*);
void setter_pos(bcf_hdr_t *, bcf1_t *, col_t *, int , mval_t*);
void setter_ref(bcf_hdr_t *, bcf1_t *, col_t *, int , mval_t*);
void setter_alt(bcf_hdr_t *, bcf1_t *, col_t *, int , mval_t*);
void setter_chrom(bcf_hdr_t *, bcf1_t *, col_t *, int , mval_t*);
void setter_id(bcf_hdr_t *, bcf1_t *, col_t *, int , mval_t*);
void setter_filter(bcf_hdr_t *, bcf1_t *, col_t *, int , mval_t*);
void setter_bed(bcf_hdr_t *, bcf1_t *, col_t *, int , mval_t*);
void setter_format(bcf_hdr_t *, bcf1_t *, col_t *, int , mval_t*);
void setter_info(bcf_hdr_t *, bcf1_t *, col_t *, int , mval_t*);

// options and handlers
struct args {
    int skip_ref;
    int print_header;
    int split_flag;
    int skip_reference;
    ccols_t *convert;
    mcache_t *cache;
    kstring_t *mempool;
    bcf_hdr_t *hdr;
};

// init args
struct args args = {
    .skip_ref = 0,
    .print_header = 1,
    .skip_reference = 0,
    .split_flag = SPLIT_DEFAULT,
    .convert = NULL,
    .cache = NULL,
    .mempool = NULL,
    .hdr = NULL,
};

ccols_t *ccols_init()
{
    ccols_t *c = (ccols_t *)malloc(sizeof(ccols_t));
    c->m = c->l = 0;
    c->cols = 0;
    c->max_unpack = 0;
    return c;
}

col_t *register_key (char *p, bcf_hdr_t *h);

ccols_t *format_string_init(char *s, bcf_hdr_t *h, int *n_sample)
{
    if (s== NULL)
	error("Empty format string!");

    ccols_t *c = ccols_init();
    int has_sample = 0, has_format = 0;
    char *p = s;
    int proc = 1;
    while (proc && *p ) {
	char *q = p;	
	while (*q && *q != ',') q++;
	if (q-p != 0) {	
	    if (*q == '\0') proc = 0;
	    else *q = '\0';
	    if (c->l == c->m) {
		c->m = c->l == 0 ? 2 : c->m << 1;
		c->cols = (col_t **)realloc(c->cols, sizeof(col_t*)*c->m);
	    }
	    col_t *t = register_key(p, h);
	    if(t == NULL)
		error("failed to registe tag %s", p);
	
	    if (t->type == is_sample) has_sample = 1;
	    else if (t->type == is_format) has_format = 1;

	    c->max_unpack |= t->unpack;
	    c->cols[c->l++] = t;
	}
	p = q+1;
    }


    if ( has_sample ) {
      *n_sample = bcf_hdr_nsamples(h);
    } else {
      if (has_format)
	warnings("no SAMPLE tag in the format string, only output first sample %s", h->samples[0]);
      *n_sample = 1;
    }

    /* int i; */
    /* for (i=0; i<c->l; ++i) */
    /* 	debug_print("key : %s", c->cols[i]->key); */
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
    
    c->key = strdup(p);
    c->id = -1;
#define same_string(a, b) (!strcmp(a, b))
    if (same_string(q, "GT") || same_string(q, "TGT")) {
	c->setter = setter_gt;
	c->type = c->type == is_unknown || c->type == is_format ? is_gt : c->type;
	c->unpack |= BCF_UN_FMT;
	c->id = bcf_hdr_id2int(h, BCF_DT_ID, "GT");
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
    else if (same_string(q, "REF")) {
	c->setter = setter_ref;
	c->type = is_ref;
    }
    else if (same_string(q, "ALT")) {
	c->setter = setter_alt;
	c->type = is_alt;
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
	    c->unpack |= BCF_UN_FMT;
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
    m->mcols = (mcache_ps_t *)calloc(m->m_samples, sizeof(mcache_ps_t));
    int i;
    for (i=0; i<m->m_samples; ++i) {
	mcache_ps_t *ps = &m->mcols[i];
	ps->n_alleles = ps->m_alleles = 0;
	ps->alvals = 0;
    }
    return m;
}

void set_matrix_cache(mcache_t *m, int n_alleles)
{
    int i, j, k;
    for (i =0; i< m->m_samples; ++i) {
	// enlarger memory cache
	mcache_ps_t *ps = &m->mcols[i];
	if (ps->m_alleles <= n_alleles) {
	    ps->m_alleles = n_alleles;
	    ps->alvals = (mcache_pa_t*)realloc(ps->alvals, ps->m_alleles *sizeof(mcache_pa_t));
	    for (j=ps->n_alleles; j< ps->m_alleles; ++j) {
		mcache_pa_t *pa = &ps->alvals[j];
		pa->n_cols = args.convert->l;
		pa->mvals = (mval_t*)calloc(pa->n_cols, sizeof(mval_t));
		for (k=0; k<pa->n_cols; ++k) {
		    pa->mvals[k].a.l = pa->mvals[k].a.m = 0;
		    pa->mvals[k].a.s = 0; 
		}
	    }
	}
	
	ps->n_alleles = n_alleles;
	for (j=0; j<ps->n_alleles; ++j) {
	    mcache_pa_t *pa = &ps->alvals[j];
	    for (k=0; k<pa->n_cols; ++k)
		pa->mvals[k].a.l = 0;	
	}
    }
}

void release_mcache(mcache_t *m)
{
    int i, j, k;
    for (i =0; i< m->m_samples; ++i) {
	mcache_ps_t *ps = &m->mcols[i];
	for (j=0; j<ps->n_alleles; ++j) {
	    mcache_pa_t *pa = &ps->alvals[j];
	    for (k=0; k<pa->n_cols; ++k) {
		if (pa->mvals[k].a.m)
		    free(pa->mvals[k].a.s);		
	    }
	    free(pa->mvals);
	}
	free(ps->alvals);
    }
    free(m->mcols);
    free(m);
}

void release_args()
{
    int i;
    for (i=0; i<args.convert->l; ++i) {
	col_t *c = args.convert->cols[i];
	free(c->key);
	free(c);
    }
    free(args.convert->cols);
    free(args.convert);
    release_mcache(args.cache);

    if (args.mempool->m)
	free(args.mempool->s);
    args.mempool->l = args.mempool->m = 0;
    args.mempool->s = 0;
}

// convert the header of output
int convert_header()
{
    ccols_t *cols = args.convert;
    int i;
    for (i=0; i<cols->l; ++i) {
	if (i) kputc('\t', args.mempool);
	else kputc('#', args.mempool);
	col_t *c = cols->cols[i];
	switch(c->type) {
	    case is_bed :
		kputs("CHROM\tSTART\tEND", args.mempool);
		break;

	    case is_unknown:
		error("unknow type!");
		
	    case is_chrom:
	    case is_pos:
	    case is_ref:
	    case is_alt:
	    case is_id:
	    case is_filter:
	    case is_info:
	    case is_gt:
	    case is_format:
	    case is_sample:
	    default:
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
    int flag = SPLIT_NONE;
    char **list = hts_readlist(s, 0, &n);
    for (i = 0; i < n; ++i) {
	if (strcmp(list[i], "NONE") == 0) { free_list(list, n); return SPLIT_NONE; }
	else if (strcmp(list[i], "ALT") == 0) { flag |= SPLIT_ALT; }
	else {
	    fprintf(stderr, "cannot recongize this tag : %s\n", list[i]);
	    exit(1);
	}
    }
    free_list(list, n);
    return flag;
}

int convert_line(bcf_hdr_t *hdr, bcf1_t *line)
{
    if (line == NULL)
	error ("null line.");
    bcf_unpack(line, args.convert->max_unpack);

    mcache_t *cache = args.cache;
    // ccols_t *cols = args.convert;
    int n_alleles = args.split_flag & SPLIT_ALT ? line->n_allele : 1;
    set_matrix_cache(cache, n_alleles);

    int i, j,k;

    for (i=0; i<cache->m_samples; ++i) { 
	// iterate samples
	mcache_ps_t *ps = &cache->mcols[i];
	//ps->n_alleles = n_alleles;
	for (j=0; j < ps->n_alleles; ++j) { // ???
	    
	    if (n_alleles > 1 && args.skip_ref && j==0) continue;
	    mcache_pa_t *pa = &ps->alvals[j];
	    int iallele = -1;
	    if (args.split_flag & SPLIT_ALT) iallele = j;
	    for (k = 0; k < pa->n_cols; ++k) {		
		col_t *col = args.convert->cols[k];
		mval_t *val = &pa->mvals[k];		
		val->type = col->type;
		val->sample_id = i;		
		/* setter function */
		col->setter(hdr, line, col, iallele, val);
		if (k) kputc('\t', args.mempool);
		kputs(val->a.s, args.mempool);
		//debug_print("%s", val->a.s);
	    } // end cols
	    kputc('\n', args.mempool);
	} // end alleles
    }  // end samples
    return 0;
}

void setter_sample(bcf_hdr_t *hdr, bcf1_t *line, col_t *c, int ale, mval_t *val)
{
    kputs(hdr->samples[val->sample_id], &val->a);
}
void setter_chrom(bcf_hdr_t *hdr, bcf1_t *line, col_t *c, int ale, mval_t *val)
{
    kputs(hdr->id[BCF_DT_CTG][line->rid].key, &val->a);
}
void setter_pos(bcf_hdr_t *hdr, bcf1_t *line, col_t *c, int ale, mval_t *val)
{
    kputw(line->pos +1, &val->a);
}
void setter_bed(bcf_hdr_t *hdr, bcf1_t *line, col_t *c, int ale, mval_t *val)
{
    uint32_t end = line->pos + strlen(line->d.allele[0]);
    ksprintf(&val->a, "%s\t%u\t%u",hdr->id[BCF_DT_CTG][line->rid].key, line->pos, end);
}
void setter_ref(bcf_hdr_t *hdr, bcf1_t *line, col_t *c, int ale, mval_t *val)
{
    kputs(line->d.allele[0], &val->a);
}
void setter_alt(bcf_hdr_t *hdr, bcf1_t *line, col_t *c, int ale, mval_t *val)
{
    if (line->n_allele == 1) {
	kputc('.', &val->a);
    } else {
	if (ale == -1) {
	    int i;
	    for (i=1; i<line->n_allele; ++i) {
		if (i>1) kputc(',', &val->a);
		kputs(line->d.allele[i], &val->a);		
	    }
	} else {
	    assert(ale < line->n_allele);
	    kputs(line->d.allele[ale], &val->a);
	}

    }
}
void setter_id(bcf_hdr_t *hdr, bcf1_t *line, col_t *c, int ale, mval_t *val)
{
    kputs(line->d.id, &val->a);
}
void setter_qual(bcf_hdr_t *hdr, bcf1_t *line, col_t *c, int ale, mval_t *val)
{
    if (bcf_float_is_missing(line->qual)) { kputc('.', &val->a); }
    else { ksprintf(&val->a, "%g", line->qual); }
}
void setter_filter(bcf_hdr_t *hdr, bcf1_t *line, col_t *c, int ale, mval_t *val)
{
    int i;
    if (line->d.n_flt) {
	for (i=0; i<line->d.n_flt; ++i) {
	    if (i) kputc(';', &val->a);
	    kputs(hdr->id[BCF_DT_ID][line->d.flt[i]].key, &val->a);
	}
    } else {
	kputc('.', &val->a);
    }
}
// example : chr pos A T
// hom-ref for homozygous reference allele, A/A
// ref-alt for heterozygosity, A/T
// hom-alt for homozygous alternatve allele, T/T
// multi-alt for multi-alternative alleles, T/G
// uncov-ref for uncovered allele and reference allele, ./A
// uncov-alt for uncovered allele and alternative allele, ./T
// uncov for uncovered, ./.
void setter_zygosity(bcf_hdr_t *hdr, bcf1_t *line, col_t *c, int ale, mval_t *val)
{
    bcf_fmt_t *fmt = bcf_get_fmt_id(line, c->id);
    if ( fmt == NULL )
	error("no found TG tag in line : %s,%d", hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1);
    int sample_id = val->sample_id;

#define BRANCH(type_t, missing, vector_end) do {\
	type_t *ptr = (type_t*)(fmt->p + sample_id*fmt->size);\
	int i;\
	int allele1 = -2; \
	int allele2 = -2;\
	assert(fmt->n <= 2);\
	allele1 = (ptr[0]>>1)-1;\
	if (fmt->n == 2) \
	    allele2 = (ptr[1]>>1) -1;\
	if (allele1 ==  
	    

}
void setter_gt(bcf_hdr_t *hdr, bcf1_t *line, col_t *c, int ale, mval_t *val)
{
    if (ale != -1)
	error ("TG, TGT only used with split-allele mode.");

    bcf_fmt_t *fmt = bcf_get_fmt_id(line, c->id);

    if (fmt == NULL)
	error ("no found TG tag in line : %s,%d", hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1);

    int sample_id = val->sample_id;
    
#define BRANCH(type_t, missing, vector_end) do {		\
	type_t *ptr = (type_t*)(fmt->p + sample_id*fmt->size);\
	int i;\
	for (i=0; i<fmt->n; ++i) {\
	    if ( i ) kputc("/|"[ptr[i]&1], &val->a);\
	    if ( !(ptr[i]>>1) ) kputc('.', &val->a); \
	    else kputs(line->d.allele[(ptr[i]>>1)-1], &val->a);	\
	}\
} while(0)

      switch(fmt->type) {
	  case BCF_BT_INT8:
	      BRANCH(int8_t, bcf_int8_missing, bcf_int8_vector_end);
	      break;

	  case BCF_BT_INT16:
	      BRANCH(int16_t, bcf_int16_missing, bcf_int16_vector_end);
	      break;

	  case BCF_BT_INT32:
	      BRANCH(int32_t, bcf_int32_missing, bcf_int32_vector_end);
	      break;

	  default:
	      error("FIXME: type %d in bcf_format_gt?", fmt->type);
      }
#undef BRANCH
}

void process_fmt_array(int iallele, kstring_t *string, int n, int type, void *data)
{

#define BRANCH(type_t, is_missing, is_vector_end, string)  do  {	\
	type_t *p = (type_t*)data;					\
	if (iallele == -1) {						\
	    int i;							\
	    for (i=0; i<n; ++i) {					\
		if (p[i] == is_vector_end) break;			\
		if (i) kputc(',', string);				\
		if (p[i] == is_missing) kputc('.', string);		\
		else kputw(p[i], string);			\
	    }								\
	} else {\
	    if (p[iallele] == is_vector_end || p[iallele] == is_missing) kputc('.', string); \
	    else kputw(p[iallele], string);				\
	}								\
} while(0)

      switch(type) {
	  case BCF_BT_CHAR:
	      do {
		  char *p = (char*)data;
		  int i;
		  if (iallele == -1) {
		      for (i=0; i<n && *p; ++i,++p) {
			  if (*p == bcf_str_missing) kputc('.', string);
			  else kputc(*p, string);
		      }
		  } else {
		      //p += iallele;
		      for (i=0; i<n && *p; ++i,++p) {
			  if (*p == bcf_str_missing) kputc('.', string);
			  else kputc(*p, string);
		      }
		      /* if (*p) kputc(*p, string); */
		      /* else kputc('.', string); */
		  }
	      } while(0);
	      break;
		
	  case BCF_BT_INT8:
	      BRANCH(int8_t, bcf_int8_missing, bcf_int8_vector_end, string);
	      break;

	  case BCF_BT_INT16:
	      BRANCH(int16_t, bcf_int16_missing, bcf_int16_vector_end, string);
	      break;

	  case BCF_BT_INT32:
	      BRANCH(int32_t, bcf_int32_missing, bcf_int32_vector_end, string);
	      break;

	  case BCF_BT_FLOAT:
	      do {
		  float *p = (float*)data;
		  int i = 0;
		  if (iallele == -1) {
		      for (i=0; i<n; ++i) {
			  if (p[i] == bcf_float_vector_end) break;
			  if (i) kputc(',', string);
			  if (p[i] == bcf_float_missing) kputc('.', string);
			  else ksprintf(string, "%g", p[i]);
		      }
		  } else {
		      assert(iallele < n);
		      if (p[iallele] == bcf_float_vector_end || p[iallele] == bcf_float_missing) kputc('.', string);
		      else ksprintf(string, "%g", p[i]);
		  }
	      } while(0);	      
	      break;

	  default:
	      error("todo: type %d", type);
      }
#undef BRANCH

}
void setter_info(bcf_hdr_t *hdr, bcf1_t *line, col_t *c, int ale, mval_t *val)
{
    assert(c->id>0);
    bcf_info_t *inf = bcf_get_info_id(line, c->id);
    if (inf == NULL) {
	kputc('.', &val->a);
	return;
    }
    if (ale == 0 && c->type == BCF_VL_A) {
	kputc('.', &val->a);
	return;
    }
    // flag tag
    if (inf->len <=0) {
	// kputc('1', &val->a);
	kputs("Yes", &val->a);
	return;
    }

    if ( inf->len == 1) {
	switch (inf->type) {
	    case BCF_BT_INT8:
		if (inf->v1.i == bcf_int8_missing ) kputc('.', &val->a);
		else kputw(inf->v1.i, &val->a);
		break;

	    case BCF_BT_INT16:
		if (inf->v1.i == bcf_int16_missing ) kputc('.', &val->a);
		else kputw(inf->v1.i, &val->a);
		break;

	    case BCF_BT_INT32:
		if (inf->v1.i == bcf_int32_missing ) kputc('.', &val->a);
		else kputw(inf->v1.i, &val->a);
		break;

	    case BCF_BT_FLOAT:
		if (bcf_float_is_missing(inf->v1.f) ) kputc('.', &val->a);
		else ksprintf(&val->a, "%g", inf->v1.f);
		break;

	    case BCF_BT_CHAR:
		if (inf->v1.i == bcf_str_missing ) kputc('.', &val->a);
		else kputc(inf->v1.i, &val->a);
		break;
		
	    default:
		error("todo: type %d\n", inf->type);

	}	
    } else {
	int iallele = ale == -1 || c->type == BCF_VL_G || c->type == BCF_VL_FIXED ? -1 : c->type == BCF_VL_R ? ale -1 : ale;
	process_fmt_array(iallele, &val->a, inf->len, inf->type, inf->vptr);
    }

}
void setter_format(bcf_hdr_t *hdr, bcf1_t *line, col_t *c, int ale, mval_t *val)
{
    bcf_fmt_t *fmt = bcf_get_fmt_id (line, c->id);
    //debug_print("id: %d", c->id);
    if (fmt == NULL) {
	kputc('.', &val->a);
	return;
    }

    if (ale == 0 && c->type == BCF_VL_A) {
	kputc('.', &val->a);
	return;
    }
    // flag tag
    if (fmt->n <= 0) {
	kputc('1', &val->a);
	return;
    }
    int iallele = ale == -1 || c->type == BCF_VL_G || c->type == BCF_VL_FIXED ? -1 : c->type == BCF_VL_R ? ale -1 : ale;
    process_fmt_array(iallele, &val->a, fmt->n, fmt->type, fmt->p + val->sample_id*fmt->size);
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
    fprintf(stderr,"\t-r, --skip-ref     skip format reference positions; suggest open this option.\n");
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

    bcf_hdr_t *header = sr->readers[0].header;
    assert(header);

    if (flag) {
	args.split_flag = init_split_flag(flag);
	free(flag);
    }
    int nsamples = 0;
    args.convert = format_string_init(format, header, &nsamples);
    free(format);

    args.mempool = (kstring_t*)malloc(sizeof(kstring_t));
    args.mempool->m = args.mempool->l = 0;
    args.mempool->s = 0;
    args.cache = mcache_init(nsamples);
    if (args.print_header)
	convert_header();

    while ( bcf_sr_next_line(sr) ) {
	bcf1_t *line = bcf_sr_get_line(sr, 0);
	if (line->rid == -1)
	    continue;
	if ( args.skip_ref == 1 && bcf_get_variant_types(line) == VCF_REF )
	    continue;
	convert_line(header, line);
	if (args.mempool->l) printf("%s",args.mempool->s);
	args.mempool->l = 0;
    }
    if (args.mempool->l) printf("%s",args.mempool->s);
    release_args();
    bcf_sr_destroy(sr);
    return 0;
}

int main(int argc, char **argv)
{
    return run(argc, argv);
}
