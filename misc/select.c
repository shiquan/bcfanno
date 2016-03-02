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
 * bcftools +select -f '%CHROM\t%POS\t%REF[\t%TGT\t%DP]' demo.vcf.gz
 *
 * bcftools +select -f '%BED\t%SAMPLE\t%REF\t%ALT\t[%DP]' demo.vcf.gz
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <ctype.h>
//#include "bcftools.h"

#define SPLIT_NONE    1          //  bcftools query mode, one bcf1_t per line
#define SPLIT_TAG     0
#define SPLIT_SAMPLE  2          //  flag for spliting by sample
#define SPLIT_ALT     4          //  flag for split by allele, each alt per line
#define SPLIT_GT      2          //  flag for split by genotype (0/0, 0/1, 1/1), it's same with split by samples ^ SPLIT_ALT
#define SPLIT_TRANS    (1<<4)    //  flag for split by TRANS name, if use this flag any tags related with transcripts should be split 
#define SPLIT_DEFAULT 0
#define SPLIT_ALL     ( SPLIT_SAMPLE | SPLIT_ALT | SPLIT_TRANS )

#define S_CHROM     1
#define S_POS       2
#define S_ID        4
#define S_REF       8
#define S_ALT       16
#define S_BED       32
#define S_QUAL      64
#define S_FILTER    128
#define S_INFO      256
#define S_FORMAT    513
#define S_SAMPLE    1024
#define S_SEP       2048
#define S_TYPE      4096
#define S_GT        8192
#define S_TGT       8192
#define S_IUPAC_GT  16384
#define S_FIRST_ALT 32768
#define S_TRANS     65536

#define IS_SEP    (S_SEP)
#define IS_SAM    (S_SAMPLE)
#define IS_ALT    (S_BED | S_ALT | S_INFO | S_FORMAT | S_TRANS)
#define IS_GT     (S_GT | S_TGT | S_IUPAC_GT)
#define IS_TRANS  (S_TRANS)

#define MEMPOOL     2621440

#define NONFLAG  -2
static void error(const char *format, ...);

/* split flag */
static int split_flag = SPLIT_DEFAULT;

bcf_hdr_t * header;

/* the fmt_t and convert_t are adapted from pd3's convert.c 
 * if we want export the sample info and allele info per line
 * we must rewrite the process_* functions in convert.c
 */
typedef struct _convert convert_t;
typedef struct _tags     tags_t;

typedef struct _fmt
{
    int id, type, is_gtf;
    char *key;
    bcf_fmt_t *fmt;
    void (*handler)(bcf1_t *, struct _fmt *, int iala, int isample, tags_t *);
}
fmt_t;

struct _convert
{
    int mfmt, nfmt;
    int max_unpack;
    fmt_t *fmt;
    char *format_str;
};

typedef struct
{
    int type;
    int m, n;    // the array size of a[], usually $m==1
    kstring_t *a; 
}
mval_t;

struct _tags
{
    int n, m; // n : fields , m : genotypesn*samples
    int k, l; // the sizes of matrix
    mval_t **trans;
};

static tags_t *tag;
static int skip_ref = 0;
static int print_header = 0;
void clear_tag(tags_t *t)
{
    int i, j, k;
    for (i=0; i<=t->m; ++i)
    {
	for (j=0; j<=t->n; ++j)
	{
	    for (k=0; k < t->trans[i][j].m; ++k)
	    {
		t->trans[i][j].a[k].l = 0;
	    }
	    t->trans[i][j].m = 0;
	}
    }
    t->m = 0;
    t->n = 0;
}
void destroy_tag(tags_t * t)
{
    int i, j, k;
    for (i=0; i<=t->n; ++i)
    {
	for (j=0; j<=t->m; ++j)
	{
	    for (k=0; k < t->trans[i][j].m; ++k)
	    {
		if (t->trans[i][j].a[k].m) free(t->trans[i][j].a[k].s);
	    }
	    if (t->trans[i][j].n) free(t->trans[i][j].a);
	}
	free(t->trans[i]);
    }
    free(t->trans);
    free(t);
}
int tags2str(tags_t * const t, kstring_t *s)
{
    int i, j, k, d;
    int l_ori = s->l;
    for (i=0; i<=t->m; ++i) // allele
    {
	d=1; k=0;
	for (; k<d;k++) // trans
	{
	    for (j=0; j<=t->n; ++j) // rows
	    {
		mval_t *f= &t->trans[i][j];

		if ((f->type & IS_TRANS) && f->m > 1)
		{
		    if (d != 1 && d != f->m) error("ERROR: TRANS tags have different transcripts!\n");
		    d=f->m;
		}
		if (f->type & IS_TRANS)
		{
		    kputs(f->a[k].s, s);
		}
		else
		{
		    kputs(f->a[0].s, s);
		}
	    }
	    //   kputc('\n', s);
	}
    }
    return s->l -l_ori;
}

static convert_t *convert = NULL;

int convert_header(kstring_t *str)
{
    int l_ori = str->l;
    int i;
    kputc('#', str);
    for ( i=0; i<convert->nfmt; ++i )
    {
	if ( convert->fmt[i].is_gtf && !(split_flag & SPLIT_SAMPLE ))
	{
	    int j = i, js, k;
	    while ( j<convert->nfmt && convert->fmt[j].is_gtf ) j++;
	    for ( js = 0; js < bcf_hdr_nsamples(header); ++js )
	    {
		for ( k=i; k<j; ++k )
		{
		    if ( convert->fmt[k].type & IS_SEP )
		    {
			if ( convert->fmt[k].key ) kputs( convert->fmt[k].key, str);
		    }
		    else
			ksprintf(str, "%s:%s", header->samples[js], convert->fmt[k].key);
		}
	    }
	    i = j -1; // skip sample fields
	    continue;
	}
	if ( convert->fmt[i].type & IS_SAM ) {
	    if ( !(split_flag & SPLIT_SAMPLE) ) error("split tag has inconsistent format\n");
	    kputs("SAMPLE", str);
	    continue;
	}
	// Fixed fields
	if ( convert->fmt[i].key )
	{
	    if (!strcmp(convert->fmt[i].key, "BED"))
	    {
		kputs("CHROM\tSTART\tSTOP", str);
	    }
	    else
		kputs(convert->fmt[i].key, str);
	}
    }
    return str->l - l_ori;
}
void init_tags(tags_t *t, bcf1_t *const line)
{
    int i, j;
    int mrow;
    mrow = split_flag & SPLIT_ALT ? line->n_allele : 1;
    if ( split_flag & SPLIT_SAMPLE ) mrow *= header->n[BCF_DT_SAMPLE];
    if (t->k == 0) t->k = convert->nfmt;
    if ( t->l < mrow )
    {
	if (t->l == 0) 	t->trans = (mval_t**)malloc(mrow*sizeof(mval_t*));
	else 
	    t->trans = (mval_t**)realloc(t->trans, mrow*sizeof(mval_t*));
	for (i=t->l; i <mrow; ++i)
	{
	    t->trans[i] = (mval_t*)calloc(t->k, sizeof(mval_t));
	    for (j=0; j<t->k; ++j)
	    {
		t->trans[i][j].m = 0;
		t->trans[i][j].n = 1;
		t->trans[i][j].type = 0;
		t->trans[i][j].a = (kstring_t*)calloc(t->trans[i][j].n, sizeof(kstring_t));
	    }
	}
	t->l = mrow;
    }
    for (i=0; i<convert->nfmt; ++i)
    {
    	//if (convert->fmt[i].id == -1) continue;
	if ( convert->fmt[i].type == S_FORMAT)
	{
	    for (j=0; j<(int)line->n_fmt; j++)
	    {
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
    //tags_t tag = { 0, 0, 0};
    init_tags(tag, line);
    int isample=-1, nsample=0;
    int row = 0;

    if ( split_flag & SPLIT_SAMPLE )
    {
	isample = 0;
	nsample = header->n[BCF_DT_SAMPLE];
    }
    for ( ; isample<nsample; ++isample )
    {
	int iala=NONFLAG;
	int lala=iala;
	if ( split_flag & SPLIT_ALT )
	{
	    bcf_fmt_t *fgt = bcf_get_fmt(header, line, "GT");
#define BRANCH(type_t, vector_end) do {					\
		type_t *ptr = (type_t*) (fgt->p + isample*fgt->size);	\
		iala = (ptr[k]>>1)-1;					\
	    } while(0)
	    for ( k=0; k<fgt->n; ++k)
	    {
		switch (fgt->type) {
		    case BCF_BT_INT8:  BRANCH(int8_t, bcf_int8_vector_end); break;
		    case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_vector_end); break;
		    case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_vector_end); break;
		    default: fprintf(stderr, "FIXME: type %d in bcf_format_gt?\n", fgt->type); abort(); break;
		}
		if (iala == -1) continue;
		if (iala != NONFLAG)
		{
		    if (lala == iala) continue;
		    lala = iala;
		}
		if ( skip_ref && iala==0) continue;
		tag->m = row++;
		for ( i=0; i<convert->nfmt; i++ )
		{
		    tag->n = i;
		    if ( convert->fmt[i].handler ) convert->fmt[i].handler(line, &convert->fmt[i], iala, isample, tag);

		}
	    }
#undef BRANCH
	}
	else
	{
	    tag->m = row++;	    
	    for ( i=0; i<convert->nfmt; i++ )
	    {
		tag->n = i;
		if ( convert->fmt[i].handler ) convert->fmt[i].handler(line, &convert->fmt[i], iala, isample, tag);
	    }
	}
    }
    if (tag->m != -1 && tag->n != -1) tags2str(tag, str);
    clear_tag(tag);
    return str->l - l_ori;
}

void free_list(char **s, int n)
{
    int i;
    for (i = 0; i < n; ++i) free(s[i]);
    free(s);
}
int init_split_flag(char *s)
{
    int n, i;
    int flag = 0;
    char **list = hts_readlist(s, 0, &n);

    for (i = 0; i < n; ++i)
    {
	if (strcmp(list[i], "NONE") == 0) { free_list(list, n); return SPLIT_NONE; }
	else if (strcmp(list[i], "ALT") == 0) flag |= SPLIT_ALT;
	else if (strcmp(list[i], "SAMPLE") == 0) flag |= SPLIT_SAMPLE;
	else if (strcmp(list[i], "TRANS") == 0) flag |= (SPLIT_TRANS|SPLIT_ALT);
	else if (strcmp(list[i], "ALL") == 0) flag |= SPLIT_ALL;
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
    //kstring_t str = { 0, 0, 0 };
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    if ( line->n_allele == 1 ) kputc('.', &pv->a[pv->m++]);
    else kputs(line->d.allele[1],  &pv->a[pv->m++]);
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_chrom(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv= &tag->trans[tag->m][tag->n];
    kputs(header->id[BCF_DT_CTG][line->rid].key, &pv->a[pv->m++]);
    pv->type = fmt->type;
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }

}
static void process_pos(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv = &tag->trans[tag->m][tag->n];
    kputw(line->pos+1, &pv->a[pv->m++]);
    pv->type = fmt->type;
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_bed(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv = &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    uint32_t end = line->pos + strlen(line->d.allele[0]);
    kputs(header->id[BCF_DT_CTG][line->rid].key, &pv->a[pv->m]);
    kputc('\t', &pv->a[pv->m]);
    kputw(line->pos, &pv->a[pv->m]);
    kputc('\t', &pv->a[pv->m]);
    kputw(end, &pv->a[pv->m]);
    pv->m++;
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_ref(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv= &tag->trans[tag->m][tag->n];
    kputs(line->d.allele[0], &pv->a[pv->m++]);
    pv->type = fmt->type;
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_id(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv= &tag->trans[tag->m][tag->n];
    kputs(line->d.id, &pv->a[pv->m++]);
    pv->type = fmt->type;
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_alt(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    if ( line->n_allele == 1 ) kputc('.', &pv->a[pv->m]);
    else if ( iala == NONFLAG )
    {
	int i;
	for (i=1; i<line->n_allele; ++i)
	{
	    if ( i>1 ) kputc(',', &pv->a[pv->m]);
	    kputs(line->d.allele[i], &pv->a[pv->m]);
	}
    }
    else kputs(line->d.allele[iala], &pv->a[pv->m]);
    pv->m++;
    if (pv->n == pv->m)
    {
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
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_filter(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv= &tag->trans[tag->m][tag->n];
    int i;
    if ( line->d.n_flt )
    {
	for (i=0; i<line->d.n_flt; i++)
	{
	    if (i) kputc(';', &pv->a[pv->m]);
	    kputs(header->id[BCF_DT_ID][line->d.flt[i]].key, &pv->a[pv->m]);
	}
    }
    else kputc('.', &pv->a[pv->m]);
    pv->m++;
    pv->type = fmt->type;
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_tgt(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    if (iala != NONFLAG ) error("%TGT only used with SPLIT_GT \n");
    mval_t *pv= &tag->trans[tag->m][tag->n];
    if ( fmt->fmt==NULL ) kputc('.', &pv->a[pv->m]);
    else
    {
	assert(fmt->fmt->type==BCF_BT_INT8);
	int l;
	int8_t *x = (int8_t*)(fmt->fmt->p + isample*fmt->fmt->size);
	for (l=0; l<fmt->fmt->n && x[l]!=bcf_int8_vector_end; ++l)
	{
	    if (l) kputc("/|"[x[l]&1], &pv->a[pv->m]);
	    if (x[l]>>1) kputs(line->d.allele[(x[l]>>1)-1], &pv->a[pv->m]);
	    else kputc('.', &pv->a[pv->m]);
	}
	if (l==0) kputc('.', &pv->a[pv->m]);
    }
    pv->m++;
    pv->type = fmt->type;
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_trans(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    mval_t *pv= &tag->trans[tag->m][tag->n];
    if ( fmt->fmt==NULL ) kputc('.', &pv->a[pv->m]);
    else if (split_flag & SPLIT_ALT)
    {
	
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
    if ( iala==NONFLAG || type==BCF_VL_FIXED || type==BCF_VL_G )
    {
	
	if ( info->len == 1 )
	{
	    switch (info->type)
	    {
	    case BCF_BT_INT8: if ( info->v1.i==bcf_int8_missing ) kputc('.', &pv->a[pv->m]); else kputw(info->v1.i, &pv->a[pv->m]); break;
	    case BCF_BT_INT16: if ( info->v1.i==bcf_int16_missing ) kputc('.', &pv->a[pv->m]); else kputw(info->v1.i, &pv->a[pv->m]); break;
	    case BCF_BT_INT32: if ( info->v1.i==bcf_int32_missing ) kputc('.', &pv->a[pv->m]); else kputw(info->v1.i, &pv->a[pv->m]); break;
	    case BCF_BT_FLOAT: if ( bcf_float_is_missing(info->v1.f) ) kputc('.', &pv->a[pv->m]); else ksprintf(&pv->a[pv->m], "%g", info->v1.f); break;
	    case BCF_BT_CHAR: kputc(info->v1.i, &pv->a[pv->m]); break;
	    default: fprintf(stderr,"todo: type %d\n", info->type); exit(1); break;
	    }
	
	}
	else
	    bcf_fmt_array(&pv->a[pv->m], info->len, info->type, info->vptr);
	goto enlar;
    }
    assert( iala>=0 );
    if ( type==BCF_VL_A || type==BCF_VL_R )
    {
	
	if ( iala==0 && type==BCF_VL_A ) { kputc('.', &pv->a[pv->m]); goto enlar; }
	int j;
	int type = bcf_hdr_id2type(header, BCF_HL_INFO, fmt->id);
	if ( type==BCF_HL_STR )
	{
	    assert(info->type==BCF_BT_CHAR);
	    char *data = (char*)info->vptr;
	    for (j=0; j<info->len && data[j]!=bcf_str_missing; ++j);
	    kputsn(data, j, &pv->a[pv->m]);
	    goto enlar;
	}

	j = type==BCF_VL_A ? iala-1 : iala ;
#define BRANCH(type_t, is_missing, is_vector_end, kprint) { \
	    type_t *p = (type_t*)info->vptr; \
	    if ( is_vector_end || is_missing ) kputc('.', &pv->a[pv->m]); \
	    else kprint; \
	}
	switch (info->type)
	{
	case BCF_BT_INT8: BRANCH(int8_t, p[j]==bcf_int8_missing, p[j]==bcf_int8_vector_end, kputw(p[j], &pv->a[pv->m])); break;
	case BCF_BT_INT16: BRANCH(int16_t, p[j]==bcf_int16_missing, p[j]==bcf_int16_vector_end, kputw(p[j], &pv->a[pv->m])); break;
	case BCF_BT_INT32: BRANCH(int32_t, p[j]==bcf_int32_missing, p[j]==bcf_int32_vector_end, kputw(p[j], &pv->a[pv->m])); break;
	case BCF_BT_FLOAT: BRANCH(float, bcf_float_is_missing(p[j]), bcf_float_is_vector_end(p[j]), ksprintf(&pv->a[pv->m], "%g", p[j])); break;
	case BCF_BT_CHAR:  BRANCH(char, p[j]==bcf_str_missing, NULL, kputc(p[j], &pv->a[pv->m])); break;
	default: fprintf(stderr, "todo: type %d\n", type); exit(1); break;
	}
#undef BRANCH
    }

enlar:
    pv->m++;
    if (pv->n == pv->m)
    {
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

    if ( iala==NONFLAG || fmt->fmt->type==BCF_VL_FIXED || fmt->fmt->type==BCF_VL_G )
    {
	bcf_fmt_array(&pv->a[pv->m], fmt->fmt->n, fmt->fmt->type, fmt->fmt->p+isample*fmt->fmt->size);
	goto fmtenlarge;
    }
    assert(iala >= 0);
    if ( fmt->fmt->type==BCF_VL_A || fmt->fmt->type==BCF_VL_R )
    {
	int j;
	j = fmt->fmt->type==BCF_VL_A ? iala-1 : iala;
#define BRANCH(type_t, is_missing, is_vector_end, data, kprint) {	\
	    type_t *p = (type_t*)data; \
	    if ( is_vector_end || is_missing ) kputc('.', &pv->a[pv->m]); \
	    else kprint; \
	}
	switch (fmt->fmt->type)
	{
	case BCF_BT_INT8: BRANCH(int8_t, p[j]==bcf_int8_missing, p[j]==bcf_int8_vector_end, fmt->fmt->p+isample*fmt->fmt->size, kputw(p[j], &pv->a[pv->m])); break;
	case BCF_BT_INT16: BRANCH(int16_t, p[j]==bcf_int16_missing, p[j]==bcf_int16_vector_end, fmt->fmt->p+isample*fmt->fmt->size, kputw(p[j], &pv->a[pv->m])); break;
	case BCF_BT_INT32: BRANCH(int32_t, p[j]==bcf_int32_missing, p[j]==bcf_int32_vector_end, fmt->fmt->p+isample*fmt->fmt->size, kputw(p[j], &pv->a[pv->m])); break;
	case BCF_BT_FLOAT: BRANCH(float, bcf_float_is_missing(p[j]), bcf_float_is_vector_end(p[j]), fmt->fmt->p+isample*fmt->fmt->size, ksprintf(&pv->a[pv->m], "%g", p[j])); break;
	case BCF_BT_CHAR:  BRANCH(char, p[j]==bcf_str_missing, NULL, fmt->fmt->p+isample*fmt->fmt->size, kputc(p[j], &pv->a[pv->m])); break;
	default: fprintf(stderr, "todo: type %d\n", fmt->fmt->type); exit(1); break;
	}
#undef BRANCH
	goto fmtenlarge;
    }

fmtenlarge:
    pv->m++;
    if (pv->n == pv->m)
    {
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
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}
static void process_sep(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    if (!fmt->key) error("FIXME: This is an empty fmt\n");
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    kputs(fmt->key,&pv->a[pv->m++]);
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (kstring_t*)realloc(pv->a, pv->n*sizeof(kstring_t));
    }
}

/* The register_tag and parse_tag functions are adapted from pd3's convert.c 
 * use S-tags instead of T-tags 
 */
fmt_t *register_tag(int type, char *key, int is_gtf)
{
    convert->nfmt++;
    if (convert->nfmt == convert->mfmt)
    {
	convert->mfmt += 10;
	convert->fmt = (fmt_t*)realloc(convert->fmt, convert->mfmt*sizeof(fmt_t));
    }
    fmt_t *fmt = &convert->fmt[ convert->nfmt-1 ];
    fmt->type = type;
    fmt->is_gtf = is_gtf;
    fmt->key = key ? strdup(key) : NULL;
    // allow non-format tags to appear amongst the format fields, split-samples mode
    if ( key )
    {
	fmt->id = bcf_hdr_id2int(header, BCF_DT_ID, fmt->key);
	if (fmt->type==S_FORMAT && !bcf_hdr_idinfo_exists(header, BCF_HL_FMT, fmt->id))
	{
	    if ( !strcmp("CHROM", key) ) fmt->type = S_CHROM;
	    else if ( !strcmp("POS", key) ) fmt->type = S_POS;
	    else if ( !strcmp("BED", key) ) fmt->type = S_BED;
	    else if ( !strcmp("ID", key) ) fmt->type = S_ID;
	    else if ( !strcmp("REF", key) ) fmt->type = S_REF;
	    else if ( !strcmp("ALT", key) ) fmt->type = S_ALT;
	    else if ( !strcmp("FIRST_ALT", key) ) fmt->type = S_FIRST_ALT;
	    else if ( !strcmp("QUAL", key) ) fmt->type = S_QUAL;
	    else if ( !strcmp("FILTER", key) ) fmt->type = S_FILTER;
	    else if ( !strcmp("SAMPLE", key) ) fmt->type = S_SAMPLE;
	    else if ( !strcmp("HGVS", key) ) fmt->type = S_TRANS;
	    else if ( fmt->id>=0 && bcf_hdr_idinfo_exists(header, BCF_HL_INFO, fmt->id))
	    {
		fmt->type = S_INFO;
	    }
	}
    }
    /* if (fmt->type == S_SEP) */
    /* { */
    /* 	fmt->id = -1; */
    /* } */
    /* else */
    /* { */
    /* 	fmt->id = bcf_hdr_id2int(header, BCF_DT_ID, fmt->key); */
    /* 	if ( fmt->id == -1 && !fmt->is_gtf) */
    /* 	{ */
    /* 	    if ( !strcmp("CHROM", key) ) fmt->type = S_CHROM; */
    /* 	    else if ( !strcmp("POS", key) ) fmt->type = S_POS; */
    /* 	    else if ( !strcmp("BED", key) ) fmt->type = S_BED; */
    /* 	    else if ( !strcmp("ID", key) ) fmt->type = S_ID; */
    /* 	    else if ( !strcmp("REF", key) ) fmt->type = S_REF; */
    /* 	    else if ( !strcmp("ALT", key) ) fmt->type = S_ALT; */
    /* 	    else if ( !strcmp("FIRST_ALT", key) ) fmt->type = S_FIRST_ALT; */
    /* 	    else if ( !strcmp("QUAL", key) ) fmt->type = S_QUAL; */
    /* 	    else if ( !strcmp("FILTER", key) ) fmt->type = S_FILTER; */
    /* 	    else if ( !strcmp("SAMPLE", key) ) fmt->type = S_SAMPLE; */
    /* 	    else  error("No such tag in the header %d\n", key); */
    /* 	} */
    /* 	/\* else *\/ */
    /* 	/\* { *\/ */
    /* 	/\*     if ( !strcmp("HGVS", key) ) fmt->type = S_TRANS; *\/ */
    /* 	/\*     else if ( !strcmp("SAMPLE", key) ) fmt->type = S_SAMPLE; *\/ */
    /* 	/\*     else if ( fmt->id >= 0 && bcf_hdr_idinfo_exists(header, BCF_HL_INFO, fmt->id) ) *\/ */
    /* 	/\*     { *\/ */
    /* 	/\* 	fmt->type = S_INFO; *\/ */
    /* 	/\* 	fprintf(stderr, "Warning: Assuming INFO %s\n", key); *\/ */
    /* 	/\*     } *\/ */
    /* 	/\*     else *\/ */
    /* 	/\*     { *\/ */
    /* 	/\* 	fmt->type = S_FORMAT; *\/ */
    /* 	/\*     } *\/ */
    /* 	/\*}*\/ */
    /* } */
    switch ( fmt->type )
    {
    case S_FIRST_ALT: fmt->handler = &process_first_alt; break;
    case S_CHROM: fmt->handler = &process_chrom; break;
    case S_POS: fmt->handler = &process_pos; break;
    case S_BED: fmt->handler = &process_bed; break;
    case S_ID: fmt->handler = &process_id; break;
    case S_REF: fmt->handler = &process_ref; break;
    case S_ALT: fmt->handler = &process_alt; break;
    case S_QUAL: fmt->handler = &process_qual; break;
    case S_FILTER: fmt->handler = &process_filter; convert->max_unpack |= BCF_UN_FLT; break;
    case S_INFO: fmt->handler = &process_info; convert->max_unpack |= BCF_UN_INFO; break;	
    case S_TRANS: fmt->handler = &process_trans; convert->max_unpack |= BCF_UN_INFO; break;
    case S_FORMAT: fmt->handler = &process_format; convert->max_unpack |= BCF_UN_FMT; break;
    case S_SAMPLE: fmt->handler = &process_sample; break;
    case S_SEP: fmt->handler = &process_sep; break;
    case S_TGT: fmt->handler = &process_tgt; convert->max_unpack |= BCF_UN_FMT; break;
    default: error("TODO: handler for type %d\n", fmt->type);
    }
    if ( key )
    {
	if ( fmt->type==S_INFO )
	{
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
    if ( is_gtf )
    {
	if ( !strcmp(str.s, "SAMPLE") ) register_tag(S_SAMPLE, "SAMPLE", is_gtf);
	else if ( !strcmp(str.s, "GT") || !strcmp(str.s, "TGT")) register_tag(S_TGT, "GT", is_gtf);
	else if ( !strcmp(str.s, "IUPACGT") ) register_tag(S_IUPAC_GT, "GT", is_gtf);
	else
	{
	    register_tag(S_FORMAT, str.s, is_gtf);
	}
    }
    else
    {
	if ( !strcmp(str.s, "CHROM") ) register_tag(S_CHROM, "CHROM", is_gtf);
	else if ( !strcmp(str.s, "POS") ) register_tag(S_POS, "POS", is_gtf);
	else if ( !strcmp(str.s, "BED") ) register_tag(S_BED, "BED", is_gtf);
	else if ( !strcmp(str.s, "POS") ) register_tag(S_CHROM, "POS", is_gtf);
	else if ( !strcmp(str.s, "REF") ) register_tag(S_REF, "REF", is_gtf);
	else if ( !strcmp(str.s, "ALT") ) register_tag(S_ALT, "ALT", is_gtf);
	else if ( !strcmp(str.s, "HGVS") ) register_tag(S_TRANS, "HGVS", is_gtf);
	else if ( !strcmp(str.s, "GT") || !strcmp(str.s, "TGT")) { split_flag |= SPLIT_SAMPLE; register_tag(S_TGT, "GT", is_gtf); }
	else if ( !strcmp(str.s, "IUPACGT") ) { split_flag |= SPLIT_SAMPLE; register_tag(S_IUPAC_GT, "GT", is_gtf); }
	else if ( !strcmp(str.s, "FUNC") ) register_tag(S_TRANS, "FUNC", is_gtf);
	else if ( !strcmp(str.s, "SAMPLE") ) { split_flag |= SPLIT_SAMPLE; register_tag(S_SAMPLE, "SAMPLE", is_gtf); }
	else if ( !strcmp(str.s, "FIRST_ALT") ) register_tag(S_FIRST_ALT, "FIRST_ALT", is_gtf);
	else if ( !strcmp(str.s, "QUAL") ) register_tag(S_QUAL, "QUAL", is_gtf);
	else if ( !strcmp(str.s, "INFO") )
	{
	    if ( *q=='/' ) error("Could not parse format string: %s\n", convert->format_str);
	    p = ++q;
	    str.l = 0;
	    while ( *q && (isalnum(*q) || *q=='_' || *q=='.') ) q++;
	    if ( q-p==0 ) error("Could not parse format string: %s\n", convert->format_str);
	    kputsn(p, q-p, &str);
	    register_tag(S_INFO, str.s, is_gtf);
	}
	else
	{
	    register_tag(S_INFO, str.s, is_gtf);
	}
    }
    free(str.s);
    return q;
}
static char *parse_sep(char *p, int is_gtf)
{
    char *q = p;
    kstring_t str = { 0, 0, 0};
    while ( *q && *q!='[' && *q!=']' && *q!='%' )
    {
	if ( *q=='\\' )
	{
	    q++;
	    if ( *q=='n' ) kputc('\n', &str);
	    else if ( *q == 't') kputc('\t', &str);
	    else kputc(*q, &str);
	}
	else kputc(*q, &str);
	q++;
    }
    if ( !str.l ) error("Could not parse format string: %s\n", convert->format_str);
    register_tag(S_SEP, str.s, is_gtf);
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
    while ( *p )
    {
	switch (*p)
	{
	case '[': is_gtf = 1; p++; break;
	case ']': is_gtf = 0; p++; break;
	case '%': p = parse_tag(p, is_gtf); break;
	default:  p = parse_sep(p, is_gtf); break;
	}
    }
}
const char *about(void)
{
    return "Select tags in pre-defined format.\n";
}

const char *usage(void)
{
    return
	"\n"
	"About : Select tags from VCF/BCF file.\n"
	"Usage:\n"
	"Standalone mode:\n"
	"\tbcfselect [Options] in.vcf.gz\n"
	"Options:\n"
	"\t-f, --format   see man page for deatils.\n"
	"\t-s, --split    split by [ALT,SAMPLE,TRANS].\n"
	"\t-p, --print-header  print the header comment.\n"
	"BCFtools plugin mode:\n"
	"\tbcftools +select [General Options] -- [Plugin Options]\n"
	"General Options:\n"
	"\trun \"bcftools plugin\" for a list of common options.\n"
	"Plugin Options: same with standalone mode.\n"
	"\n"
	"Website:\n";
}

kstring_t *mempool = NULL;


int run(int argc, char**argv)
{
    struct option const long_opts[] =
	{
	    {"format", required_argument, NULL, 'f'},
	    {"skip-ref", no_argument, NULL, 'r'},
	    {"split", required_argument, NULL, 's'},
	    {"print-header", no_argument, NULL, 'p'},
	    {0, 0, 0, 0}
	};
    char c;
    char *format = NULL, *flag = NULL;
    while ((c = getopt_long(argc, argv, "f:s:rph?", long_opts, NULL)) >= 0)
    {
	switch (c)
	{
	case 'f': format = strdup(optarg); break;
	case 'r': skip_ref = 1; break;
	case 's': flag = strdup(optarg); break;
	case 'p': print_header = 1; break;
	case 'h':
	case '?':
	default: error("%s", usage()); break;
	}
    }
    bcf_srs_t *sr = bcf_sr_init();
    sr->require_index = 0;
    /* if ( optind==argc || !strcmp(argv[optind],"-")) */
    /* { */
    /* 	if ( !isatty(fileno((FILE*)stdin)) ) fname = "-"; */
    /* 	else error("%s", usage()); */
    /* } */
    /* else */
    /* { */
    /* 	fname = argv[optind]; */
    /* } */
    /* if ( !bcf_sr_add_reader(sr, fname) ) error("Failed to open %s : %s\n", fname, bcf_sr_strerror(sr->errnum)); */
    if (!argv[0]) error("%s", usage());
    if (!bcf_sr_add_reader(sr, argv[0])) error("Failed to open %s : %s\n", argv[optind], bcf_sr_strerror(sr->errnum));
    header = sr->readers[0].header;

    if (flag)
    {
	split_flag = init_split_flag(flag);
	free(flag);
    }
    if (format)
    {
	convert_init(format);
	free(format);
    }
    else
    {
	error("Must set a format string with -f \n");
    }
    
    mempool = (kstring_t*)malloc(sizeof(kstring_t));
    mempool->m = mempool->l = 0;
    mempool->s = NULL;

    tag = (tags_t*)malloc(sizeof(tags_t));
    tag->l = tag->k = 0;

    if (print_header) convert_header(mempool);

    while (bcf_sr_next_line(sr))
    {
	bcf1_t *line = bcf_sr_get_line(sr, 0);
	convert_line(line, mempool);
	if (mempool->l > MEMPOOL)
	{
	    fprintf(stdout, "%s", mempool->s);
	    mempool->l = 0;
	}
    }
    if (mempool->l) { fprintf(stdout, "%s", mempool->s); mempool->l = 0; }
    bcf_sr_destroy(sr);
    destroy_tag(tag);
    if (mempool->m) free(mempool->s);
    free(mempool);
    return 0;
}


#ifdef _SELECT_MAIN
int main(int argc, char **argv)
{
    return run(argc-1, argv+1);
}
#endif

#ifndef BCFTOOLS_VERSION
#define BCFTOOLS_VERSION "1.2"

static void version(const char **bcftools_version, const char **htslib_version)
{
    *bcftools_version = BCFTOOLS_VERSION;
    *htslib_version = hts_version();
}

static void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

#endif
