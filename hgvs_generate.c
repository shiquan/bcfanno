#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "utils.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/khash.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "htslib/faidx.h"
#include "hgvs.h"

char str_init[2];

KHASH_MAP_INIT_STR(name, char*)

typedef khash_t(name) namehash_type;

typedef char *(*copy_seqs_func)(const char*, unsigned long n);

struct name_list {
    int l, m;
    char **a;
};

static void genepred_line_clear(struct genepred_line *line)
{
    if ( line->clear == 1 )
	return;
    line->clear = 1;
    if (line->clear == 0) {
	free(line->chrom);
	free(line->name1);
	if ( line->name2 )
	    free(line->name2);
	free(line->exons);
	free(line->dna_ref_offsets);
    }
    memset(line, 0, sizeof(struct genepred_line));
}
static void hgvs_core_clear(struct hgvs_core *c)
{
    c->str.l = 0;
    c->l_name = 0;
    c->l_type = 0;
    c->type = func_region_unknown;
}
static void hgvs_clear(struct hgvs *name)
{
    int i;
    for ( i = 0; i < name->l; ++i )
	hgvs_core_clear(&name->a[i]);				   
    name->l = 0;
}
void hgvs_cache_clear(struct hgvs_cache *cache)
{
    int i;
    for ( i = 0; i < cache->l; ++i )
	hgvs_clear(&cache->a[i]);
    cache->l = 0;
}
void hgvs_cache_destroy(struct hgvs_cache *cache)
{
    int i, j;
    for ( i = 0; i < cache->i; ++i ) {
	struct hgvs *name = &cache->a[i];
	for ( j = 0; j < name->i; ++j ) {
	    kstring_t *str = &name->a[j].str;
	    if ( str->m ) 
		free(str->s);	    
	}
	free(name->a);	    
    }
    free(cache->a);
}
namehash_type *init_gene_name(const char *name)
{
    int i, n = 0;
    khiter_t k;
    int ret;
    namehash_type *hash;
    char **names = hts_readlist(name, 1, &n);

    if ( n == 0 )
	return NULL;
    hash = kh_init(name);

    for ( i = 0; i < n; ++n ) {
       k = kh_put(name, hash, names[i], &ret);
    }
    return hash;
}

tbx_t *load_genepred_file(const char *fn)
{
    return tbx_index_load(fn);
}
faidx_t *load_refseq_file(const char *fn)
{
    return fai_load(fn);
}
uint32_t bcf_calend(bcf1_t *line)
{
    return line->pos + line->rlen;
}
int hgvs_bcf_header_add_gene(bcf_hdr_t *hdr)
{
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, "Gene");
    if (id == -1) {
	bcf_hdr_append(hdr, "##INFO=<ID=Gene,Number=1,Type=String,Description=\"Gene names\">");
	bcf_hdr_sync(hdr);
	id = bcf_hdr_id2int(hdr, BCF_DT_ID, "Gene");
	assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));
    }
    return id;
}
int hgvs_bcf_header_add_trans(bcf_hdr_t *hdr)
{
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, "Transcript");
    if (id == -1) {
	bcf_hdr_append(hdr, "##INFO=<ID=Transcript,Number=G,Type=String,Description=\"Transcript names\">");
	bcf_hdr_sync(hdr);
	id = bcf_hdr_id2int(hdr, BCF_DT_ID, "Transcript");
	assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));
    }
    return id;
}
int hgvs_bcf_header_add_dna(bcf_hdr_t *hdr)
{
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, "HGVSDNA");
    if (id == -1) {
	bcf_hdr_append(hdr, "##INFO=<ID=HGVSDNA,Number=G,Type=String,Description=\"HGVS nomenclature for the description of DNA sequence variants\">");
	bcf_hdr_sync(hdr);
	id = bcf_hdr_id2int(hdr, BCF_DT_ID, "HGVSDNA");
	assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));
    }
    return id;
}

static struct genepred_format refgene_formats = {
    .chrom = 2,
    .name1 = 1,
    .name2 = 12,
    .strand = 3,
    .txstart = 4,
    .txend = 5,
    .cdsstart = 6,
    .cdsend = 7,
    .exoncount = 8,
    .exonstarts = 9,
    .exonends = 10,
};
    
static struct genepred_format genepred_formats = {    
    .name1 = 0,
    .chrom = 1,
    .strand = 2,
    .txstart = 3,
    .txend = 4,
    .cdsstart = 5,
    .cdsend = 6,
    .exoncount = 7,
    .exonstarts = 8,
    .exonends = 9,
    .name2 = 10,
};

void genepred_prase_core(kstring_t *string, struct genepred_line *line, struct genepred_format *type)
{
    int nfields = 0;
    genepred_line_clear(line);
    
    // accept any gene_pred-like format, like refGene, ensGene, gene_pred, see our manual how to get these databases
    // split the string by tab
    int *splits = ksplit(string, 0, &nfields);
    // chromosome name
    line->chrom = strdup(string->s + splits[type->chrom]);
    // usually gene names or ensemble gene id
    line->name1 = strdup(string->s + splits[type->name1]);
    // strand, char '+' or '-'    
    line->strand = memcmp(string->s + splits[type->strand], "+", 1) ? '-' : '+';
    // trans start
    line->txstart = atoi(string->s + splits[type->txstart]);
    // trans end
    line->txend = atoi(string->s + splits[type->txend]);
    // cds start, for mRNA cds start should greater than txstart, for ncRNA cdsstart should equal to txend
    line->cdsstart = atoi(string->s + splits[type->cdsstart]);
    // cds end, for ncRNA cdsend == cdsstart == txend
    line->cdsend = atoi(string->s + splits[type->cdsend]);;
    // exon number
    line->exoncount = atoi(string->s + splits[type->exoncount]);   
    // for some genePred file, no name2 specified.
    line->name2 = strdup(string->s + splits[type->name2]);
    
    // the exons region in gene_pred file is like  1,2,3,4,5. try to init the exon_pair[] 
    // and exon_offset_pair[] by exoncount 
    char *ss = string->s + splits[type->exonstarts];
    char *se = string->s + splits[type->exonends];
    
    line->exons = (struct exon_pair*)calloc(line->exoncount, sizeof(struct exon_pair));
    line->dna_ref_offsets = (struct exon_offset_pair*)calloc(line->exoncount, sizeof(struct exon_offset_pair));
    free(splits);
    int i;
    char *ss1, *se1;
    for ( i = 0; i < line->exoncount; ++i ) {
	// start
	ss1 = ss;
	while(*ss1 && *ss1 != ',') ss1++;
	ss1[0] = '\0';
	line->exons[i].start = atoi(ss);
	ss = ++ss1; // skip ','
	// end
	se1 = se;
	while(*se && *se1 != ',') se1++;
	se1[0] = '\0';
	line->exons[i].end = atoi(se);	
	se = ++se1; // skip ','
    }
    line->clear = 0;
}
void genepred_line_praser(kstring_t *string, struct genepred_line *line, struct genepred_format *type)
{
    // prase string into genepred line struct
    genepred_prase_core(string, line, type);
    // calculate the length of function regions, for plus strand, forward length is the length of UTR5, and
    // backward length is the length of UTR3, for minus strand 
    int ref_len = 0;
    int read_len = 0;
    int forward_len = 0;
    int backward_len = 0;
    int is_coding = line->cdsstart == line->cdsend ? 0 : 1;
    int i;
    for ( i = 0; i < line->exoncount; ++i ) {
	// in genepred format, start is 0 based, end is 1 based
	int exon_len = line->exons[i].end - line->exons[i].start;
	// add exon length to dna reference length
	ref_len += exon_len;	
	// only mRNA has UTRs
	if (is_coding == 0) continue;

	if (line->exons[i].end < line->cdsstart) {
	    forward_len += exon_len;
	    continue;
	}
	// first cds 
	if (line->cdsstart > line->exons[i].start) {
	    forward_len += line->cdsstart - line->exons[i].start;
	    continue;
	} 
	// come to end regions
	if (line->cdsend <= line->exons[i].start) {
	    backward_len += exon_len;
	    continue;
	}
	// last cds
	if (line->cdsend < line->exons[i].end) {
	    backward_len += line->exons[i].end - line->cdsend;
	}		   
    }

    // read_len is the coding reference length for coding transcript or length of noncoding transcript 
    read_len = ref_len - forward_len - backward_len;

    // Coding DNA reference:
    //                   -3                                            *3
    //                    -2                                          *2
    // -93   -45 -44       -1 1                187 188           351 *1   *96 *97    *223
    // |_______|____________|====================|=================|________|_______|
    //        / \            ATG                / \             TGA        / \
    //       /   \                             /   \                      /   \
    //      /gtga,,,g                         /gta,,,act                 /gtag,,,cc
    //  -45+1        -44-1                187+1        188-1         *96+1        *97-1
    //   -45+2      -44-2                  187+2      188-2           *96+2      *97-2
    //    -45+3    -44-3                    187+3    188-3             *96+3    *97-3
    // 
    // 
    // Noncoding DNA reference:
    // 
    // 1        49 50             280 281                 540 541             667
    // |__________|__________________|_______________________|_________________|
    // 
     
    int l1 = 0, l2 = line->exoncount - 1;
    // count forward
    int forward_offset = 0;
    int backward_offset = 0;
    for ( ; l1 < line->exoncount && forward_len; ++l1 ) {
	int32_t exon_len = line->exons[l1].end - line->exons[l1].start;

	line->dna_ref_offsets[l1].start = (forward_len<<OFFSET_BITWISE) | (line->strand == '+' ? REG_UTR5 : REG_UTR3);

	if ( forward_len > exon_len ) {
	    forward_len -= exon_len;	
	    line->dna_ref_offsets[l1].end = ((forward_len+1)<<OFFSET_BITWISE) |(line->strand == '+' ? REG_UTR5 : REG_UTR3);
	    continue;
	}
	if (line->strand == '+') {
	    forward_offset = exon_len - forward_len;
	} else {
	    forward_offset = read_len + forward_len - exon_len +1; // for minus strand count from backward
	}
	line->dna_ref_offsets[l1].end = (forward_offset<<OFFSET_BITWISE)| REG_CODING;
	forward_len = 0;
	++l1;
	break;
    }

    // count backward
    for ( ; l2 > 0 && backward_len; --l2 ) {
	int32_t exon_len = line->exons[l2].end - line->exons[l2].start;

	line->dna_ref_offsets[l2].end = (backward_len<<OFFSET_BITWISE) |(line->strand == '+' ? REG_UTR3 : REG_UTR5);

	if (backward_len > exon_len) {
	    backward_len -= exon_len;
	    line->dna_ref_offsets[l2].start = ((backward_len+1)<<OFFSET_BITWISE) |(line->strand=='+' ? REG_UTR3 : REG_UTR5);
	    continue;
	} 
	if (line->strand == '+') {
	    backward_offset = read_len + backward_len - exon_len + 1;
	} else {
	    backward_offset = exon_len - backward_len;
	}
	line->dna_ref_offsets[l2].start = (backward_offset<<OFFSET_BITWISE) | REG_CODING;
	backward_len = 0;
	--l2;
	break;
    }

    // count inter regions
    if (line->strand == '+') {
	if (is_coding)
	    read_len = backward_offset -1;       
	int l;
	for ( l = l2; l >= l1;  l-- ) {
	    int32_t exon_len = line->exons[l].end - line->exons[l].start;
	    line->dna_ref_offsets[l].end = (read_len<<OFFSET_BITWISE) |	(is_coding ? REG_CODING : REG_NONCODING);
	    read_len -= exon_len;
	    line->dna_ref_offsets[l].start = ((read_len+1)<< OFFSET_BITWISE) | (is_coding ? REG_CODING : REG_NONCODING);
	}
    } else {
	if (is_coding) {
	    read_len = forward_offset -1;
	}
	int l;
	for ( l = l1; l <= l2; l++ ) {
	    int32_t exon_len = line->exons[l].end - line->exons[l].start;
	    line->dna_ref_offsets[l].start = (read_len<<OFFSET_BITWISE) |(is_coding ? REG_CODING : REG_NONCODING);
	    read_len -= exon_len;
	    line->dna_ref_offsets[l].end = ((read_len+1)<<OFFSET_BITWISE) |(is_coding ? REG_CODING : REG_NONCODING);
	}
    }    
}

void generate_dbref_database(struct genepred_line *line)
{
    if (line->clear == 1) return;
    int i, j;
    for ( i = 0; i < line->exoncount; ++i ) {
    	kstring_t temp[2] = { KSTRING_INIT, KSTRING_INIT };
    	int types[2];
    	types[0] = line->dna_ref_offsets[i].start & REG_MASK;
    	types[1] = line->dna_ref_offsets[i].end & REG_MASK;
    	int j;
	// [start, end] 
    	for ( j = 0; j < 2; ++j ) {
    	    kstring_t *temp1 = &temp[j];
    	    int type = types[j];
    	    switch (type) {
    		case REG_UTR5 :
    		    kputc('-', temp1);
    		    break;
    		case REG_UTR3:
    		    kputc('*', temp1);
    		    break;
    		case REG_CODING:
    		    kputs("c.", temp1);
    		    break;
    		case REG_NONCODING:
    		    kputs("n.", temp1);
    		    break;
    		default:
    		    error("Unknown type : %d", type);
    	    }
    	}
    	kputw((line->dna_ref_offsets[i].start>>OFFSET_BITWISE), &temp[0]);
    	kputw((line->dna_ref_offsets[i].end>>OFFSET_BITWISE), &temp[1]);

	// format: CHROM,START,END,STRAND, GENE, TRANSCRIPT, EXON, START_LOC, END_LOC 
    	free(temp[0].s);
    	free(temp[1].s);
    }
}
static void push_mempool(struct genepred_memory_pool *pool, kstring_t *str, struct genepred_format *type)
{
    if (pool->m == pool->l) {
	pool->m = pool->m == 0 ? 2 : pool->m << 1;
	pool->a = (struct genepred_line *)realloc(pool->a, sizeof(struct genepred_line)*pool->m);
    }
    // i should always greater than l. if i == l increase i to init a new line for future use. go abort if i < l
    if (pool->i == pool->l) {
	pool->a[pool->i++].clear = -1;
    }
    // the prase func must return a point or go abort
    genepred_line_praser(str, &pool->a[pool->l], type);
    pool->l++;
}
static void update_mempool(struct genepred_memory_pool *pool)
{    
    if (pool->l == 0)
	return;
    // try to loop this pool and find the edges
    pool->start = pool->a[0].txstart;
    pool->end = pool->a[0].txend;
    int i;    
    for ( i = 1; i < pool->l; ++i ) {
	if ( pool->a[i].txstart < pool->start )
	    pool->start = pool->a[i].txstart;
	if ( pool->a[i].txend > pool->end )
	    pool->end = pool->a[i].txend;
    }
}
void hgvs_fill_memory(struct genepred_memory_pool *pool, htsFile *fp, tbx_t *tbx, int rid, uint32_t start, uint32_t end, struct genepred_format *type)
{
    hts_itr_t *itr = tbx_itr_queryi(tbx, rid, start, end);
    kstring_t str = KSTRING_INIT;
    if (rid == -1) return;
    
    while ( tbx_itr_next(fp, tbx, itr, &str) >= 0) {
	push_mempool(pool, &str, type);
	str.l = 0;
    }
    pool->rid = rid;
    if ( str.m )
	free(str.s);
    tbx_itr_destroy(itr);
    update_mempool(pool);    
}
static void hgvs_memory_clear(struct genepred_memory_pool *pool)
{
    int i;
    for ( i = 0; i < pool->i; ++i)
	genepred_line_clear(&pool->a[i]);
    
    pool->l = 0;
    // i and m should be kept for future use
    pool->start = 0;
    pool->end = 0;
    pool->rid = -1;
    hgvs_cache_destroy(&pool->cache);
}
void namehash_destroy(void *_hash)
{
    namehash_type *hash = (namehash_type*)_hash;
    khiter_t k;
    for ( k = 0; k < kh_end(hash); ++k ) {
	if ( kh_exist(hash, k) ) 
	    free((char*)kh_key(hash, k));
    }
    kh_destroy(name, hash);
}
struct refgene_options *refgene_opts_init()
{
    struct refgene_options *opts = (struct refgene_options*)malloc(sizeof(struct refgene_options));
    opts->genepred_fname = 0;
    opts->refseq_fname = 0;
    opts->fp = 0;
    opts->genepred_idx = 0;
    opts->check_refseq = 0;
    opts->refseq_fai = 0;
    opts->trans_list_fname = 0;
    opts->genes_list_fname = 0;
    opts->screen_by_genes = 0;
    opts->screen_by_transcripts = 0;
    opts->genehash = 0;
    opts->transhash = 0;
    return opts;    
}
void refgene_opts_destroy(struct refgene_options *opt)
{
    hts_close(opt->fp);
    tbx_destroy(opt->genepred_idx);
    if ( check_refseq == 1 )
	fai_destroy(opt->refseq_fai);
    if ( opt->screen_by_genes == 1 )
	namehash_destroy(opt->genehash);
    if ( opt->screen_by_transcripts == 1 )
	namehash_destroy(opt->transhash);
    hgvs_memory_clear(&opt->pool);
    free(opt);
}
void hgvs_des_destory(struct hgvs_des *des)
{
    if ( des->ref_length > 0 )
	free(des->ref);
    if ( des->alt_length > 0 )
	free(des->alt);
    free(des);
}
static char *rev_seqs(const char *dna_seqs, unsigned long n)
{
    static unsigned char rev_seqs_matrix[256] = {
	'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
	'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
	'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
	'T','G','C','A','N','N','N','N','N','N','N','N','N','N','N','N',
	'N','T','N','G','N','N','N','C','N','N','N','N','N','N','N','N',
	'N','N','N','N','A','N','N','N','N','N','N','N','N','N','N','N',
	'N','T','N','G','N','N','N','C','N','N','N','N','N','N','N','N',
	'N','N','N','N','A','N','N','N','N','N','N','N','N','N','N','N',
	'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
	'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
	'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
	'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
	'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
	'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
	'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
	'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N'
    };
    
    if (n == 0) return NULL;
    char *rev = (char*)calloc(n+1, sizeof(char));
    int i;
    for ( i = 0; i < n; ++i )
	rev[i] = rev_seqs_matrix[dna_seqs[n-i-1]];    
    rev[n] = '\0';
    return rev;
}
void hgvs_des_reverse(struct hgvs_des *des, uint8_t strand)
{
    if ( des->type == var_type_ref )
	return;
    if ( des->strand == strand )
	return;
    if ( des->strand == 0 && strand == '+' )
	return;
    des->strand = strand;
    char *ref = des->ref;
    char *alt = des->alt;
    int32_t start = des->end;
    des->end = des->start;
    des->start = start;
    des->ref = rev_seqs(des->ref, des->ref_length);
    des->alt = rev_seqs(des->alt, des->alt_length);
    if ( des->ref_length > 0 )
	free(ref);
    if ( des->alt_length > 0 )
	free(alt);
}
static struct hgvs_des *describe_variants(const char *ref, const char *alt, int _pos)
{
    struct hgvs_des *des = (struct hgvs_des*)malloc(sizeof(struct hgvs_des));
    const char *a = alt;
    const char *r = ref;
    // pos may differ from _pos, but _pos will not change in this function
    int pos = _pos;
    // skip same string from start
    while (*a && *r && toupper(*a) == toupper(*r)) { a++; r++; pos++; }
    
    if ( !a[0] && !r[0] ) {
	des->type = var_type_ref;
	return des;
    }
    
    des->strand = 0; // 0 for unknown, '+' for plus, '-' for minus
    // if ref and alternative allele are 1 base, take it as snp
    if ( *a && *r && !a[1] && !r[1] ) {
	// mpileip may output X allele, treat as ref
	if ( *a == '.' || *a == 'X' || *a == '*') {
	    des->type = var_type_ref;
	    return des;
	}
	// for most cases, it is a snp
	des->type = var_type_snp;
	des->start = des->end = pos;
	des->ref_length = des->alt_length = 1;
	des->ref = strndup(r, 1);
	des->alt = strndup(a, 1);
	return des;
    }
    // if alternate allele longer than ref, take it as insertion
    if ( *a && !*r ) {
	des->type = var_type_ins;
	while ( *a ) a++;
	des->start = pos;
	des->ref_length = 0;
	des->ref = str_init;
	des->end = pos +1;
	des->alt_length = (a-alt) - (r-ref);
	des->alt = strndup(a-des->alt_length, des->alt_length);
	return des;
    }
    // if ref allele longer than alt, should be deletion
    if ( !*a && *r) {
	
	des->type = var_type_dels;
	while ( *r ) r++;
	des->start = pos;
	des->ref_length = (r-ref) - (a-alt);
	des->ref = strndup(r-des->ref_length, des->ref_length);
	des->end = pos + des->ref_length -1;
	des->alt_length = 0;
	des->alt = str_init;
	return des;
    }

    // trim tails if ends are same
    const char *ae = a;
    const char *re = r;
    while ( ae[1] ) ae++;
    while ( re[1] ) re++;
    while ( re > r && ae > a && toupper(*re) == toupper(*ae) ) {
	re--;
	ae--;
    }

    // check a and e in first step, so "re==r && ae == a" would not happen here 
    if ( ae == a) {
	des->type = var_type_dels;
	des->start = pos;
	des->ref_length = re -r;
	des->ref = strndup(r, re-r);
	des->end = pos + des->ref_length -1;
	des->alt_length = 0;
	des->alt = str_init;
	return des;
    }
    if (re == r) {
	des->type = var_type_ins;
	des->start = pos;
	des->end = pos+1;
	des->ref_length = 0;
	des->ref = str_init;
	des->alt_length = ae - a;
	des->alt = strndup(a, ae-a);
	return des;
    }
    // delins
    des->type = var_type_delins;
    des->start = pos;
    des->ref_length = re-r;
    des->ref = strndup(r, re-r);
    des->alt_length = ae-a;
    des->alt = strndup(a, ae-a);
    des->end = pos + des->ref_length -1;
    return des;
}
static void check_cnv_hgvs_des(struct hgvs_des *des, struct genepred_line *line, faidx_t *fai, htsFile *fp)
{
    if (des->type == var_type_ins) {

    } else if (des->type == var_type_dels) {

    }
}
static void find_exons_loc(uint32_t pos, int exoncount, struct exon_pair *pair, int *l1, int *l2 )
{
    *l1 = 0;
    *l2 = exoncount*2 -1;
    // init          l1
    //       start   |       |       |
    //       end     |       |       |
    //                               l2
    // result:               l2
    //               |       |       |
    //               |       |       |
    //                       l1    
    while ( *l1 < *l2 ) {
	if ( *l2 - *l1 == 1 ) break;
	int l = *l1 + *l2;
	l = l & 1 ? l/2 + 1 : l/2;
	uint32_t iter = l & 1 ?  pair[l/2].end : pair[l/2].start;
	if ( iter > pos )
	    *l2 = l;
	else
	    *l1 = l;
    } 
}
// convert pos to hgvs pos string
static enum func_region_type pos_convert(int32_t pos, int exoncount, int strand, struct exon_pair *pair, struct exon_offset_pair *locs, int *is_coding, char **cpos)
{
    // purpose: find the most nearest edge for variant position
    // init:   l1                     l2
    //         |----|   |----|  |-----|  
    // End:         |  .|
    //             l1   l2
    // for iter l1 and l2, if even iter come from start array, if odd iter come from end array        

    int l1, l2;    
    find_exons_loc(pos, exoncount, pair, &l1, &l2);
    assert(l2 - l1 == 1);
    
    enum func_region_type ftype = func_region_unknown;
    kstring_t str = KSTRING_INIT; 
    // offset[1,2] are the length between pos and near edges
    uint32_t pos_start, pos_end, loc_start, loc_end;
    uint8_t type_start, type_end;
    if ( l1 & 1 ) {
	pos_start = pair[l1/2].end;
	pos_end = pair[l2/2].end;
	loc_start = locs[l1/2].end >> OFFSET_BITWISE;
	loc_end = locs[l2/2].end >> OFFSET_BITWISE;
	type_start = locs[l1/2].end & REG_MASK;
	type_end =  locs[l2/2].end & REG_MASK;
    } else {
	pos_start = pair[l1/2].start +1;
	pos_end = pair[l2/2].end;
	loc_start = locs[l1/2].start >> OFFSET_BITWISE;
	loc_end = locs[l2/2].end >> OFFSET_BITWISE;
	type_start = locs[l1/2].start & REG_MASK;
	type_end =  locs[l2/2].end & REG_MASK;
    }
    if (pos < pos_start || pos > pos_end) {
	kputc('?', &str);
	*is_coding = 0;
	*cpos = str.s;
	return ftype;
    }

    int32_t offset_start = pos - pos_start;
    int32_t offset_end = pos_end - pos;    
    // if pos in Intron, check the nearest edge
    if ( l1 & 1 ) {	
	// the nestest edge is start of exon l2/2
	uint32_t loc;
	uint8_t type;
	int offset;
	// find the most nearest edge
	// cap to end
	if ( offset_start > offset_end ) {
	    loc = loc_end;
	    type = type_end;
	    offset = strand == '+' ? -offset_end : offset_end;
	} else {
	    loc = loc_start;
	    type = type_start;
	    offset = strand == '+' ? offset_start : -offset_start;
	}

	if ( type & REG_NONCODING ) {
	    *is_coding = 0;
	    kputw(loc, &str);
	} else {
	    *is_coding = 1;
	    if (type & REG_UTR5) {
		kputc('-', &str);
	    } else if (type & REG_UTR3) {
		kputc('*', &str);
	    }
	    kputw(loc, &str);
	}
	if (offset > 0) {
	    kputc('+', &str);
	}
	kputw(offset, &str);		 
	*cpos = str.s;
	return offset_start < 4 || offset_end < 4 ? func_region_split_sites :func_region_intron;
    }
    // if pos in exon
    // pos_[start, end] is the locs of near edges,
    //    pos_start pos       pos_end
    //    |         |         |
    //    |-------------------|
    //    offset_start offset_end
    int loc;
    uint8_t type;
    if ( type_start == type_end ) {
	loc = loc_end > loc_start ?  loc_end - offset_end : loc_start - offset_start;
	type = type_start;
	goto generate_cpos;
    }
    // edge stituation. if start and end come to different function regions.
    if ( !(type_start & (REG_UTR3 | REG_UTR5 | REG_CODING)) )
	error("Unknown type. type_start : %d", type_start);
    if ( !(type_end & (REG_UTR3 | REG_UTR5 | REG_CODING)) )	
	error("Unknown type. type_end : %d", type_end);
    
    if ( type_start & REG_UTR5) {
	// should only be plus strand
	loc = loc_start - offset_start;
	type = type_start;
	if ( type_end & REG_CODING ) {
	    if ( loc <= 0 ) {
		loc = -loc +1;
		type = type_end;
	    }
	} else if ( type_end == REG_UTR3 ) {
	    // the cds region inside one exon, closed with UTRs
	    if (loc <= 0) {
		int loc1 = loc_end - offset_end;
		if (loc1 <= 0) {
		    loc = -loc + 1;
		    type = REG_CODING;
		} else {
		    loc = loc1;
		    type = type_end;
		}			    
	    }		    
	}
	goto generate_cpos;
    }
    if ( type_start & REG_CODING ) {	    
	if ( type_end & REG_UTR3) {
	    // strand plus
	    loc = loc_end - offset_end;
	    type = type_end;
	    if ( loc <= 0 ) {
		loc = loc_start + offset_start;
		type = type_start;
	    }
	} else if (type_end & REG_UTR5) {
	    // strand minus
	    loc = loc_start - offset_start;
	    type = type_start;
	    if ( loc <= 0 ) {
		loc = -loc + 1;
		loc = type_end;
	    }
	}
	goto generate_cpos;
    }
    
    if (type_start == REG_UTR3) {
	// should only be minus strand
	if ( type_end == REG_CODING ) {
	    loc = loc_start - offset_start;
	    type = type_start;
	    if (loc <= 0) {
		loc = -loc + 1;
		type = type_end;
	    }		    
	} else if (type_end == REG_UTR5) {
	    // cds region inside one exon
	    loc = loc_start - offset_start;
	    type = type_start;
	    if ( loc <= 0 ) {
		int loc1 = loc_end - offset_end;
		if (loc1 <= 0) {
		    loc = -loc + 1;
		    type = REG_CODING;			    
		} else {
		    loc = loc1;
		    type = type_end;
		}
	    }
	}
	goto generate_cpos;
    }
    // if not found the type of start 
    error("This is a impossible stituation. type_start : %d", type_start);

  generate_cpos:
    assert(loc > 0);    
    ftype = func_region_cds;	
    if (type & REG_NONCODING) {
	*is_coding = 0;
    } else {
	*is_coding = 1;
	if (type & REG_UTR5) {
	    kputc('-', &str);
	    ftype = func_region_utr5;
	} else if (type & REG_UTR3) {
	    kputc('*', &str);
	    ftype = func_region_utr3;
	}
    }
    kputw(loc, &str);
    *cpos = str.s;
    return ftype;
}
static void generate_hgvs_core(struct genepred_line *line, struct hgvs_des *des, struct hgvs_core *c, int *valid)
{
    hgvs_core_clear(c);
    *valid = 0;
    if (des->start > line->txend || des->end < line->txstart)
	return;    
    uint32_t start_var = (uint32_t)des->start;
    uint32_t end_var = (uint32_t)des->end;
    kstring_t *string = &c->str;
    // put names into string, usually transcript 
    // if no transcript name, use gene name then
    char *name = line->name1 == NULL || line->name1[0] == '.' ? line->name2 : line->name1;    
    kputs(name, string);
    c->l_name = string->l;    
    assert(c->l_name);
    kputc(':', string);

    // put hgvs pos into string
    char *pos1;
    int is_coding = -1;
    c->type = pos_convert(des->start, line->exoncount, line->strand, line->exons, line->dna_ref_offsets, &is_coding, &pos1);
    if (is_coding == 1)
	kputs("c.", string);
    else if (is_coding == 0)
	kputs("n.", string);
    else
	error("Unknown transcript type! %s", name);
    kputs(pos1, string);
    // put offsets into string
    c->l_type = string->l -1;
    copy_seqs_func func = line->strand == '+' ? strndup : rev_seqs;
    char *ref = des->ref_length == 0 ? NULL : func(des->ref, des->ref_length);
    char *alt = des->alt_length == 0 ? NULL : func(des->alt, des->alt_length);
    int is_coding1;
    if (des->type == var_type_snp) {
	ksprintf(string, "%s>%s", ref, alt);
    } else if (des->type == var_type_dels) {
	// for one base deletion, format like NM_0001:c.123del
	if (des->ref_length == 1) {
	    ksprintf(string, "del%s", des->ref);
	} else {
	    // for two bases or more, format like NM_0001:c.123_125del
	    kputc('_', string);
	    char *pos2;
	    pos_convert(des->end, line->exoncount, line->strand, line->exons, line->dna_ref_offsets, &is_coding1, &pos2);
	    kputs(pos2, string);
	    free(pos2);
	    kputs("del", string);
	}
    } else if (des->type == var_type_ins) {
	// insertion format should be NM_0001:c.123_124insXX
	kputc('_', string);
	char *pos2;
	pos_convert(des->end, line->exoncount, line->strand, line->exons, line->dna_ref_offsets, &is_coding1, &pos2);
	kputs(pos2, string);
	free(pos2);
	// if insertion sequences is small show sequences directly. for long sequences show number.
	if ( des->alt_length < 20)
	    ksprintf(string, "ins%s", alt);
	else
	    ksprintf(string, "ins%d", des->alt_length);
    } else if (des->type == var_type_delins) {
	// delins format should be NM_0001:c.123_125delinsXX
	if ( des->ref_length > 1 ) {
	    kputc('_', string);
	    char *pos2;
	    pos_convert(des->end, line->exoncount, line->strand, line->exons, line->dna_ref_offsets, &is_coding1, &pos2);
	    kputs(pos2, string);
	    free(pos2);
	}
	if ( des->alt_length < 20 )
	    ksprintf(string, "delins%s", alt);
	else
	    ksprintf(string, "delins%d", des->alt_length);
    } else {
	error("unknow type : %d", des->type);
    }
    
#ifdef DEBUG_MODE
    debug_print("%d : %s", des->start,string->s);
#endif
    free(pos1);
    // it is valid
    *valid = 1;    
    if ( des->ref_length != 0 )
	free(ref);
    if ( des->alt_length != 0 )
	free(alt);    
}
void name_list_push(struct name_list *names, char *name)
{
    int i;
    for ( i = 0; i < names->l; ++i ) {
	if ( strcmp(names->a[i], name) == 0 )
	    return;
    }
    if ( i == names->l ) {
	if ( names->m == names->l ) {
	    names->m = names->m == 0 ? 2 : names->m<<1;
	    names->a = (char**)realloc(names->a, sizeof(char*)*names->m);
	}
	names->a[names->l++] = (char*)strdup(name);
    }
}
// -- assuming type of hgvs.name tag is string and number is R. this is mandontary for VCF file 
// @pool        cached memory pool for gene predictions records
// @line        input variant record
// @n_names1    number of gene names, usually all the transcripts overlaped in this regions only
// belong to one gene, but for structure variants, this variants may span more than
// one. Also, there might be some genes located in the introns of another genes. It
// is rarely happened, but we should consider of it.
// @names1      name of gene names
// @n_ales      related allele numbers (R)
// @n_hgvs      number of hgvs names strings, several transcript for one allele 
// should be cache into one string like NM_xxxx:c.xx; NC_xxxx:n.xxx; ...; 
// @hgvs  name of hgvs names strings
// @n_names2    number of transcripts
// @names2      name of transcripts
void generate_hgvs(struct genepred_memory_pool *pool, const bcf1_t *line, int *n_names2, char **names2, int *n_ale)
{
    *n_names2 = 0;
    *n_ale = line->n_allele -1;
    if (*n_ale == 0)
	return;
    //bcf_dec_t *d = &line->d;
    struct hgvs_cache * cache = &pool->cache;
    // simple name structure inited for caching gene names and de-duplicate 
    // struct  names = KSTRING_INIT;
    if ( cache->m < line->n_allele ) {
	cache->m = line->n_allele;
	cache->a = (struct hgvs*)realloc(cache->a, sizeof(struct hgvs)*cache->m);
    }
    cache->l = 0;
    // roadmap : snv -> del -> ins -> dup ->delins
    struct name_list names = { 0, 0, 0};
    int i, j;    
    for ( i = 1; i < line->n_allele; ++i ) {
	struct hgvs *name = &cache->a[cache->l];
	if (cache->l == cache->i) {
	    name->l = name->m = name->i = 0;
	    name->a = NULL;
	    cache->i++;
	}
	cache->l++;
	hgvs_clear(name);
	struct hgvs_des *des = describe_variants(line->d.allele[0], line->d.allele[i], line->pos +1);
	for ( j = 0; j < pool->l; ++j ) {
	    int valid = 0;
	    if (name->l == name->m) {
		name->m = name->m == 0 ? 2 : name->m<<1;
		name->a = (struct hgvs_core*)realloc(name->a, sizeof(struct hgvs_core)*name->m);
	    }
	    if (name->i == name->l) {
		kstring_t *string = &name->a[name->i].str;
		string->l = string->m = 0;
		string->s = NULL;
		name->i++;
	    }
	    struct hgvs_core *core= &name->a[name->l];
	    struct genepred_line *gl = &pool->a[j];
	    hgvs_des_reverse(des, gl->strand);
	    generate_hgvs_core(gl, des, core, &valid);
	    if (valid == 0)
		continue;
	    name_list_push(&names, gl->name2);
	    name->l++;
	}
	hgvs_des_destory(des);
	
    }

    kstring_t str = KSTRING_INIT;
    int n = 0;
    for (i = 0; i < names.l; i++) {
	if ( i ) kputc(',', &str);
	kputs(names.a[i], &str);
	free(names.a[i]);
    }
    if ( names.m )
	free(names.a);
    *n_names2 = names.l;
    *names2 = str.l ? str.s : NULL;
}
// number : R, type : string
static void setter1_hgvs(bcf_hdr_t *hdr, bcf1_t *line, char *string) 
{
    // the data construct along with alt alleles, so there is no need to remap the order of alleles
    bcf_update_info_string(hdr, line, "HGVSDNA", string);
}
// number : 1, type : string
static void setter1_gene_names(bcf_hdr_t *hdr, bcf1_t *line, char *string)
{
    bcf_update_info_string(hdr, line, "Gene", string);
}
// number : R, type : string
static void setter1_transcripts_names(bcf_hdr_t *hdr, bcf1_t *line, char *string)
{
    bcf_update_info_string(hdr, line, "Transcript", string);
}
void setter_hgvs_string(bcf_hdr_t *hdr, bcf1_t *line, const char *key, char *string)
{
    if (key == NULL) return;
    if (strcmp(key, "Gene") == 0) {
	setter1_gene_names(hdr, line, string);	
    } else if (strcmp(key, "HGVSDNA") == 0) {
	setter1_hgvs(hdr, line, string);		
    } else if (strcmp(key, "Transcript") == 0) {
	setter1_transcripts_names(hdr, line, string);
    } else {
	error("Unrecongnized tag %s, only Gene, HGVSDNA, and Transcript are supported for now.", key);
    }    
}
char *generate_transcript_string(struct hgvs_cache *cache, int n)
{
    kstring_t str = KSTRING_INIT;
    int i, j;    
    for ( i = 0; i < n; ++i ) {
	if ( i )
	    kputc(',', &str);
	// foreach allele
	struct hgvs *name = &cache->a[i];
	for ( j = 0; j < name->l; ++j ) {
	    struct hgvs_core *core = &name->a[j];
	    if (core->l_name)
		kputsn(core->str.s, core->l_name-1, &str);
	    else
		kputc('.', &str);
	    kputc('|', &str);
	}
    }
    return str.s;
}
char *generate_hgvsvarnomen_string(struct hgvs_cache *cache,int n)
{
    kstring_t str = KSTRING_INIT;
    int i, j;    
    for ( i = 0; i < n; ++i ) {
	if ( i )
	    kputc(',', &str);
	// foreach allele
	struct hgvs *name = &cache->a[i];
	for ( j = 0; j < name->l; ++j ) {
	    struct hgvs_core *core = &name->a[j];
	    if (core->l_name)
		kputs(core->str.s, &str);
	    else
		kputc('.', &str);
	    kputc('|', &str);
	}
    }
    return str.s;
}

// hgvs_anno() only used to annotate bcf/vcf standalone, all the valid tags-added functions will be called in 
// this function, to annotate the vcf more specification, use setter_genepred_* functions instead of it. 
void hgvs_anno(struct genepred_memory_pool *pool, htsFile *fp, tbx_t *tbx, bcf_hdr_t *hdr, bcf1_t *line, struct genepred_format *type)
{
    // retrieve the regions this variants located first, check the memory pool and update the pool if the regions is out of
    // cached positions. just skip if the variant type of line is a ref. 
    if (line->pos <= 0)
	return;
    if (bcf_get_variant_types(line) == VCF_REF)
	return;

    uint32_t end = bcf_calend(line);
    int id = tbx_name2id(tbx, bcf_hdr_id2name(hdr, line->rid));
    if ( pool->rid == -1 || pool->rid != line->rid || pool->start > line->pos+1 || pool->end < line->pos+1 ) {
	// fill mempool, skip if no record in the memory pool
	pool->l = 0;
	hgvs_fill_memory(pool, fp, tbx, id, line->pos, end, type);
	if (pool->l == 0) return;
    }
    int n_names2 = 0;
    char *names2 = NULL;
    int n_ale = 0;
    
    generate_hgvs(pool, line, &n_names2, &names2, &n_ale);
    if (n_ale == 0)
	return;
    char *trans_string = generate_transcript_string(&pool->cache, n_ale);
    char *hgvs_string = generate_hgvsvarnomen_string(&pool->cache, n_ale);
#ifdef DEBUG_MODE
    debug_print("trans : %s", trans_string);
    debug_print("hgvs : %s", hgvs_string);
#endif
    if ( hgvs_string )
    	setter_hgvs_string(hdr, line, "Transcript", trans_string);
    if ( trans_string )
    	setter_hgvs_string(hdr, line, "HGVSDNA", hgvs_string);
    if ( n_names2 )
	setter_hgvs_string(hdr, line, "Gene", names2);
    free(names2);
    free(trans_string);
    free(hgvs_string);
}

#ifdef _HGVS_MAIN

#include <sys/time.h>
#include <sys/stat.h>

static double get_time() {
    struct timeval tv;
    if ( gettimeofday(&tv, 0) != 0 ) 
	error("Failed to get time of day.");
    return (double)tv.tv_sec + 1.0e-6 *(double)tv.tv_usec;
}

struct args {
    const char *genepred_file;
    const char *refseq_file;
    const char *out_file;
    const char *out_type;
    const char *input_fname;
    htsFile *fp;
    htsFile *fout;
    tbx_t *gpidx;
    faidx_t *faidx;
    namehash_type *gname;
    struct genepred_memory_pool pool;
};
struct args args = {
    .genepred_file = 0,
    .refseq_file = 0,
    .out_file = 0,
    .out_type = 0,
    .input_fname = 0,
    .gpidx = 0,
    .faidx = 0,
    .gname = 0,
};

void clear_args(struct args *args)
{
    tbx_destroy(args->gpidx);
    if (args->faidx) fai_destroy(args->faidx);
    if (args->gname) kh_destroy(name, args->gname);
    hgvs_memory_clear(&args->pool);
    hts_close(args->fout);
    hts_close(args->fp);
}
int usage(int argc, char **argv)
{
    fprintf(stderr, "Usage: %s [-h] -refseq refseq.fa.gz -data gene_pred.txt.gz|refgene.txt.gz -O [u|v|z|b] -o output_file input_fname.vcf.gz \n", argv[0]);
    return 0;
}
const char *hts_bcf_wmode(int file_type)
{
    if ( file_type == FT_BCF ) return "wbu";    // uncompressed BCF
    if ( file_type & FT_BCF ) return "wb";      // compressed BCF
    if ( file_type & FT_GZ ) return "wz";       // compressed VCF
    return "w";                                 // uncompressed VCF
}

int main(int argc, char **argv)
{
    int i;        
    for (i = 1; i< argc;) {
	const char *a = argv[i++];
	if (strcmp(a, "-h") == 0) 
	    return usage(argc, argv);

	const char **arg_var = 0;
	if (strcmp(a, "-o") == 0 && args.out_file == 0)
	    arg_var = &args.out_file;
	else if (strcmp(a, "-O") == 0 && args.out_type == 0)
	    arg_var = &args.out_type;
	else if (strcmp(a, "-data") == 0 && args.genepred_file == 0)
	    arg_var = &args.genepred_file;
	else if (strcmp(a, "-refseq") == 0 && args.refseq_file == 0)
	    arg_var = &args.refseq_file;

	if (arg_var != 0) {
	    if (i == argc)
		error("Missing arg after %s", a);
	    *arg_var = argv[i++];
	    continue;
	}

	if (args.input_fname == 0) {
	    args.input_fname = a;
	    continue;
	}
	error("Unknown arg : %s", a);
    }
    // assuming input file is stdin, only vcf, gzipped vcf, bcf file is accept,
    // err msg will output if file type unrecongnized
    if (args.input_fname == 0 && (!isatty(fileno(stdin)))) args.input_fname = "-";
    if (args.input_fname == 0)
	error("No input file ! Use -h for more informations.");
 
    if (args.genepred_file == 0)
	error("Reasons :\n"
	      "-data gene preditions database is needed.\n"
	      "You could download refGene.txt.gz or ensGene.txt.gz from UCSC websites and sort and indexed by tabix.");


    htsFile *fgp = hts_open(args.genepred_file, "r");
    if (fgp == NULL)
	error("Failed to open %s", args.genepred_file);
    args.gpidx = load_genepred_file(args.genepred_file);
    if (args.gpidx == NULL)
	error("Failed to load index of %s", args.genepred_file);

    // args.faidx = load_refseq_file(args.refseq_file); 
    // if (args.faidx == NULL) 
    // 	error("Failed to load index of %s", args.refseq_file); 
    
    args.fp = hts_open(args.input_fname, "r");
    if (args.fp == NULL)
	error("Failed to open %s.", args.input_fname);
    htsFormat type = *hts_get_format(args.fp);
    if (type.format  != vcf && type.format != bcf)
	error("Unsupported input format! %s", args.input_fname);
    int out_type = FT_VCF;
    if (args.out_type != 0) {
	switch (args.out_type[0]) {
	    case 'b':
		out_type = FT_BCF_GZ; break;
	    case 'u':
		out_type = FT_BCF; break;
	    case 'z':
		out_type = FT_VCF_GZ; break;
	    case 'v':
		out_type = FT_VCF; break;
	    default :
		error("The output type \"%d\" not recognised\n", out_type);
	};
    }

    // init gene predictions memory pool 
    struct genepred_memory_pool *pool = &args.pool;
    pool->m = pool->l = 0;
    pool->rid = -1;
    pool->start = pool->end = 0;
    pool->a = NULL;
    pool->cache.i = pool->cache.l = pool->cache.m = 0;
    pool->cache.a = NULL;
    double c0 = get_time();
    LOG_print("Init ...");
    bcf_hdr_t *hdr = bcf_hdr_read(args.fp);    
    LOG_print("Prase header.");
    // duplicate header file and sync new tag in the output header .
    // assuming output is stdout in Vcf format in default. 
    args.fout = args.out_file == 0 ? hts_open("-", hts_bcf_wmode(out_type)) : hts_open(args.out_file, hts_bcf_wmode(out_type));
    bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);	
    hgvs_bcf_header_add_gene(hdr);
    hgvs_bcf_header_add_trans(hdr);
    hgvs_bcf_header_add_dna(hdr);
#ifndef DEBUG_MODE   
    bcf_hdr_write(args.fout, hdr_out);
#endif
    // init gene_pred or refgene database, hold tabix index cache and faidx cache in memory 
    bcf1_t *line = bcf_init();
    while ( bcf_read(args.fp, hdr, line) == 0 ) {     
	hgvs_anno(pool, fgp, args.gpidx, hdr_out, line, &genepred_formats);
#ifndef DEBUG_MODE   
	bcf_write1(args.fout, hdr_out, line);
#endif
    }
    LOG_print("Clear ...");
    clear_args(&args);
    hts_close(fgp);
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    bcf_hdr_destroy(hdr_out);
    double c1 = get_time();
    fprintf(stderr, "Run time: %.2fs\n", c1 -c0);
    return 0;
}
#endif
