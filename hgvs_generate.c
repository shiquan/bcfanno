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
#include "anno_hgvs.h"

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
static void hgvs_cache_clear(struct hgvs_cache *cache)
{
    int i;
    for ( i = 0; i < cache->l; ++i )
	hgvs_clear(&cache->a[i]);
    cache->l = 0;
}
static void hgvs_cache_destroy(struct hgvs_cache *cache)
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
 }
static namehash_type *init_gene_name(const char *name)
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
static tbx_t *load_genepred_file(const char *fn)
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

static struct genepred_format *type = &genepred_formats;

void set_format_refgene()
{
    type = &refgene_formats;
}
void set_format_genepred()
{
    type = &genepred_formats;
}
static void genepred_prase_core(kstring_t *string, struct genepred_line *line)
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
void genepred_line_praser(kstring_t *string, struct genepred_line *line)
{
    // prase string into genepred line struct
    genepred_prase_core(string, line);
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
    	    int t = types[j];
    	    switch ( t ) {
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
    		    error("Unknown type : %d", t);
    	    }
    	}
    	kputw((line->dna_ref_offsets[i].start>>OFFSET_BITWISE), &temp[0]);
    	kputw((line->dna_ref_offsets[i].end>>OFFSET_BITWISE), &temp[1]);

	// format: CHROM,START,END,STRAND, GENE, TRANSCRIPT, EXON, START_LOC, END_LOC 
    	free(temp[0].s);
    	free(temp[1].s);
    }
}
static void push_mempool(struct genepred_memory *pool, kstring_t *str)
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
    genepred_line_praser(str, &pool->a[pool->l]);
    pool->l++;
}
static void update_mempool(struct genepred_memory *pool)
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
static int hgvs_fill_memory(struct refgene_options *opts, bcf1_t *line)
{
    // get the chromosome id in tabix indexed genepred database
    int id = tbx_name2id(opts->genepred_idx, bcf_hdr_id2name(opts->hdr_out, line->rid));
    // no chromosome found, return 1
    if (id == -1)
	return 1;    

    int32_t end = bcf_calend(line);    
    hts_itr_t *itr = tbx_itr_queryi(opts->genepred_idx, id, line->pos+1, end);
    struct genepred_memory *pool = &opts->buffer;
    kstring_t string = KSTRING_INIT;
    
    while ( tbx_itr_next(opts->fp, opts->genepred_idx, itr, &string) >= 0) {
	push_mempool(pool, &string);
	string.l = 0;
    }
    pool->rid = id;
    if ( string.m )
	free(string.s);
    tbx_itr_destroy(itr);
    update_mempool(pool);
    return 0;
}
static void hgvs_memory_clear(struct genepred_memory *pool)
{
    int i;
    for ( i = 0; i < pool->i; ++i)
	genepred_line_clear(&pool->a[i]);
    
    pool->l = 0;
    // i and m should be kept for future use
    pool->start = 0;
    pool->end = 0;
    pool->rid = -1;
}
static void namehash_destroy(void *_hash)
{
    namehash_type *hash = (namehash_type*)_hash;
    khiter_t k;
    for ( k = 0; k < kh_end(hash); ++k ) {
	if ( kh_exist(hash, k) ) 
	    free((char*)kh_key(hash, k));
    }
    kh_destroy(name, hash);
}

void refgene_opts_destroy(struct refgene_options *opt)
{
    hts_close(opt->fp);
    tbx_destroy(opt->genepred_idx);
    if ( opt->check_refseq == 1 )
	fai_destroy(opt->refseq_fai);
    if ( opt->screen_by_genes == 1 )
	namehash_destroy(opt->genehash);
    if ( opt->screen_by_transcripts == 1 )
	namehash_destroy(opt->transhash);
}
static void hgvs_des_destory(struct hgvs_des *des)
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
static void hgvs_des_reverse(struct hgvs_des *des, uint8_t strand)
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
    error("This is an impossible stituation. type_start : %d", type_start);

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
static int generate_hgvs_core(struct genepred_line *line, struct hgvs_des *des, struct hgvs_core *c)
{
    if (des->start > line->txend || des->end < line->txstart)
	error("Variants is not in this transcript. des->start: %u, des->end: %d, genepred_line: %s\t%s\t%u\t%u",
	      des->start, des->end, line->chrom, line->txstart, line->txend);
    // clear hgvs_core memory before init
    hgvs_core_clear(c);
    
    uint32_t start_var = (uint32_t)des->start;
    uint32_t end_var = (uint32_t)des->end;
    kstring_t *string = &c->str;
    // put names into string, usually transcript 
    // if no transcript name, use gene name then
    char *name = line->name1;
    char *name2 = line->name2;
    if (name != NULL)
	kputs(name, string);
    c->l_name = string->l;
    if (name2 == NULL) {
	c->l_name2 = 0;
    } else {
	kputc('(', string);
	kputs(name2, string);
	kputc(')', string);
    }
    if (c->l_name == 0) 
	error ("Failed to prase genepred database. %s\t%u\t%u", line->chrom, line->txstart, line->txend);
    kputc(':', string);

    // put hgvs pos into string
    char *pos1;
    int is_coding = -1;
    c->type = pos_convert(des->start, line->exoncount, line->strand, line->exons, line->dna_ref_offsets, &is_coding, &pos1);
    // put the type of transcript
    if (is_coding == 1)
	kputs("c.", string);
    else if (is_coding == 0)
	kputs("n.", string);
    else
	error("Unknown transcript type! %s in %s\t%u\t%u", name, line->chrom, line->txstart, line->txend);    
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
    if ( des->ref_length != 0 )
	free(ref);
    if ( des->alt_length != 0 )
	free(alt);
    return 0;
}

static void name_list_push(struct name_list *names, char *name)
{
    int i;
    for ( i = 0; i < names->l; ++i ) {
	if ( strcmp(names->a[i], name) == 0 )
	    return;
    }
    if ( names->m == names->l ) {
	names->m = names->m == 0 ? 2 : names->m<<1;
	names->a = (char**)realloc(names->a, sizeof(char*)*names->m);
    }
    names->a[names->l++] = (char*)strdup(name);
}
// -- assuming type of hgvs.name tag is string and number is R. this is mandontary for VCF file 
// genepred_memory        cached memory pool for gene predictions records
// line        input variant record in bcf1_t format
// hgvs_cache will init or update in this function
static void generate_hgvs(struct refgene_options *opts, const bcf1_t *line)
{
    struct genepred_memory *buffer = &opts->buffer;
    struct hgvs_cache *cache = &opts->cache;
    int n_alleles = line->n_allele -1;
    if (n_alleles == 0)
	error("This is a reference position, should not come here. %s:%d", bcf_seqname(opts->hdr_out, line), line->pos+1);
    
    // simple name structure inited for caching gene names and de-duplicate 
    // struct  names = KSTRING_INIT;
    if ( cache->m < n_alleles ) {
	cache->m = n_alleles;
	cache->a = (struct hgvs*)realloc(cache->a, sizeof(struct hgvs)*cache->m);
    }
    cache->l = 0;
    // roadmap : snv -> del -> ins -> dup ->delins
    int i, j;    
    for ( i = 1; i < line->n_allele; ++i ) {
	// init hgvs name structure for each allele
	struct hgvs *name = &cache->a[cache->l];
	// if cache is not inited, init and increase i flag
	if (cache->l == cache->i) {
	    name->l = name->m = name->i = 0;
	    name->a = NULL;
	    cache->i++;
	}
	cache->l++;
	// clear name structure
	hgvs_clear(name);

	// get variant description
	struct hgvs_des *des;
	des = describe_variants(line->d.allele[0], line->d.allele[i], line->pos +1);
	// loop genepred memory pool and generate hgvs names
	for ( j = 0; j < buffer->l; ++j ) {

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
	    struct genepred_line *gl = &buffer->a[j];
	    // check strand of transcript and reverse the description if it is minus strand
	    hgvs_des_reverse(des, gl->strand);
	    // generate hgvs name core
	    generate_hgvs_core(gl, des, core);	    
	}
	// destroy hgvs descriptions structure
	hgvs_des_destory(des);
    }
    
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
static void setter_hgvs_string(bcf_hdr_t *hdr, bcf1_t *line, const char *key, char *string)
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
static char *generate_transcript_string(struct hgvs_cache *cache)
{
    kstring_t str = KSTRING_INIT;
    int i, j;    
    for ( i = 0; i < cache->l; ++i ) {
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
static char *generate_hgvsvarnomen_string(struct hgvs_cache *cache)
{
    kstring_t str = KSTRING_INIT;
    int i, j;    
    for ( i = 0; i < cache->l; ++i ) {
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
static char *get_gene_string(struct hgvs_cache *cache)
{
    int i, j;
    kstring_t str = KSTRING_INIT;
    for (i = 0; i < cache->l; ++i) {
	struct hgvs *hgvs = &cache->a[i];
	for (j = 0; j < hgvs->l; ++j) {
	    struct hgvs_core *core = &hgvs->a[j];
	    if (core->l_name  && core->l_name2 > 0) {
		kstring_t *string = &core->str;
		assert(string->s[core->l_name]  == '(');
		kputsn(string->s + core->l_name + 1, core->l_name2 - core->l_name -1, &str);
		return str.s;
	    }
	}
    }
    return NULL;
}
int hgvs_names_init(struct refgene_options *opts, bcf1_t *line)
{
    struct hgvs_cache *cache = &opts->cache;
    // check the hgvs names inited already ?
    hgvs_cache_clear(cache);
    // check the buffer is filled or not

    // fill memory pool, only if no transcriptc cover this variant return 1, else return 0
    if ( hgvs_fill_memory(opts, line) != 0 )
	return 1;
    //
    generate_hgvs(line, opts);
    
    return 0;
}

int anno_refgene_core(bcf1_t *line, struct refgene_options *opts)
{
    // retrieve the regions this variants located first, check the memory pool and update the pool if the regions is out of
    // cached positions. just skip if the variant type of line is a ref.     
    if ( line->pos < 0 )
	return 1;
    // skip reference positions, usually vcfanno will check the type of position in the first step, get error if see ref here
    if ( bcf_get_variant_types(line) == VCF_REF )
	error("This is a reference position. Should not come here. %s : %d", bcf_seqname(opts->hdr_out, line), line->pos+1);
    // fill memory pool and init hgvs names structure, return 1 if no transcript found
    if ( hgvs_names_init(opts, line) != 0 )
	return 1;
    // setter each column into vcf line
    int i;
    for ( i = 0; i < opts->n_cols; ++i ) {
	struct anno_col *col = &opts->cols[i];
	col->setter->hgvs(opts, line, col);
    }
    // if success, return 0
    return 0;
}
static char *hgvs_get_names_string(struct hgvs_cache *cache, const char *key)
{
    if ( strcmp(key, "Gene") == 0 ) {
	 return get_gene_string(cache);
    } else if ( strcmp(key, "Transcript") == 0 ) {
	return generate_transcript_string(cache);
    } else if ( strcmp(key, "HGVSDNA") == 0) {
	return generate_hgvsvarnomen_string(cache);
    } else {
	error("Failed to recongnise key %s", key);
    }
}
int setter_hgvs_names(struct refgene_options *opts, bcf1_t *line, struct anno_col *col)
{
    uint32_t end = bcf_calend(line);
    // hgvs_names_init(opts, line);
    char *string = hgvs_get_names_string(&opts->cache, col->hdr_key);
    setter_hgvs_string(opts->hdr_out, line, col->hdr_key, string);
    return 0;
}
int hgvs_cols_prase(struct refgene_options *opts, const char *column)
{
    kstring_t string = KSTRING_INIT;
    kputs(column, &string);
    int nfields = 0;
    int *splits = ksplit(&string, ',', &nfields);
    opts->n_cols = nfields;
    opts->cols = (struct anno_col*)calloc(opts->n_cols, sizeof(struct anno_col));
    int i;
    for (i = 0; i < opts->n_cols; ++i ) {
	struct anno_col *col = &opts->cols[i];
	char *ss = string.s + splits[i];
	if (ss[0] == '+') {
	    col->replace = REPLACE_MISSING;
	    ss++;
	} else if (ss[0] == '-') {
	    col->replace = REPLACE_EXISTING;
	    ss++;
	} else {
	    col->replace = REPLACE_ALL;
	}
	col->hdr_key = strdup(ss);	
	if ( strcmp(col->hdr_key, "Gene") == 0 ) {
	    col->icol = hgvs_bcf_header_add_gene(opts->hdr_out);	    
	} else if ( strcmp(col->hdr_key, "Transcript") == 0 ) {
	    col->icol = hgvs_bcf_header_add_trans(opts->hdr_out);	    
	} else if ( strcmp(col->hdr_key, "HGVSDNA") == 0 ) {
	    col->icol = hgvs_bcf_header_add_hgvs(opts->hdr_out);
	} else {
	    error("Failed to prase column tag %s", col->hdr_key);
	}
	col->number =  bcf_hdr_id2length(opts->hdr_out, BCF_HL_INFO, col->icol);
    }
    return 0;
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

struct refgene_options opts = {
    .genepred_fname = 0,
    .refseq_fname = 0,
    .fp = 0,
    .check_refseq = 0,
    .hdr_out = 0,
    .trans_list_fname = 0,
    .genes_list_fname = 0,
    .n_cols = 0,
    .cols = 0,
    .buffer = GENEPRED_MEMORY_INIT,
    .cache = HGVS_CACHE_INIT,    
};
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
    const char *out_type = 0;
    const char *out_fname = 0;
    const char *columns = "Gene,Transcript,HGVSDNA";
    for (i = 1; i< argc;) {
	const char *a = argv[i++];
	if (strcmp(a, "-h") == 0) 
	    return usage(argc, argv);

	const char **arg_var = 0;
	if (strcmp(a, "-o") == 0 && optsout_file == 0)
	    arg_var = &out_fname;
	else if (strcmp(a, "-O") == 0 && args.out_type == 0)
	    arg_var = &out_type;
	else if (strcmp(a, "-data") == 0 && args.genepred_file == 0)
	    arg_var = &opts.genepred_fname;
	else if (strcmp(a, "-refseq") == 0 && args.refseq_file == 0)
	    arg_var = &opts.refseq_fname;

	if (arg_var != 0) {
	    if (i == argc)
		error("Missing arg after %s", a);
	    *arg_var = argv[i++];
	    continue;
	}

	if (input_fname == 0) {
	    input_fname = a;
	    continue;
	}
	error("Unknown arg : %s", a);
    }
    // assuming input file is stdin, only vcf, gzipped vcf, bcf file is accept,
    // err msg will output if file type unrecongnized
    if (input_fname == 0 && (!isatty(fileno(stdin))))
	input_fname = "-";
    if (input_fname == 0)
	error("No input file ! Use -h for more informations.");
 
    if (opts.genepred_file == 0)
	error("Reasons :\n"
	      "-data gene preditions database is needed.\n"
	      "You could download refGene.txt.gz or ensGene.txt.gz from UCSC websites and sort and indexed by tabix.");

    hgvs_cols_prase(&opts, columns);
    opts.fp = hts_open(opts.genepred_fname, "r");
    if (opts.fp == NULL)
	error("Failed to open %s", opts.genepred_fname);
    opts.genepred_idx = load_genepred_file(opts.genepred_fname);
    if (opts.genepred_idx == NULL)
	error("Failed to load index of %s", opts.genepred_fname);
    
    htsFile *fp = hts_open(input_fname, "r");
    if (fp == NULL)
	error("Failed to open %s.", input_fname);
    htsFormat type = *hts_get_format(fp);
    if (type.format  != vcf && type.format != bcf)
	error("Unsupported input format! %s", input_fname);
    int out_type = FT_VCF;
    if (out_type != 0) {
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

    double c0 = get_time();
    LOG_print("Init ...");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);    
    LOG_print("Prase header.");
    // duplicate header file and sync new tag in the output header .
    // assuming output is stdout in Vcf format in default. 
    htsFile *fout = out_fname == 0 ? hts_open("-", hts_bcf_wmode(out_type)) : hts_open(out_fname, hts_bcf_wmode(out_type));
    opts.hdr_out = bcf_hdr_dup(hdr); 

    hgvs_bcf_header_add_trans(hdr);
    hgvs_bcf_header_add_dna(hdr);
#ifndef DEBUG_MODE   
    bcf_hdr_write(args.fout, hdr_out);
#endif
    set_format_genepred();
    // init gene_pred or refgene database, hold tabix index cache and faidx cache in memory 
    bcf1_t *line = bcf_init();
    while ( bcf_read(fp, hdr, line) == 0 ) {     
	anno_refgene_core(line, &opts);
#ifndef DEBUG_MODE   
	bcf_write1(fout, hdr_out, line);
#endif
    }
    LOG_print("Clear ...");
    hts_close(opts.genepred_idx);
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    bcf_hdr_destroy(hdr_out);
    refgene_opts_destroy(&opts);
    double c1 = get_time();
    fprintf(stderr, "Run time: %.2fs\n", c1 -c0);
    return 0;
}
#endif
