#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <zlib.h>
#include "utils.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/khash.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "htslib/faidx.h"

#define KSTRING_INIT {0, 0, 0}

char str_init[2];

KHASH_MAP_INIT_STR(list, char*)

typedef khash_t(list) glist_t;

typedef char *(*copy_seqs_func)(const char*, unsigned long n);
		      
//  for each transcript, the region span several exon partial regions, use exon_pair[] stand each exons
// and exon_offset_pair[] is the offset / coding coordiante consider of cdsstart 
struct exon_pair {
    uint32_t start;
    uint32_t end;
};

// 
// The dna reference offset for DNA coding and noncoding transcript consist of two parts. the offset value
// and the function region construct a 32-bits value.
// coding/nocoding coordinate
// 32                        4321
// |_________________________||||
// ONLY one-bit below is accept per value
// 
// offset 1:  is_noncoding
// offset 2:  is_coding
// offset 3:  is_utr5
// offset 4:  is_utr3
// 

#define OFFSET_BITWISE   4
#define REG_NONCODING  1
#define REG_CODING     2
#define REG_UTR5       4
#define REG_UTR3       8
#define REG_MASK       0xF

struct exon_offset_pair {
    uint32_t start;
    uint32_t end;
};

struct gene_pred_line {
    // mark the memory is allocated or not inside this struct, clear == 0 is empty, clear == 1 recall clear_gene_pred_line to free it.
    int clear;    
    char *chrom;
    uint32_t txstart;
    uint32_t txend;
    char strand;
    char *name1;
    char *name2; // transcript name or null for some gene_pred file
    uint32_t cdsstart;
    uint32_t cdsend;
    uint32_t exoncount;
    struct exon_pair *exons;
    struct exon_offset_pair *dna_ref_offsets;
};
enum func_region_type {
    func_region_unknown,
    func_region_split_sites,
    func_region_cds,
    func_region_intron,
    func_region_utr5,
    func_region_utr3,
    func_region_intergenic,
};
//  hgvs_core keeps transcript/protein name, the format of hgvs name is construct by
//                                    l_name   l_type
//                                    |        |
//  [transcript|protein|gene|ensemble]:"prefix".postion...
//
//  l_name          l_type
//  |               |
//  [chrom]:"prefix".postion...
struct hgvs_names_core {    
    uint16_t l_name; // name offset in the data cahce, for chrom l_name == 0;
    uint16_t l_type; // the byte after type offset, usually offset of '.'
    kstring_t str;
    enum func_region_type type;
};

// Dynamic allocate and init technology.
// For struct { int l, m, i; void *a }, a is the cached array of predefined type. l for used length, m for max length,
// i for inited length. And i always >= l, and m always >= i. If l == m, the cache array should be reallocated by a new
// memory size. This structure used to reuse the cache complex structures other than points.
struct hgvs_names {
    int l, m, i; 
    struct hgvs_names_core *a;
};
struct hgvs_names_cache {
    int l, m, i; // l == n_allele -1
    struct hgvs_names *a;
};
struct name_list {
    int l, m;
    char **a;
};

struct gene_pred_memory_pool {
    int rid;
    uint32_t start;
    uint32_t end;
    int l, m, i; // l for used length, m for max length, i for inited length
    struct gene_pred_line *a;
    struct hgvs_names_cache cache;
};
void gene_pred_line_clear(struct gene_pred_line *l)
{
    if (l->clear == 1) return;
    if (l->clear == 0) {
	free(l->chrom);
	free(l->name1);
	free(l->name2);
	free(l->exons);
	free(l->dna_ref_offsets);
    }
    l->chrom = NULL;
    l->name1 = NULL;
    l->name2 = NULL;
    l->exons = NULL;
    l->dna_ref_offsets = NULL;
    l->txstart = l->txend = 0;
    l->cdsstart = l->cdsend = 0;
    l->exoncount = 0;
    l->clear = 1;
}
void hgvs_names_core_clear(struct hgvs_names_core *c)
{
    c->str.l = 0;
    c->l_name = 0;
    c->l_type = 0;
    c->type = func_region_unknown;
}
void hgvs_names_clear(struct hgvs_names *hn)
{
    int i;
    for (i = 0; i < hn->l; ++i) {
	hgvs_names_core_clear(&hn->a[i]);				   
    }
    hn->l = 0;
}
void hgvs_names_cache_clear(struct hgvs_names_cache *cache)
{
    int i;
    for (i = 0; i < cache->l; ++i)
	hgvs_names_clear(&cache->a[i]);
    cache->l = 0;
}
void hgvs_names_cache_destroy(struct hgvs_names_cache *cache)
{
    int i, j;
    for (i = 0; i < cache->i; ++i) {
	struct hgvs_names *hn = &cache->a[i];
	for (j = 0; j < hn->i; ++j) {
	    if (hn->a[j].str.m) free(hn->a[j].str.s);	    
	}
	free(hn->a);	    
    }
    free(cache->a);
}
glist_t *init_gene_list(const char *list)
{
    int i, n = 0;
    char **names = hts_readlist(list, 1, &n);
    if (n== 0) return NULL;
    glist_t *glist = kh_init(list);
    int ig, k;
    for (i=0; i<n; ++n) {
       k = kh_put(list, glist, names[i], &ig);
    }
    return glist;
}

tbx_t *load_gene_pred_file(const char *fn)
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
struct format_type {
    int chrom_col;
    int name1_col;
    int name2_col; // name2 could be empty;
    int strand_col;
    int txstart_col;
    int txend_col;
    int cdsstart_col;
    int cdsend_col;
    int exoncount_col;
    int exonstarts_col;
    int exonends_col;
};
const struct format_type refgene_format_cols = {
    .chrom_col = 2,
    .name1_col = 1,
    .name2_col = 12,
    .strand_col = 3,
    .txstart_col = 4,
    .txend_col = 5,
    .cdsstart_col = 6,
    .cdsend_col = 7,
    .exoncount_col = 8,
    .exonstarts_col = 9,
    .exonends_col = 10,
};
const struct format_type genepred_format_cols = {
    .chrom_col = 1,
    .name1_col = 0,
    .name2_col = 11,
    .strand_col = 2,
    .txstart_col = 3,
    .txend_col = 4,
    .cdsstart_col = 5,
    .cdsend_col = 6,
    .exoncount_col = 7,
    .exonstarts_col = 8,
    .exonends_col = 9,
};
struct anno_data_file_handler {
    struct format_type type;
    const char *fn;
    htsFile *fp;
    tbx_t *tbx;
};
// test the format of anno dataset, usually genePred or refgene
void test_annodata_file_type(struct anno_data_file_handler *handler)
{
    
}
void gene_pred_line_praser(kstring_t *str, struct gene_pred_line *line)
{
    int nfields = 0;
    int *splits = ksplit(str, 0, &nfields);
    if (line->clear != 1) gene_pred_line_clear(line);
    line->clear = 0;

    // accept any gene_pred-like format, like refGene, ensGene, gene_pred
    // assert(nfields == 16);

    char *s = str->s;
    char *chrom = s + splits[refgene_format_cols.chrom_col];
    char *name1 = s + splits[refgene_format_cols.name1_col];
    char *strand = s + splits[refgene_format_cols.strand_col];
    char *txstart = s + splits[refgene_format_cols.txstart_col];
    char *txend = s + splits[refgene_format_cols.txend_col];
    char *cdsstart = s + splits[refgene_format_cols.cdsstart_col];
    char *cdsend = s + splits[refgene_format_cols.cdsend_col];
    char *exoncount = s + splits[refgene_format_cols.exoncount_col];   
    char *name2 = NULL;

    line->chrom = strdup(chrom);
    line->name1 = strdup(name1);
    line->strand = memcmp(strand, "+", 1) ? '-' : '+';
    line->txstart = atoi(txstart);
    line->txend = atoi(txend);
    line->cdsstart = atoi(cdsstart);
    line->cdsend = atoi(cdsend);
    line->exoncount = atoi(exoncount);
    // for some genePred file, no name2 specified.
    line->name2 = NULL;
    if (refgene_format_cols.name2_col > 0) {
	name2 = s + splits[refgene_format_cols.name2_col];
	line->name2 = strdup(s + splits[12]);
    }
    // the exons region in gene_pred file is like  1,2,3,4,5. try to init the exon_pair[] 
    // and exon_offset_pair[] by exoncount 
    char *ss = s + splits[9];
    char *se = s + splits[10];
    line->exons = (struct exon_pair*)calloc(line->exoncount, sizeof(struct exon_pair));
    line->dna_ref_offsets = (struct exon_offset_pair*)calloc(line->exoncount, sizeof(struct exon_offset_pair));
    free(splits);
    int i;
    char *ss1, *se1;

    // calculate the length of function regions 
    int dna_reference_length = 0;
    int dna_read_reference_length = 0;
    int dna_forward_reference_length = 0;
    int dna_backward_reference_length = 0;
    int is_coding_reference = line->cdsstart == line->cdsend ? 0 : 1;
    
    for (i = 0; i < line->exoncount; ++i) {
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
	
	// add exon length to dna reference length
	dna_reference_length += line->exons[i].end - line->exons[i].start;

	if (is_coding_reference == 0) continue;

	if (line->exons[i].end < line->cdsstart) {
	    dna_forward_reference_length += line->exons[i].end - line->exons[i].start;
	} else {
	    // first cds 
	    if (line->cdsstart > line->exons[i].start) {
		dna_forward_reference_length += line->cdsstart - line->exons[i].start;
	    } else {
		// come to end regions 
		if (line->cdsend <= line->exons[i].start) {
		    dna_backward_reference_length += line->exons[i].end - line->exons[i].start;
		} else {
		    if (line->cdsend < line->exons[i].end) {
			dna_backward_reference_length += line->exons[i].end - line->cdsend;
		    }		   
		}
	    }
	}
    } // end for

    // dna_read_reference_length is the coding reference length for coding transcript or length of noncoding transcript 
    dna_read_reference_length = dna_reference_length - dna_forward_reference_length - dna_backward_reference_length;

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
    int dna_read_reference_forward_offset = 0;
    int dna_read_reference_backward_offset = 0;
    for (; l1<line->exoncount && dna_forward_reference_length; ++l1) {
	int32_t length_exon = line->exons[l1].end - line->exons[l1].start;
	line->dna_ref_offsets[l1].start =
	    (dna_forward_reference_length << OFFSET_BITWISE) |
	    (line->strand == '+' ? REG_UTR5 : REG_UTR3);

	if (dna_forward_reference_length > length_exon) {
	    dna_forward_reference_length -= length_exon;		
	    line->dna_ref_offsets[l1].end =
		((dna_forward_reference_length+1) << OFFSET_BITWISE) |
		(line->strand == '+' ? REG_UTR5 : REG_UTR3);
	} else {
	    if (line->strand == '+') {
		dna_read_reference_forward_offset = length_exon - dna_forward_reference_length;
	    } else {
		dna_read_reference_forward_offset = dna_read_reference_length + dna_forward_reference_length - length_exon +1; // for minus strand count from backward
	    }
	    line->dna_ref_offsets[l1].end = (dna_read_reference_forward_offset<<OFFSET_BITWISE)| REG_CODING;
	    dna_forward_reference_length = 0;
	    ++l1;
	    break;
	}
    }

    // count backward
    for (; l2 > 0 && dna_backward_reference_length; --l2) {

	int32_t length_exon = line->exons[l2].end - line->exons[l2].start;

	line->dna_ref_offsets[l2].end =
	    (dna_backward_reference_length << OFFSET_BITWISE) |
	    (line->strand == '+' ? REG_UTR3 : REG_UTR5);

	if (dna_backward_reference_length > length_exon) {
	    dna_backward_reference_length -= length_exon;
	    line->dna_ref_offsets[l2].start =
		((dna_backward_reference_length+1) << OFFSET_BITWISE) |
		(line->strand == '+' ? REG_UTR3 : REG_UTR5);
	} else {
	    if (line->strand == '+') {
		dna_read_reference_backward_offset = dna_read_reference_length + dna_backward_reference_length - length_exon + 1;
	    } else {
		dna_read_reference_backward_offset = length_exon - dna_backward_reference_length;
	    }
	    line->dna_ref_offsets[l2].start = (dna_read_reference_backward_offset<<OFFSET_BITWISE) | REG_CODING;
	    dna_backward_reference_length = 0;
	    --l2;
	    break;
	}
    }

    // count inter regions
    if (line->strand == '+') {
	if (is_coding_reference)
	    dna_read_reference_length = dna_read_reference_backward_offset -1;       
	int l;
	for (l=l2; l>=l1;  l--) {
	    int32_t length_exon = line->exons[l].end - line->exons[l].start;
	    line->dna_ref_offsets[l].end =
		(dna_read_reference_length << OFFSET_BITWISE) |
		(is_coding_reference ? REG_CODING : REG_NONCODING);

	    dna_read_reference_length -= length_exon;

	    line->dna_ref_offsets[l].start = ((dna_read_reference_length+1)<< OFFSET_BITWISE) |
		(is_coding_reference ? REG_CODING : REG_NONCODING);
	}
    } else {
	if (is_coding_reference) {
	    dna_read_reference_length = dna_read_reference_forward_offset -1;
	}
	int l;
	for (l=l1; l<=l2; l++) {
	    int32_t length_exon = line->exons[l].end - line->exons[l].start;
	    
	    line->dna_ref_offsets[l].start = (dna_read_reference_length << OFFSET_BITWISE) |
		(is_coding_reference ? REG_CODING : REG_NONCODING);

	    dna_read_reference_length -= length_exon;
	    line->dna_ref_offsets[l].end = ((dna_read_reference_length+1) << OFFSET_BITWISE) |
		(is_coding_reference ? REG_CODING : REG_NONCODING);
	}
    }
}
void generate_dbref_database(struct gene_pred_line *line)
{
    if (line->clear == 1) return;
    int i, j;
    for (i=0; i<line->exoncount; ++i) {
    	kstring_t temp[2] = {KSTRING_INIT, KSTRING_INIT};
    	int types[2];
    	types[0] = line->dna_ref_offsets[i].start & REG_MASK;
    	types[1] = line->dna_ref_offsets[i].end & REG_MASK;
    	int j;
	// [start, end] 
    	for (j=0; j<2; ++j) {
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
void push_mempool(struct gene_pred_memory_pool *pool, kstring_t *str)
{
    if (pool->m == pool->l) {
	pool->m = pool->m == 0 ? 2 : pool->m << 1;
	pool->a = (struct gene_pred_line *)realloc(pool->a, sizeof(struct gene_pred_line)*pool->m);
    }
    // i should always greater than l. if i == l increase i to init a new line for future use. go abort if i < l
    if (pool->i == pool->l) {
	pool->a[pool->i++].clear = -1;
    }
    // the prase func must return a point or go abort
    gene_pred_line_praser(str, &pool->a[pool->l]);
    pool->l++;
}
void update_mempool(struct gene_pred_memory_pool *pool)
{    
    if (pool->l == 0) return;
    // try to loop this pool and find the edges
    pool->start = pool->a[0].txstart;
    pool->end = pool->a[0].txend;
    int i;    
    for (i=1; i<pool->l; ++i) {
	// if (pool->a[i] == NULL) continue;
	if (pool->a[i].txstart < pool->start) pool->start = pool->a[i].txstart;
	if (pool->a[i].txend > pool->end) pool->end = pool->a[i].txend;
    }
}
void fill_mempool(struct gene_pred_memory_pool *pool, htsFile *fp, tbx_t *tbx, int rid, uint32_t start, uint32_t end)
{
    hts_itr_t *itr = tbx_itr_queryi(tbx, rid, start, end);
    kstring_t str = KSTRING_INIT;
    if (rid == -1) return;
    
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
	push_mempool(pool, &str);
	str.l = 0;
    }
    pool->rid = rid;
    if (str.m) free(str.s);
    tbx_itr_destroy(itr);
    update_mempool(pool);    
}
void memory_pool_clear(struct gene_pred_memory_pool *pool)
{
    int i;
    for (i=0; i<pool->l; ++i) {
	gene_pred_line_clear(&pool->a[i]);
    }
    pool->l = 0;
    // i and m should be kept for future use
    pool->start = 0;
    pool->end = 0;
    pool->rid = -1;
    hgvs_names_cache_destroy(&pool->cache);
}

// HGVS nomenclature : 
// DNA recommandations *
// - substitution variant, 
// Format: “prefix”“position_substituted”“reference_nucleoride””>”new_nucleoride”, e.g. g.123A>G
// - deletion variant, 
// Format: “prefix”“position(s)_deleted”“del”, e.g. g.123_127del
// - duplication variant,
// Format: “prefix”“position(s)_duplicated”“dup”, e.g. g.123_345dup
// - insertion variant,
// Format: “prefix”“positions_flanking”“ins”“inserted_sequence”, e.g. g.123_124insAGC
// - inversion variant,
// Format: “prefix”“positions_inverted”“inv”, e.g. g.123_345inv
// - conversion variant,
// Format: “prefix”“positions_converted”“con”“positions_replacing_sequence”, e.g. g.123_345con888_1110
// - deletion-insertion variant,
// Format: “prefix”“position(s)_deleted”“delins”“inserted_sequence”, e.g. g.123_127delinsAG
// - repeated sequences variant,
// Format (repeat position): “prefix”“position_repeat_unit””["”copy_number””]”, e.g. g.123_125[36]
// 
// [ reference : http:// www.HGVS.org/varnomen ]
enum hgvs_variant_type {
    var_type_ref = 0,
    var_type_snp,
    var_type_dels,
    var_type_ins,
    var_type_delins,
    var_type_copy, // dup for ins
    var_type_complex, 
    var_type_unknow,
};

struct hgvs_names_description {
    // if type ==  var_type_ref, hgvs name generater will skip to construct a name, only type will be inited, 
    // DONOT use any other value in this struct then 
    enum hgvs_variant_type type;
    uint8_t strand; // for plus strand start <= end, for minus start >= end
    int32_t start; // position on ref
    // for var_type_snp end == start, for var_type_dels end > start, for var_type_ins end = start +1, 
    // for delvar_type_ins end > start 
    int32_t end;
    // for var_type_snp ref_length == 1, for var_type_dels ref_length == length(dels),
    // for var_type_ins ref_length == 0, for var_type_delins ref_length == length(var_type_dels) 
    int ref_length;
    char *ref;
    // for var_type_snp alt_length == 1, for var_type_dels alt_length == 0, for var_type_ins alt_length ==
    // length(var_type_ins), for delvar_type_ins alt_length == length(var_type_ins) 
    int alt_length;
    char *alt;
    // the copy number only valid if variants type is var_type_copy, for most case, my algrithm only check
    // mark var_type_dels or var_type_ins in the first step, and check the near sequences from refseq 
    // databases and if there is copy number changed (more or less), the variants will update to var_type_copy 
    int copy_number;
    int copy_number_ori;
};
void hgvs_names_description_destory(struct hgvs_names_description *des)
{
    if (des->ref_length) free(des->ref);
    if (des->alt_length) free(des->alt);
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
    for (i = 0; i < n; ++i) {
	//debug_print("i : %d, n : %d", i, n);
	//debug_print("rev : %c, seq : %c,",rev[i], dna_seqs[n-i-1]);
	rev[i] = rev_seqs_matrix[dna_seqs[n-i-1]];
	//debug_print("i : %d, rev : %c, seq : %c, n : %d", i, rev[i], dna_seqs[n-i-1], n);
    }
    rev[n] = '\0';
    return rev;
}
void hgvs_names_description_reverse(struct hgvs_names_description *des, uint8_t strand)
{
    if (des->type == var_type_ref) return;
    if (des->strand == strand) return;
    if (des->strand == 0 && strand == '+') return;
    des->strand = strand;
    char *ref = des->ref;
    char *alt = des->alt;
    int32_t start = des->end;
    des->end = des->start;
    des->start = start;
    //debug_print("ref : %s, %d, alt : %s, %d", des->ref, des->ref_length, des->alt, des->alt_length);
    des->ref = rev_seqs(des->ref, des->ref_length);
    des->alt = rev_seqs(des->alt, des->alt_length);
    if (des->ref_length) free(ref);
    if (des->alt_length) free(alt);
}
struct hgvs_names_description *describe_variants(const char *ref, const char *alt, int _pos)
{
    struct hgvs_names_description *des = (struct hgvs_names_description*)malloc(sizeof(struct hgvs_names_description));
    const char *a = alt;
    const char *r = ref;
    int pos = _pos; // pos may differ from _pos, but _pos will not change in this strnduption
    while (*a && *r && toupper(*a) == toupper(*r)) { a++; r++; pos++; }
    if (!a[0] && !r[0]) {
	des->type = var_type_ref;
	return des;
    }
    des->strand = 0; // 0 for unknown, '+' for plus, '-' for minus
    if (*a && *r && !a[1] && !r[1] ) {
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

    if (*a && !*r) {
	des->type = var_type_ins;
	while ( *a ) a++;
	des->start = pos;
	des->ref_length = 0;
	des->ref = str_init;
	des->end = pos +1;
	des->alt_length = (a-alt) - (r-ref);
	des->alt = strndup(a-des->alt_length, des->alt_length);
	return des;
    } else if (!*a && *r) {
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

    const char *ae = a;
    const char *re = r;
    while ( ae[1] ) ae++;
    while ( re[1] ) re++;
    while ( re > r && ae > a && toupper(*re) == toupper(*ae)) {
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
    } else if (re == r) {
	des->type = var_type_ins;
	des->start = pos;
	des->end = pos+1;
	des->ref_length = 0;
	des->ref = str_init;
	des->alt_length = ae - a;
	des->alt = strndup(a, ae-a);
	return des;
    }

    des->type = var_type_delins;
    des->start = pos;
    des->ref_length = re-r;
    des->ref = strndup(r, re-r);
    des->alt_length = ae-a;
    des->alt = strndup(a, ae-a);
    des->end = pos + des->ref_length -1;
    return des;
}
void check_cnv_hgvs_names_descriptions(struct hgvs_names_description *des, struct gene_pred_line *line, faidx_t *fai, htsFile *fp)
{
    if (des->type == var_type_ins) {

    } else if (des->type == var_type_dels) {

    }
}
void find_exons_loc(uint32_t pos, int exoncount, struct exon_pair *pair, int *l1, int *l2 )
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
    while(*l1 < *l2) {
	if ( *l2 - *l1 == 1 ) break;
	/* uint32_t start = (*l1) & 1 ? pair[(*l1)/2].end : pair[(*l1)/2].start; */
	/* uint32_t end =  (*l2) & 1 ? pair[(*l2)/2].end : pair[(*l2)/2].start; */
	/* debug_print("start : %u, end : %u, pos : %u", start, end, pos); */
	/* assert(pos >= start && pos <= end); */
	int l = *l1 + *l2;
	l = l & 1 ? l/2 + 1 : l/2;
	uint32_t iter = l & 1 ?  pair[l/2].end : pair[l/2].start;
	if ( iter > pos ) *l2 = l;
	else *l1 = l;
    } 
}
enum func_region_type pos_convert(int32_t pos, int exoncount, int strand, struct exon_pair *pair, struct exon_offset_pair *locs, int *is_coding, char **cpos)
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
    if (l1 & 1) {	
	// the nestest edge is start of exon l2/2
	uint32_t loc;
	uint8_t type;
	int offset;
	// find the most nearest edge
	if (offset_start > offset_end) { // cap to end
	    loc = loc_end;
	    type = type_end;
	    offset = strand == '+' ? -offset_end : offset_end;
	} else {
	    loc = loc_start;
	    type = type_start;
	    offset = strand == '+' ? offset_start : -offset_start;
	}

	if (type & REG_NONCODING) {
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
    if (type_start == type_end) {
	loc = loc_end > loc_start ?  loc_end - offset_end : loc_start - offset_start;
	type = type_start;
    } else {

	if (type_start == REG_UTR5) { // should only be plus strand
	    loc = loc_start - offset_start;
	    type = type_start;
	    if (type_end == REG_CODING) {
		if (loc <= 0) {
		    loc = -loc +1;
		    type = type_end;
		}
	    } else if(type_end == REG_UTR3) { // the cds region inside one exon, closed with UTRs
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
	    } else {
		// impossible
		error("This is a impossible stituation. type_end : %d", type_end);
	    }		
	} else if (type_start == REG_CODING) {

	    if (type_end == REG_UTR3) { // strand plus
		loc = loc_end - offset_end;
		type = type_end;
		if (loc <= 0) {
		    loc = loc_start + offset_start;
		    type = type_start;
		}
	    } else if (type_end == REG_UTR5) { // strand minus
		loc = loc_start - offset_start;
		type = type_start;
		if (loc <= 0) {
		    loc = -loc + 1;
		    loc = type_end;
		}
	    } else {
		error("This is a impossible stituation. type_end : %d", type_end);
	    }
	} else if (type_start == REG_UTR3) { // should only be minus strand
	    if (type_end == REG_CODING) {
		loc = loc_start - offset_start;
		type = type_start;
		if (loc <= 0) {
		    loc = -loc + 1;
		    type = type_end;
		}		    
	    } else if (type_end == REG_UTR5) { // cds region inside one exon
		loc = loc_start - offset_start;
		type = type_start;
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
	    } else {
		error("This is a impossible stituation. type_end : %d", type_end);
	    }
	} else {
	    error("This is a impossible stituation. type_start : %d", type_start);
	}
    }
    assert(loc>0);
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
void generate_hgvs_names_core(struct gene_pred_line *gp, struct hgvs_names_description *des, struct hgvs_names_core *c, int *valid)
{
    hgvs_names_core_clear(c);
    *valid = 0;
    if (des->start > gp->txend || des->end < gp->txstart) return;    
    uint32_t start_var = (uint32_t)des->start;
    uint32_t end_var = (uint32_t)des->end;
    kstring_t *str = &c->str;
    // if no transcript name, use gene name then
    char *name = gp->name1 == NULL || gp->name1[0] == '.' ? gp->name2 : gp->name1;
    kputs(name, str);
    c->l_name = str->l;    
    assert(c->l_name);
    kputc(':', str);
    char *pos1;
    int is_coding = -1;

    c->type = pos_convert(des->start, gp->exoncount, gp->strand, gp->exons, gp->dna_ref_offsets, &is_coding, &pos1);
    if (is_coding == 1) kputs("c.", str);
    else if (is_coding == 0) kputs("n.", str);
    else error("Unknown transcript type! %s", name);
    kputs(pos1, str);

    c->l_type = str->l -1;
    copy_seqs_func func = gp->strand == '+' ? strndup : rev_seqs;
    char *ref = des->ref_length == 0 ? NULL : func(des->ref, des->ref_length);
    char *alt = des->alt_length == 0 ? NULL : func(des->alt, des->alt_length);
    int is_coding1;
    if (des->type == var_type_snp) {
	ksprintf(str, "%s>%s", ref, alt);
    } else if (des->type == var_type_dels) {
	if (des->ref_length == 1) {
	    ksprintf(str, "del%s", des->ref);
	} else {
	    kputc('_', str);
	    char *pos2;
	    pos_convert(des->end, gp->exoncount, gp->strand, gp->exons, gp->dna_ref_offsets, &is_coding1, &pos2);
	    kputs(pos2, str);
	    free(pos2);
	    kputs("del", str);
	}
    } else if (des->type == var_type_ins) {
	kputc('_', str);
	char *pos2;
	pos_convert(des->end, gp->exoncount, gp->strand, gp->exons, gp->dna_ref_offsets, &is_coding1, &pos2);
	kputs(pos2, str);
	free(pos2);	
	if ( des->alt_length < 20) ksprintf(str, "ins%s", alt);
	else ksprintf(str, "ins%d", des->alt_length);
    } else if (des->type == var_type_delins) {
	if (des->ref_length > 1) {
	    kputc('_', str);
	    char *pos2;
	    pos_convert(des->end, gp->exoncount, gp->strand, gp->exons, gp->dna_ref_offsets, &is_coding1, &pos2);
	    kputs(pos2, str);
	    free(pos2);		    
	}
	if ( des->alt_length < 20) ksprintf(str, "delins%s", alt);
	else ksprintf(str, "delins%d", des->alt_length);
    } else {
	error("unknow type : %d", des->type);
    }
#ifdef DEBUG_MODE
    debug_print("%d : %s", des->start,str->s);
#endif
    free(pos1);
    // it is valid
    *valid = 1;
    
    if (des->ref_length != 0) free(ref);
    if (des->alt_length != 0) free(alt);    
}
void name_list_push(struct name_list *names, char *name)
{
    int i;
    for (i = 0; i < names->l; ++i) {
	if (strcmp(names->a[i], name) == 0)
	    return;
    }
    if (i == names->l) {
	if (names->m == names->l) {
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
// @names1      list of gene names
// @n_ales      related allele numbers (R)
// @n_hgvs      number of hgvs names strings, several transcript for one allele 
// should be cache into one string like NM_xxxx:c.xx; NC_xxxx:n.xxx; ...; 
// @hgvs_names  list of hgvs names strings
// @n_names2    number of transcripts
// @names2      list of transcripts
void generate_hgvs_names(struct gene_pred_memory_pool *pool, const bcf1_t *line, int *n_names2, char **names2, int *n_ale)
{
    *n_names2 = 0;
    *n_ale = line->n_allele -1;
    if (*n_ale == 0) return;
    //bcf_dec_t *d = &line->d;
    struct hgvs_names_cache * cache = &pool->cache;
    // a simple list structure inited for caching gene names and de-duplicate 
    struct name_list names = KSTRING_INIT;
    if (cache->m < line->n_allele) {
	cache->m = line->n_allele;
	cache->a = (struct hgvs_names*)realloc(cache->a, sizeof(struct hgvs_names)*cache->m);
    }
    cache->l = 0;
    // roadmap : snv -> del -> ins -> dup ->delins
    int i, j;    
    for (i = 1; i < line->n_allele; ++i) {
	if (cache->l == cache->i) {
	    cache->a[cache->i].l = cache->a[cache->i].m = cache->a[cache->i].i = 0;
	    cache->a[cache->i].a = NULL;
	    cache->i++;
	}
	struct hgvs_names *hn = &cache->a[cache->l++];
	hgvs_names_clear(hn);
	struct hgvs_names_description *des = describe_variants(line->d.allele[0], line->d.allele[i], line->pos +1);
	for (j = 0; j < pool->l; ++j) {
	    int valid = 0;
	    if (hn->l == hn->m) {
		hn->m = hn->m == 0 ? 2 : hn->m << 1;
		hn->a = (struct hgvs_names_core*)realloc(hn->a, sizeof(struct hgvs_names_core)*hn->m);
	    }
	    if (hn->i == hn->l) {
		kstring_t *str = &hn->a[hn->i].str;
		str->l = str->m = 0;
		str->s = NULL;
		hn->i++;
	    }
	    struct hgvs_names_core *core= &hn->a[hn->l];
	    struct gene_pred_line *gp = &pool->a[j];
	    hgvs_names_description_reverse(des, gp->strand);
	    generate_hgvs_names_core(gp, des, core, &valid);
	    if (valid == 0) continue;
	    name_list_push(&names, gp->name2);
	    hn->l++;
	}
	hgvs_names_description_destory(des);
	
    }

    kstring_t str = KSTRING_INIT;
    int n = 0;
    for (i = 0; i < names.l; i++) {
	if ( i ) kputc(',', &str);
	kputs(names.a[i], &str);
	free(names.a[i]);
    }
    if (names.m) free(names.a);
    *n_names2 = names.l;
    *names2 = str.l ? str.s : NULL;
}
// number : R, type : string
static void setter1_hgvs_names(bcf_hdr_t *hdr, bcf1_t *line, char *string) 
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
	setter1_hgvs_names(hdr, line, string);		
    } else if (strcmp(key, "Transcript") == 0) {
	setter1_transcripts_names(hdr, line, string);
    } else {
	error("Unrecongnized tag %s, only Gene, HGVSDNA, and Transcript are supported for now.", key);
    }    
}
char *generate_transcript_string(struct hgvs_names_cache *cache, int n)
{
    kstring_t str = KSTRING_INIT;
    int i, j;    
    for (i=0; i<n; ++i) {
	if (i) kputc(',', &str);
	// foreach allele
	struct hgvs_names *hn = &cache->a[i];
	for ( j = 0; j < hn->l; ++j ) {
	    struct hgvs_names_core *core = &hn->a[j];
	    if (core->l_name) kputsn(core->str.s, core->l_name-1, &str);
	    else kputc('.', &str);
	    kputc('|', &str);
	}
    }
    return str.s;
}
char *generate_hgvsvarnomen_string(struct hgvs_names_cache *cache,int n)
{
    kstring_t str = KSTRING_INIT;
    int i, j;    
    for (i=0; i<n; ++i) {
	if (i) kputc(',', &str);
	// foreach allele
	struct hgvs_names *hn = &cache->a[i];
	for ( j = 0; j < hn->l; ++j ) {
	    struct hgvs_names_core *core = &hn->a[j];
	    if (core->l_name) kputs(core->str.s, &str);
	    else kputc('.', &str);
	    kputc('|', &str);
	}
    }
    return str.s;
}

// anno_hgvs_core() only used to annotate bcf/vcf standalone, all the valid tags-added functions will be called in 
// this function, to annotate the vcf more specification, use setter_genepred_* functions instead of it. 
void anno_hgvs_core(struct gene_pred_memory_pool *pool, htsFile *fp, tbx_t *tbx, bcf_hdr_t *hdr, bcf1_t *line)
{
    // retrieve the regions this variants located first, check the memory pool and update the pool if the regions is out of
    // cached positions. just skip if the variant type of line is a ref. 
    if (line->pos <= 0) return;
    if (bcf_get_variant_types(line) == VCF_REF) return;

    uint32_t end = bcf_calend(line);
    int id = tbx_name2id(tbx, bcf_hdr_id2name(hdr, line->rid));
    if (pool->rid == -1 || pool->rid != line->rid || pool->start > line->pos+1 || pool->end < line->pos+1) {
	// fill mempool, skip if no record in the memory pool
	pool->l = 0;
	fill_mempool(pool, fp, tbx, id, line->pos, end);
	if (pool->l == 0) return;
    }
    int n_names2 = 0;
    char *names2 = NULL;
    int n_ale = 0;
    
    generate_hgvs_names(pool, line, &n_names2, &names2, &n_ale);
    if (n_ale == 0) return;
    char *trans_string = generate_transcript_string(&pool->cache, n_ale);
    char *hgvs_string = generate_hgvsvarnomen_string(&pool->cache, n_ale);
    //debug_print("trans : %s", trans_string);
    //debug_print("hgvs : %s", hgvs_string);

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

#include <sys/time.h>
static double get_time() {
    struct timeval tv;
    if ( gettimeofday(&tv, 0) != 0 ) 
	error("Failed to get time of day.");
    return (double)tv.tv_sec + 1.0e-6 *(double)tv.tv_usec;
}
struct args {
    const char *gene_pred_file;
    const char *refseq_file;
    const char *out_file;
    const char *out_type;
    const char *input_fname;
    htsFile *fp;
    htsFile *fout;
    tbx_t *gpidx;
    faidx_t *faidx;
    glist_t *glist;
    struct gene_pred_memory_pool pool;
} args = {
    .gene_pred_file = 0,
    .refseq_file = 0,
    .out_file = 0,
    .out_type = 0,
    .input_fname = 0,
    .gpidx = 0,
    .faidx = 0,
    .glist = 0,
};

void clear_args(struct args *args)
{
    tbx_destroy(args->gpidx);
    if (args->faidx) fai_destroy(args->faidx);
    if (args->glist) kh_destroy(list, args->glist);
    memory_pool_clear(&args->pool);
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
	else if (strcmp(a, "-data") == 0 && args.gene_pred_file == 0)
	    arg_var = &args.gene_pred_file;
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
 
    if (args.gene_pred_file == 0)
	error("Reasons :\n"
	      "-data gene preditions database is needed.\n"
	      "You could download refGene.txt.gz or ensGene.txt.gz from UCSC websites and sort and indexed by tabix.");

    // if (args.refseq_file == 0) 
    // 	error("Reasons :\n" 
    // 	      "-refseq refseq.fa.gz file is need.\n" 
    // 	      "This is the transcripts reference sequences in fasta format." 
    // 	      "Notice the transcripts names should be consistance with gene preditions databases."); 

    htsFile *fgp = hts_open(args.gene_pred_file, "r");
    if (fgp == NULL)
	error("Failed to open %s", args.gene_pred_file);
    args.gpidx = load_gene_pred_file(args.gene_pred_file);
    if (args.gpidx == NULL)
	error("Failed to load index of %s", args.gene_pred_file);

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
    struct gene_pred_memory_pool *pool = &args.pool;
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
    int id = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "Gene");
    if (id == -1) {
	bcf_hdr_append(hdr_out, "##INFO=<ID=Gene,Number=1,Type=String,Description=\"Gene names\">");
	bcf_hdr_sync(hdr_out);
	id = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "Gene");
	assert(bcf_hdr_idinfo_exists(hdr_out, BCF_HL_INFO, id));
    }
    id = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "Transcript");
    if (id == -1) {
	bcf_hdr_append(hdr_out, "##INFO=<ID=Transcript,Number=G,Type=String,Description=\"Transcript names\">");
	bcf_hdr_sync(hdr_out);
	id = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "Transcript");
	assert(bcf_hdr_idinfo_exists(hdr_out, BCF_HL_INFO, id));
    }
    id = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "HGVSDNA");
    if (id == -1) {
	bcf_hdr_append(hdr_out, "##INFO=<ID=HGVSDNA,Number=G,Type=String,Description=\"HGVS nomenclature for the description of DNA sequence variants\">");
	bcf_hdr_sync(hdr_out);
	id = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "HGVSDNA");
	assert(bcf_hdr_idinfo_exists(hdr_out, BCF_HL_INFO, id));
    } 
#ifndef DEBUG_MODE   
    bcf_hdr_write(args.fout, hdr_out);
#endif
    // init gene_pred or refgene database, hold tabix index cache and faidx cache in memory 
    bcf1_t *line = bcf_init();
    while ( bcf_read(args.fp, hdr, line) == 0 ) {     
	anno_hgvs_core(pool, fgp, args.gpidx, hdr_out, line);
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
