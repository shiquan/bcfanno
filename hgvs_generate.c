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

KHASH_MAP_INIT_STR(list, char*)
typedef khash_t(list) glist_t;

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
 
#define DNA_REF_OFFSETS_BITS   4
#define FUNC_REGION_NONCODING  1
#define FUNC_REGION_CODING     2
#define FUNC_REGION_UTR5       4
#define FUNC_REGION_UTR3       8
#define FUNC_REGION_MASK       0xF

struct exon_offset_pair {
    uint32_t start;
    uint32_t end;
};

struct gene_predictions_line {
    // mark the memory is allocated or not inside this struct, clear == 0 is empty, clear == 1 recall clear_gene_predictions_line to free it.
    int clear;
    
    char *chrom;
    uint32_t txstart;
    uint32_t txend;
    char strand;
    char *name1;
    char *name2; // transcript name or null for some gene_predictions file
    uint32_t cdsstart;
    uint32_t cdsend;
    uint32_t exoncount;
    struct exon_pair *exons;
    struct exon_offset_pair *dna_ref_offsets;
};

struct gene_predictions_memory_pool {
    int tid;
    uint32_t start;
    uint32_t end;
    int l, m;
    struct gene_predictions_line *a;
};

void clear_gene_predictions_line(struct gene_predictions_line *);

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

tbx_t *load_gene_predictions_file(const char *fn)
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

void gene_predictions_prase_line(kstring_t *str, struct gene_predictions_line *line)
{
    int nfields = 0;
    int *splits = ksplit(str, 0, &nfields);
    if (line->clear == 1) clear_gene_predictions_line(line);
    // accept any gene_predictions-like format, like refGene, ensGene, gene_predictions
    assert(nfields == 16);
    line->chrom = strdup(str->s + splits[2]);
    line->name1 = strdup(str->s + splits[1]);
    line->strand = memcmp(str->s + splits[3], "+", 1) ? '-' : '+';
    line->txstart = atoi(str->s + splits[4]);
    line->txend = atoi(str->s + splits[5]);
    line->cdsstart = atoi(str->s + splits[6]);
    line->cdsend = atoi(str->s + splits[7]);
    line->exoncount = atoi(str->s + splits[8]);
    line->name2 = strdup(str->s + splits[12]);
    // the exons region in gene_predictions file is like  1,2,3,4,5. try to init the exon_pair[] 
    // and exon_offset_pair[] by exoncount 
    char *ss = str->s + splits[9];
    char *se = str->s + splits[10];
    line->exons = (struct exon_pair*)calloc(line->exoncount, sizeof(struct exon_pair));
    line->dna_ref_offsets = (struct exon_offset_pair*)calloc(line->exoncount, sizeof(struct exon_offset_pair));

    int i;
    char *ss1, *se1;

    // calculate the length of function regions 
    int dna_reference_length = 0;
    int dna_read_reference_length = 0;
    int dna_forward_reference_length = 0;
    int dna_backward_reference_length = 0;
    int is_coding_reference = line->cdsstart == line->cdsend ? 0 : 1;
    
    for (i=0; i<line->exoncount; ++i) {
	ss1 = ss;
	while(*ss1 && *ss1 != ',') ss1++;
	ss1[0] = '\0';
	line->exons[i].start = atoi(ss);
	ss = ++ss1;
	se1 = se;
	while(*se && *se1 != ',') se1++;
	se1[0] = '\0';
	line->exons[i].end = atoi(se);
	se = ++se1;
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
			dna_backward_reference_length += line->exons[i].end - line->cdsend+1;
		    }
		}
	    }
	}
    } // end for

    // dna_read_reference_length is the coding reference length for coding transcript or length of noncoding transcript 
    dna_read_reference_length = dna_reference_length - dna_forward_reference_length - dna_backward_reference_length;

    // Coding DNA reference:
    // -3                                            *3
    // -2                                          *2
    // -93   -45 -44       -1 1                187 188           351 *1   *96 *97    *223
    // |_______|____________|====================|=================|________|_______|
    // / \            ATG                / \             TGA        / \
    // /   \                             /   \                      /   \
    // /gtga,,,g                         /gta,,,act                 /gtag,,,cc
    // -45+1        -44-1                187+1        188-1         *96+1        *97-1
    // -45+2      -44-2                  187+2      188-2           *96+2      *97-2
    // -45+3    -44-3                    187+3    188-3             *96+3    *97-3
    // 
    // 
    // Noncoding DNA reference:
    // 
    // 1        49 50             280 281                 540 541             667
    // |__________|__________________|_______________________|_________________|
    // 
     
    int l1=0, l2=line->exoncount -1;
    // count forward
    int dna_read_reference_forward_offset = 0;
    int dna_read_reference_backward_offset = 0;
    for (; l1<line->exoncount && dna_forward_reference_length; ++l1) {
	int32_t length_exon = line->exons[l1].end - line->exons[l1].start;
	line->dna_ref_offsets[l1].start =
	    (dna_forward_reference_length << DNA_REF_OFFSETS_BITS) |
	    (line->strand == '+' ? FUNC_REGION_UTR5 : FUNC_REGION_UTR3);

	if (dna_forward_reference_length > length_exon) {
	    dna_forward_reference_length -= length_exon;		
	    line->dna_ref_offsets[l1].end =
		((dna_forward_reference_length+1) << DNA_REF_OFFSETS_BITS) |
		(line->strand == '+' ? FUNC_REGION_UTR5 : FUNC_REGION_UTR3);
	} else {
	    if (line->strand == '+') {
		dna_read_reference_forward_offset = length_exon - dna_forward_reference_length;
	    } else {
		dna_read_reference_forward_offset = dna_read_reference_length + dna_forward_reference_length - length_exon +1; // for minus strand count from backward
	    }
	    line->dna_ref_offsets[l1].end = (dna_read_reference_forward_offset<<DNA_REF_OFFSETS_BITS)| FUNC_REGION_CODING;
	    dna_forward_reference_length = 0;
	    ++l1;
	    break;
	}
    }

    // count backward
    for (; l2 > 0 && dna_backward_reference_length; --l2) {
	int32_t length_exon = line->exons[l2].end - line->exons[l2].start;
	line->dna_ref_offsets[l2].start =
	    (dna_backward_reference_length << DNA_REF_OFFSETS_BITS) |
	    (line->strand == '+' ? FUNC_REGION_UTR3 : FUNC_REGION_UTR5);
	if (dna_backward_reference_length > length_exon) {
	    dna_backward_reference_length -= length_exon;
	    line->dna_ref_offsets[l2].end =
		((dna_backward_reference_length+1) << DNA_REF_OFFSETS_BITS) |
		(line->strand == '+' ? FUNC_REGION_UTR3 : FUNC_REGION_UTR5);
	} else {
	    if (line->strand == '+') {
		dna_read_reference_backward_offset = dna_read_reference_length + dna_backward_reference_length - length_exon + 1;
	    } else {
		dna_read_reference_backward_offset = length_exon - dna_backward_reference_length;
	    }
	    line->dna_ref_offsets[l2].end = (dna_read_reference_backward_offset<<DNA_REF_OFFSETS_BITS) | FUNC_REGION_CODING;
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
		(dna_read_reference_length << DNA_REF_OFFSETS_BITS) |
		(is_coding_reference ? FUNC_REGION_CODING : FUNC_REGION_NONCODING);

	    dna_read_reference_length -= length_exon;

	    line->dna_ref_offsets[l].start = ((dna_read_reference_length+1)<< DNA_REF_OFFSETS_BITS) |
		(is_coding_reference ? FUNC_REGION_CODING : FUNC_REGION_NONCODING);
	}
    } else {
	if (is_coding_reference) {
	    dna_read_reference_length = dna_read_reference_forward_offset -1;
	}
	int l;
	for (l=l1; l<=l2; l++) {
	    int32_t length_exon = line->exons[l].end - line->exons[l].start;
	    
	    line->dna_ref_offsets[l].start = (dna_read_reference_length << DNA_REF_OFFSETS_BITS) |
		(is_coding_reference ? FUNC_REGION_CODING : FUNC_REGION_NONCODING);

	    dna_read_reference_length -= length_exon;
	    line->dna_ref_offsets[l].end = ((dna_read_reference_length+1) << DNA_REF_OFFSETS_BITS) |
		(is_coding_reference ? FUNC_REGION_CODING : FUNC_REGION_NONCODING);
	}
    }
}
void generate_dbref_database(struct gene_predictions_line *line)
{
    if (line->clear == 1) return;
    int i, j;
    for (i=0; i<line->exoncount; ++i) {
    	kstring_t temp[2] = {KSTRING_INIT, KSTRING_INIT};
    	int types[2];
    	types[0] = line->dna_ref_offsets[i].start & FUNC_REGION_MASK;
    	types[1] = line->dna_ref_offsets[i].end & FUNC_REGION_MASK;
    	int j;
    	for (j=0; j<2; ++j) {
    	    kstring_t *temp1 = &temp[j];
    	    int type = types[j];
    	    switch (type) {
    		case FUNC_REGION_UTR5 :
    		    kputc('-', temp1);
    		    break;
    		case FUNC_REGION_UTR3:
    		    kputc('*', temp1);
    		    break;
    		case FUNC_REGION_CODING:
    		    kputs("c.", temp1);
    		    break;
    		case FUNC_REGION_NONCODING:
    		    kputs("n.", temp1);
    		    break;
    		default:
    		    error("Unknown type : %d", type);
    	    }
    	}
    	kputw((line->dna_ref_offsets[i].start>>DNA_REF_OFFSETS_BITS), &temp[0]);
    	kputw((line->dna_ref_offsets[i].end>>DNA_REF_OFFSETS_BITS), &temp[1]);

	// format: CHROM,START,END,STRAND, GENE, TRANSCRIPT, EXON, START_LOC, END_LOC 
    	fprintf(
	    stderr, "%s\t%u\t%u\t%c\t%s\t%s\tEX%d\t%s\t%s\n",
	    line->chrom, line->exons[i].start, line->exons[i].end, line->strand,
	    line->name1, line->name2, i+1, temp[0].s, temp[1].s
	    );
    	free(temp[0].s);
    	free(temp[1].s);
    }
}

void push_mempool(struct gene_predictions_memory_pool *pool, kstring_t *str)
{
    if (pool->m == pool->l) {
	pool->m = pool->m == 0 ? 2 : pool->m << 1;
	pool->a = (struct gene_predictions_line *)realloc(pool->a, sizeof(struct gene_predictions_line)*pool->m);
    }
    // the prase func must return a point or go abort 
    gene_predictions_prase_line(str, &pool->a[pool->l]);
    pool->l++;
}
void update_mempool(struct gene_predictions_memory_pool *pool)
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
void fill_mempool(struct gene_predictions_memory_pool *pool, htsFile *fp, tbx_t *tbx, int tid, uint32_t start, uint32_t end)
{
    hts_itr_t *itr = tbx_itr_queryi(tbx, tid, start, end);
    kstring_t str = KSTRING_INIT;
    if (tid == -1) return;
    
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
	push_mempool(pool, &str);
	str.l = 0;
    }
    pool->tid = tid;
    if (str.m) free(str.s);
    tbx_itr_destroy(itr);
    update_mempool(pool);    
}
void clear_gene_predictions_line(struct gene_predictions_line *line)
{
    if (line->clear == 1) return;
    free(line->chrom);
    free(line->name1);
    free(line->name2);
    free(line->exons);
    free(line->dna_ref_offsets);
    line->txstart = line->txend = 0;
    line->cdsstart = line->cdsend = 0;
    line->exoncount = 0;
    line->clear = 1;
}
void clear_memory_pool(struct gene_predictions_memory_pool *pool)
{
    int i;
    for (i=0; i<pool->m; ++i)
	clear_gene_predictions_line(&pool->a[i]);
    pool->l = 0;
    pool->m = 0;
    free(pool->a);
    pool->a = NULL;
    pool->start = 0;
    pool->end = 0;
    pool->tid = -1;
}

// HGVS nomenclature : 
// DNA recommandations *
// - substitution variant, 
// Format: “prefix”“position_substituted”“reference_nucleotide””>”new_nucleotide”, e.g. g.123A>G
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
    var_type_ref,
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
struct hgvs_names_description *describe_variants(const char *ref, const char *alt, int _pos)
{
    struct hgvs_names_description *des = (struct hgvs_names_description*)malloc(sizeof(struct hgvs_names_description));
    const char *a = alt;
    const char *r = ref;
    int pos = _pos; // pos may differ from _pos, but _pos will not change in this function
    while (*a && *r && toupper(*a) == toupper(*r)) { a++; r++; pos++; }

    if (!a[0] && !r[0]) {
	des->type = var_type_ref;
	return des;
    }
    
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
	des->end = pos;
	des->alt_length = (a-alt) - (r-ref);
	des->alt = strndup(a-des->alt_length, des->alt_length);
	return des;
    } else if (!*a && *r) {
	des->type = var_type_dels;
	while ( *r ) r++;
	des->start = pos;
	des->ref_length = (r-ref) - (a->alt);
	des->end = pos + des->ref_length -1;
	des->alt_length = 0;
	des->alt = 0;
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
	des->end = pos + des->ref_length -1;
	des->alt_length = 0;
	des->alt = 0;
	des->ref = strndup(r, re-r);
	return des;
    } else if (re == r) {
	des->type = var_type_ins;
	des->start = pos;
	des->end = pos;
	des->ref_length = 0;
	des->ref = 0;
	des->alt_length = ae - a;
	des->alt = strndup(a, ae-a);
	return des;
    }

    var->type = var_type_delins;
    des->start = pos;
    des->ref_length = re-r;
    des->ref = strndup(r, re-r);
    des->alt_length = ae-a;
    des->alt = strndup(a, ae-a);
    des->end = pos + des->ref_length -1;    
    return des;
}
void check_cnv_hgvs_names_descriptions(struct hgvs_names_description *des, struct gene_predictions_line *line, faidx_t *fai, htsFile *fp)
{
    if (des->type == var_type_ins) {

    } else if (des->type == var_type_dels) {

    }
}

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
    uint8_t *data;
};

struct hgvs_names_compact {
    int l, m;
    struct hgvs_names_core *a;
};

struct simple_list {
    struct simple_list *next;
    void *data;
};
struct hgvs_names_core *hgvs_names_core_init(void)
{
    struct hgvs_names_core *c = (struct hgvs_names_core*)malloc(sizeof(struct hgvs_names_core));
    c->l_name = 0;
    c->l_type = 0;
    c->data = NULL;
    return c;
}
void clean_hgvs_names_core(struct hgvs_names_core *c)
{
    if (c->data != NULL)
	free(c->data);
    c->data = NULL;
    c->l_name = 0;
    c->l_type = 0;    
}

void clean_hgvs_names_compact(struct hgvs_names_compact *compact)
{
    int i;
    for (i=0; i<compact->l; ++i)
	clean_hgvs_names_core(&compact->a[i]);
    
    if (compact->m) free(compact->a);
    compact->m = 0;
    compact->l = 0;
}

void generate_hgvs_names_core(struct gene_predictions_line *line, struct hgvs_names_description *des, struct hgvs_names_core *c, int *valid)
{
    clean_hgvs_names_core(c);
    *valid = 0;
    if (des->start > line->txend || des->end < line->txstart) return;

    int exIn_id; // exon / intron id
    for (exIn_id = 0; exIn_id < 


    

    
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
// should be compact into one string like NM_xxxx:c.xx; NC_xxxx:n.xxx; ...; 
// @hgvs_names  list of hgvs names strings
// @n_names2    number of transcripts
// @names2      list of transcripts
struct hgvs_names_compact *generate_hgvs_names(struct gene_predictions_memory_pool *pool, bcf1_t *line, int *n_names1, char *names1, int *n_ale)
{
    *n_names1 = 0;
    *n_ale = line->n_allele -1;
    if (*n_ale == 0) return NULL;
    bcf_dec_t *d = &line->d;

    struct hgvs_names_compact * compact = (struct hgvs_names_compact*)calloc(*n_ale, sizeof(struct hgvs_names_compact));
    
    // a simple list structure inited for caching gene names and de-duplicate 
    struct simple_list *list = (struct simple_list*)malloc(sizeof(struct simple_list));
    struct simple_list *head = list;
    list->next = NULL;
    list->data = NULL;
    
    // roadmap : snv -> del -> ins -> dup ->delins
    int i, j;
    for (i=1; i<line->n_allele; ++i) {

	if (list->data == NULL) {
	    list->data = (void*)strdup(line->name1);
	} else {
	    struct simple_list **ll = &head;
	    while (*ll) {
		if (memcmp((char*)(*ll)->data, line->name1, strlen(line->name1))) *ll = (*ll)->next;
		else break;
	    }
	    if (*ll) {
		list->next = (struct simple_list*)malloc(sizeof(struct simple_list));
		list = list->next;
		list->data = (void*)strdup(line->name1);
		list->next = NULL;
	    }		    
	}
	struct hgvs_names_compact *c1 = &compact[i-1];
	struct hgvs_names_description *des = describe_variants(d->allele[0], d->allele[i], line->pos +1);
	
	for (j=0; j<pool->l; ++j) {
	    int valid = 0;
	    if (c1->l == c1->m) {
		c1->m = c1->m == 0 ? 2 : c1->m <<1;
		c1->a = (struct hgvs_names_core*)realloc(c1->a, sizeof(struct hgvs_names_core)*c1->m);
	    }
	    struct hgvs_names_core *c = &c1->a[c1->l];
	    generate_hgvs_names_core(&pool->a[i], des, c, &valid);
	    if (valid == 0) continue;
	    c1->l++;
	}
    }
    list = head;
    kstring_t str = KSTRING_INIT;
    int n = 0;
    while (list) {
	kputs(list->data, &str);
	kputc(',', &str);
	n++;	    
	head = list;
	list = list->next;
	free(head->data);
	free(head);
    }
    *n_names1 = n;
    names1 = str.s;         
}
// number : R, type : string
static void setter1_hgvs_names(bcf_hdr_t *hdr, bcf1_t *line, char *string) 
{
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, "HGVSDNA");
    if (id == -1) {
	bcf_hdr_append(hdr, "##INFO=<ID=HGVSDNA,Number=G,Type=String,Description=\"HGVS nomenclature for the description of DNA sequence variants.\">");
	bcf_hdr_sync(hdr);
    }
    // the data construct along with alt alleles, so there is no need to remap the order of alleles
    bcf_update_info_string(hdr, line, "HGVSDNA", string);
}
// number : 1, type : string
static void setter1_gene_names(bcf_hdr_t *hdr, bcf1_t *line, char *string)
{
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, "Gene");
    if (id == -1) {
	bcf_hdr_append(hdr, "##INFO=<ID=Gene,Number=1,Type=String,Description=\"Gene names.\">");
	bcf_hdr_sync(hdr);
    }
    bcf_update_info_string(hdr, line, "Gene", string);
}
// number : R, type : string
static void setter1_transcripts_names(bcf_hdr_t *hdr, bcf1_t *line, char *string)
{
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, "Transcript");
    if (id == -1) {
	bcf_hdr_append(hdr, "##INFO=<ID=Transcript,Number=G,Type=String,Description=\"Transcript names.\">");
	bcf_hdr_sync(hdr);
    }
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
char *generate_transcript_string(struct hgvs_names_compact *compact, int n)
{
    kstring_t str = KSTRING_INIT;
    int i, j;    
    for (i=0; i<n; ++l) {
	if (i) kputc(',', &str);
	// foreach allele
	struct hgvs_names_compact *compact1 = &compact[i];
	for ( j = 0; j < compact1->l; ++j ) {
	    struct hgvs_names_core *line = &compact1->a[j];
	    if (line->l_name) kputsn(line->data, line->l_name-1, &str);
	    else kputc('.', &str);
	    kputs("; ", &str);
	}
    }
    return str.s;
}
char *generate_hgvsvarnomen_string(struct hgvs_names_compact *compact,int n)
{
    kstring_t str = KSTRING_INIT;
    int i, j;    
    for (i=0; i<n; ++l) {
	if (i) kputc(',', &str);
	// foreach allele
	struct hgvs_names_compact *compact1 = &compact[i];
	for ( j = 0; j < compact1->l; ++j ) {
	    struct hgvs_names_core *line = &compact1->a[j];
	    if (line->l_name) kputs(line->data, &str);
	    else kputc('.', &str);
	    kputs("; ", &str);
	}
    }
    return str.s;
}
// anno_hgvs_core() only used to annotate bcf/vcf standalone, all the valid tags-added functions will be called in 
// this function, to annotate the vcf more specification, use setter_genepred_* functions instead of it. 
void anno_hgvs_core(struct gene_predictions_memory_pool *pool, htsFile *fp, tbx_t *tbx, bcf_hdr_t *hdr, bcf1_t *line)
{
    // retrieve the regions this variants located first, check the memory pool and update the pool if the regions is out of
    // cached positions. just skip if the variant type of line is a ref. 
    if (line->pos <= 0) return;
    if (bcf_get_variant_types(line) == VCF_REF) return;

    uint32_t end = bcf_calend(line);
    int id = tbx_name2id(tbx, bcf_hdr_id2name(hdr, line->rid));
    // fill mempool, skip if no record in the memory pool 
    fill_mempool(pool, fp, tbx, id, line->pos, end);
    if (pool->l == 0) return;
    int n_names1 = 0;
    char *names1 = NULL;
    int n_ale = 0;
    
    struct hgvs_names_compact *c = generate_hgvs_names(pool, line, &n_names1, names1, &n_ale);
    if (n_ale == 0) return;
    char *trans_string = generate_transcript_string(c, n_ale);
    char *hgvs_string = generate_hgvsvarnomen_string(c, n_ale);
    
    setter_hgvs_string(hdr, line, "Transcript", trans_string);
    setter_hgvs_string(hdr, line, "HGVSDNA", hgvs_string);
    setter_hgvs_string(hdr, line, "Gene", names1);
    int i;
    for ( i = 0; i < n_ale; ++i )
	clean_hgvs_names_compact(&c[i]);
    free(c);
}

#include <sys/time.h>
static double get_time() {
    struct timeval tv;
    if ( gettimeofday(&tv, 0) != 0 ) 
	error("Failed to get time of day.");
    return (double)tv.tv_sec + 1.0e-6// (double)tv.tv_usec;
}

struct args {
    const char *gene_predictions_file;
    const char *refseq_file;
    const char *out_file;
    const char *out_type;
    const char *input_fname;
    tbx_t *gpidx;
    faidx_t *faidx;
    glist_t *glist;
    struct gene_predictions_memory_pool pool;
};

void clear_args(struct args *args)
{
    tbx_destroy(args->gpidx);
    if (args->faidx) fai_destroy(args->faidx);
    if (args->glist) kh_destroy(list, args->glist);
    clear_memory_pool(&args->pool);
}
int usage(int argc, char **argv)
{
    fprintf(stderr, "Usage: %s [-h] -refseq refseq.fa.gz -data gene_predictions.txt.gz|refgene.txt.gz -O [u|v|z|b] -o output_file input_fname.vcf.gz \n", argv[0]);
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
    static struct args args = {
	.gene_predictions_file = 0,
	.refseq_file = 0,
	.out_file = 0,
	.out_type = 0,
	.input_fname = 0,
	.gpidx = 0,
	.faidx = 0,
	.glist = 0,
    };
        
    for (int i = 1; i< argc;) {
	const char *a = argv[i++];
	if (strcmp(a, "-h") == 0) 
	    return usage(argc, argv);

	const char **arg_var = 0;
	if (strcmp(a, "-o") == 0 && args.out_file == 0)
	    arg_var = &args.out_file;
	else if (strcmp(a, "-O") == 0 && args.out_type == 0)
	    arg_var = &args.out_type;
	else if (strcmp(a, "-data") == 0 && args.gene_predictions_file == 0)
	    arg_var = &args.gene_predictions_file;
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

    if (args.gene_predictions_file == 0)
	error("Reasons :\n"
	      "-data gene preditions database is needed.\n"
	      "You could download refGene.txt.gz or ensGene.txt.gz from UCSC websites and sort and indexed by tabix.");

    // if (args.refseq_file == 0) 
    // 	error("Reasons :\n" 
    // 	      "-refseq refseq.fa.gz file is need.\n" 
    // 	      "This is the transcripts reference sequences in fasta format." 
    // 	      "Notice the transcripts names should be consistance with gene preditions databases."); 

    htsFile *fgp = hts_open(args.gene_predictions_file, "r");
    if (fgp == NULL)
	error("Failed to open %s", args.gene_predictions_file);
    args.gpidx = load_gene_predictions_file(args.gene_predictions_file);
    if (args.gpidx == NULL)
	error("Failed to load index of %s", args.gene_predictions_file);

    // args.faidx = load_refseq_file(args.refseq_file); 
    // if (args.faidx == NULL) 
    // 	error("Failed to load index of %s", args.refseq_file); 
    
    // assuming input file is stdin, only vcf, gzipped vcf, bcf file is accept,
       err msg will output if file type unrecongnized
    if (args.input_fname == 0 && (!isatty(fileno(stdin)))) args.input_fname = "-";
    debug_print("input_fname : %s", args.input_fname);
    htsFile *fp = hts_open(args.input_fname, "r");
    if (fp == NULL)
	error("Failed to open %s.", args.input_fname);
    htsFormat type = *hts_get_format(fp);
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
    struct gene_predictions_memory_pool *mempool = &args.pool;
    mempool->m = mempool->l = 0;
    mempool->tid = -1;
    mempool->start = mempool->end = 0;
    mempool->a = NULL;
    
    double c0 = get_time();    
    bcf_hdr_t *hdr = bcf_hdr_read(fp);    

    // duplicate header file and sync new tag in the output header .
    // assuming output is stdout in Vcf format in default. 
    bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);	
    htsFile *fout = args.out_file == 0 ? hts_open("-", hts_bcf_wmode(out_type)) : hts_open(args.out_file, hts_bcf_wmode(out_type));
    
    // init gene_predictions or refgene database, hold tabix index cache and faidx cache in memory 
    bcf1_T *line = bcf_init();
    while ( bcf_read(fp, hdr, line) == 0 ) {     
	anno_hgvs_core(mempool, fgp, args.gpidx, hdr, line);	
	bcf_write1(fout, hdr_out, line);
    }

    clear_args(&args);
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    bcf_hdr_destroy(hdr_out);
    hts_close(fout);
    hts_close(fp);
    double c1 = get_time();
    fprintf(stderr, "Run time: %.2fs\n", c1 -c0);
    return 0;
}
