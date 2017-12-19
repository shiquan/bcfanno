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
#include "sequence.h"
// variants in genome reference positions in raw vcf file,
// recommand hgvs nomenclature to describe variants,]
// http://varnomen.hgvs.org
char str_init[2];

KHASH_MAP_INIT_STR(name, char*)

typedef khash_t(name) namehash_type;

typedef char *(*copy_seqs_func)(const char*, unsigned long n);

// lite structure for cache refseq 
static struct refseq_cache {
    char *name;
    int start;
    int end;
    int l_seq;
    faidx_t *idx;
    char *seq;
} refseq_cache = {
    .name = NULL,
    .start = 0,
    .end = 0,
    .l_seq = 0,
    .idx = NULL,
    .seq = NULL,
};
                      
static void hgvs_core_clear(struct hgvs_core *c)
{
    c->str.l = 0;
    c->l_name = 0;
    c->l_type = 0;
    memset(&c->type, 0, sizeof(struct var_func_type));
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

    for ( i = 0; i < n; ++i ) {
       k = kh_put(name, hash, names[i], &ret);
    }
    return hash;
}
faidx_t *load_refseq_file(const char *fn)
{
    return fai_load(fn);
}
int32_t bcf_calend(bcf1_t *line)
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
	bcf_hdr_append(hdr, "##INFO=<ID=Transcript,Number=1,Type=String,Description=\"Transcript names\">");
	bcf_hdr_sync(hdr);
	id = bcf_hdr_id2int(hdr, BCF_DT_ID, "Transcript");
	assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));
    }
    return id;
}
int hgvs_bcf_header_add_hgvsdna(bcf_hdr_t *hdr)
{
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, "HGVSDNA");
    if (id == -1) {
	bcf_hdr_append(hdr, "##INFO=<ID=HGVSDNA,Number=A,Type=String,Description=\"HGVS nomenclature for the description of DNA sequence variants\">");
	bcf_hdr_sync(hdr);
	id = bcf_hdr_id2int(hdr, BCF_DT_ID, "HGVSDNA");
	assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));
    }
    return id;
}
// function regions
// Intergenic, Promoter, UTR5, CDS{n}, Intron{n}, UTR3, Noncoding Exon{n}
int hgvs_bcf_header_add_exid(bcf_hdr_t *hdr)
{
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, "ExIn_id");
    if (id == -1) {
	bcf_hdr_append(hdr, "##INFO=<ID=ExIn_id,Number=1,Type=String,Description=\"Exon or intron id on transcripts.\">");
	bcf_hdr_sync(hdr);
	id = bcf_hdr_id2int(hdr, BCF_DT_ID, "ExIn_id");
	assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));
    }
    return id;
}
// variantion types
// synonymous, missense, inframe insertion, inframe deletion, stop gained, stop lost, stop retained, splice donor, splice acceptor, referenece
// ref: http://asia.ensembl.org/info/genome/variation/consequences.jpg
int hgvs_bcf_header_add_vartype(bcf_hdr_t *hdr)
{
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, "VarType");
    if (id == -1) {
	bcf_hdr_append(hdr, "##INFO=<ID=VarType,Number=A,Type=String,Description=\"Variant type.\">");
	bcf_hdr_sync(hdr);
	id = bcf_hdr_id2int(hdr, BCF_DT_ID, "VarType");
	assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));
    }
    return id;
}
char *describe_description(struct hgvs_des *des)
{
    kstring_t string = KSTRING_INIT;
    kputw(des->start, &string);
    kputc(',', &string);
    kputw(des->end, &string);
    return string.s;
}

static int hgvs_fill_memory(struct refgene_options *opts, bcf1_t *line)
{
    // get the chromosome id in tabix indexed genepred database
    assert(opts->genepred_idx);
    int rid = -1;
    const char *name = bcf_hdr_id2name(opts->hdr_out, line->rid);
    rid = tbx_name2id(opts->genepred_idx, name);
    // no chromosome found, return 1
    if (rid == -1)
	return 1;    
    struct gp_mempool *pool = &opts->buffer;
    int32_t end = bcf_calend(line);
    if (pool->rid == rid && line->pos+1 >= pool->start && end <= pool->end)
    	return 0;
    gp_clear_mempool(pool);
    hts_itr_t *itr = tbx_itr_queryi(opts->genepred_idx, rid, line->pos, end);

    kstring_t string = KSTRING_INIT;    
    while ( tbx_itr_next(opts->fp, opts->genepred_idx, itr, &string) >= 0) {
	gp_push_mempool(pool, &string);
	string.l = 0;
    }
    pool->rid = rid;
    if ( string.m )
	free(string.s);
    tbx_itr_destroy(itr);
    gp_update_mempool(pool);
    return 0;
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

int refgene_options_destroy(struct refgene_options *opt)
{
    if ( opt->refgene_is_inited == 0 )
        return 1;
    hts_close(opt->fp);
    tbx_destroy(opt->genepred_idx);
    if ( opt->check_refseq == 1 )
	fai_destroy(opt->refseq_fai);
    if ( opt->screen_by_genes == 1 )
	namehash_destroy(opt->genehash);
    if ( opt->screen_by_transcripts == 1 )
	namehash_destroy(opt->transhash);
    gp_clear_mempool(&opt->buffer);
    hgvs_cache_clear(&opt->cache);
    if ( refseq_cache.name )
        free(refseq_cache.name);
    if ( refseq_cache.seq )
        free(refseq_cache.seq);
    return 0;
}
static void hgvs_des_destory(struct hgvs_des *des)
{
    if (des->type == var_type_ref || des->type == var_type_nonref)
	return;
    if ( des->ref_length > 0 )
	free(des->ref);
    if ( des->alt_length > 0 )
	free(des->alt);
    free(des);
}

static void hgvs_des_reverse(struct hgvs_des *des, uint8_t strand)
{
    if ( des->type == var_type_ref || des->type == var_type_nonref)
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
    // for GATK users, <NON_REF> allele will be treat as ref
    if ( strcmp(alt, "<NON_REF>") == 0 ) {
      des->type = var_type_nonref;
      return des;
    }
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
    if (ae == a && re == r) {
        des->type = var_type_snp;
        des->start = pos;
        des->ref_length =  1;
        des->ref = strndup(r, 1);
        des->alt_length = 1;
        des->alt = strndup(a, 1);
        return des;
    }
    // check a and e in first step, so "re==r && ae == a" would not happen here 
    if ( ae == a) {
	des->type = var_type_dels;
	des->start = pos;
	des->ref_length = re -r + 1;
	des->ref = strndup(r, des->ref_length);
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
	des->alt_length = ae - a+1;
	des->alt = strndup(a, des->alt_length);
	return des;
    }
    // delins
    des->type = var_type_delins;
    des->start = pos;
    des->ref_length = re-r+1;
    des->ref = strndup(r, des->ref_length);
    des->alt_length = ae-a+1;
    des->alt = strndup(a, des->alt_length);
    des->end = pos + des->ref_length -1;
    return des;
}
// todo : check indels is in trf regions or not
static int check_repeat(struct hgvs_des *des, struct genepred *line, faidx_t *idx, htsFile *fp)
{
    if (des->type == var_type_ins) {

    } else if (des->type == var_type_dels) {

    }
    return 0;
}

static int indel2repeat(struct hgvs_des *des, struct genepred *line, faidx_t *idx, htsFile *fp)
{
    return check_repeat(des, line, idx, fp);
}

static char *retrieve_refseq_block(char *name, int start, int end)
{
    if ( refseq_cache.name == NULL )
        goto update_name;

    if ( strcmp(name, refseq_cache.name) == 0) {
        if (start == refseq_cache.start && end == refseq_cache.end)
            return refseq_cache.seq;
        else
            goto update_block;
    }
    free(refseq_cache.name);
    
  update_name:
    refseq_cache.name = strdup(name);

  update_block:
    if ( refseq_cache.seq )
        free(refseq_cache.seq);        
    refseq_cache.seq = faidx_fetch_seq(refseq_cache.idx, name, start-1, end -1, &refseq_cache.l_seq);
    if (refseq_cache.l_seq == 0 || refseq_cache.seq == NULL)
        return NULL;
    refseq_cache.start = start;
    refseq_cache.end = end;
    return refseq_cache.seq;    
}
// 
static int describe_vartype(struct genepred *line, int l1, int l2, struct hgvs_des *des, struct hgvs_core *c)
{
    struct var_func_type *type = &c->type;    
    int is_strand = line->strand == '+';
    type->start_flag = is_strand ? l1 : line->exoncount * 2 - l2 -1;
    
    if ( c->type.func == func_region_intron) {
        type->vartype = var_is_intron;
        return 0;
    }
    // for noncoding
    if ( line->cdsstart == line->cdsend ) {
        type->vartype = var_is_noncoding;
        return 0;
    }
    

    int upstream_offset = 0;
    int downstream_offset = 0;
    uint32_t loc_start = 0, loc_end = 0;
    
    if ( is_strand ) {
        upstream_offset = des->start - read_start(line->exons, l1/2);
        downstream_offset = read_end(line->exons, l2/2) - des->start;
        loc_start = read_start(line->loc, l1/2);
        loc_end  = read_end(line->loc, l2/2);
    } else {
        upstream_offset = read_end(line->exons, l2/2) - des->start;
        downstream_offset = des->start - read_start(line->exons, l1/2);
        loc_start  = read_end(line->loc, l2/2);
        loc_end  = read_start(line->loc, l1/2);
    }    

    assert(upstream_offset>=0 && downstream_offset >= 0);
    assert(loc_end > loc_start);

    // for intron block
    if ( type->start_flag & 1) {
        if ( upstream_offset < 3)
            type->vartype2 = var_is_splice_donor;
        else if ( downstream_offset < 3)
            type->vartype2 = var_is_splice_acceptor;
        //else
        type->vartype = var_is_intron;
        return 0;
    }
   
    int var_loc = upstream_offset + loc_start;    
    int coding_end_loc = line->reference_length - line->backward_length;
    
    if ( var_loc <= line->forward_length ) {
        type->vartype = var_is_utr5;
        return 0;
    }
    if ( var_loc > coding_end_loc ) {
        type->vartype = var_is_utr3;
        return 0;
    }
    
    int coding_start = 0, coding_end = 0;
    if ( loc_start > line->forward_length ) {
        var_loc -= loc_start;
    } else {
        loc_start = line->forward_length + 1;
    }
    // coding_start =  loc_start > line->forward_length ? loc_start : line->forward_length + 1;
    coding_end = loc_end > coding_end_loc ? coding_end_loc : loc_end;
    int block_length = coding_end - coding_start + 1;
    char *seq = retrieve_refseq_block(line->name1, coding_start, coding_end);
    if (seq == NULL)
        return 1;
    type->vartype = check_var_type(seq, block_length, var_loc, des->ref, des->ref_length, des->alt, des->alt_length);
    return 0;
}

// purpose: find the most nearest edge for variant position
// init:   l1                     l2
//         |----|   |----|  |-----|  
// End:         |  .|
//             l1   l2
// for iter l1 and l2, if even iter come from start array, if odd iter come from end array        
// 
// init          l1
//       start   |       |       |
//       end     |       |       |
//                               l2
// result:               l2
//               |       |       |
//               |       |       |
//                       l1
static int generate_hgvs_location(struct genepred *line, struct hgvs_des *des, struct hgvs_core *c)
{

#define POS_LOCATE(_l1, _l2, pos) do {\
        int l;\
        _l1 = 0;\
        _l2 = line->exoncount*2 -1;\
        while ( _l1 < _l2 ) {\
            if ( _l2 - _l1 == 1)\
                break;\
            l = _l1 + _l2;\
            l = l/2 + (l&1);\
            uint32_t ip = line->exons[l&1][l/2];  \
            if ( ip > pos)\
                _l2 = l;\
            else\
                _l1 = l;\
        }\
    } while(0)
    
#define PRINT_LOC(_t, _s, _loc, _func) do {     \
        if ( _t & REG_NONCODING ) {\
            _func = func_region_noncoding;\
        } else if ( _t & REG_UTR5) { \
            _func = func_region_utr5;\
            kputc('-', _s);\
        } else if ( _t & REG_UTR3) {\
            _func = func_region_utr3;\
            kputc('*', _s);\
        }\
        kputw(_loc, _s);\
    } while (0)

    int l1, l2;
    // [left|right]_loc is the utr/cds location of exon block edge, should always greater than 0.    
    int left_loc, right_loc;
    // [left|right]_offset is the length between position and block edge, should always greater than or equal to 0.
    int left_offset, right_offset; 
    // [left|right]_type is the function type of block edge, should be NONCODING, CODING, UTR3 and UTR5
    uint8_t left_type, right_type; 
    // location and type of the variant
    int loc = 0;
    uint8_t type = 0;

    int i = 0;
    kstring_t temp = KSTRING_INIT;

    for ( ; i < 2; i++) {
        uint32_t pos = i == 0 ? des->start : des->end;
        // find the block which variant located
        POS_LOCATE(l1, l2, pos);
        assert(l2 - l1 == 1);
        // left and right offset are the length between variant and near block edges
        left_offset = pos - line->exons[l1&1][l1/2];
        right_offset = line->exons[l2&1][l2/2] - pos;

        left_loc = line->dna_ref_offsets[l1&1][l1/2] >> TYPEBITS;
        right_loc = line->dna_ref_offsets[l2&1][l2/2] >> TYPEBITS;
        left_type = line->dna_ref_offsets[l1&1][l1/2] & REG_MASK;
        right_type = line->dna_ref_offsets[l2&1][l2/2] & REG_MASK;
    
        assert (left_offset >=0 || right_offset >=0);       

        // l1/2 should be the start of the function region, l1&2 == 1 for intron, for exon/cds/utr
        // regions l1&2 should always be 0
        if ( i == 1 )
            kputc('_', &temp);
        
        enum func_region_type func =  func_region_unknown;
        
        if ( l1 & 1 ) { // intron
            int offset;
            if ( left_offset > right_offset ) {
                loc = right_loc;
                type = right_type;                
                offset = line->strand == '+' ? -1 * right_offset : right_offset;
            } else {
                loc = left_loc;
                type = left_type;                
                offset = line->strand == '+' ? left_offset : -1 * left_offset;
            }
            PRINT_LOC(type, &temp, loc, func);
            if ( offset > 0 )
                kputc('+', &temp);
            kputw(offset, &temp);
            func = func_region_intron;
        } else {
            // if pos in exon
            // pos_[start, end] is the locs of near edges,
            //    pos_start pos       pos_end
            //    |         |         |
            //    |-------------------|
            //    offset_start offset_end
            func = func_region_cds;
            type = left_type;
            if ( left_type == right_type ) {
                loc = right_loc > left_loc ? right_loc - right_offset : left_loc - left_offset;
            } else {
                // if the block cross two different function regions, like utr5 to cds1

                // check the type of block edges
                if ( !(left_type & (REG_UTR3 | REG_UTR5 | REG_CODING)) )
                    error("Unknown type. left_type : %d", left_type);
                if ( !(right_type & (REG_UTR3 | REG_UTR5 | REG_CODING)) )	
                    error("Unknown type. right_type : %d", right_type);
                
                if ( left_type & REG_UTR5 ) {
                    // should only happens for plus strand
                    loc = left_loc - left_offset;
                    if ( loc <= 0 ) {
                        type = right_type;
                        if ( right_type & REG_CODING ) {
                            loc = 1 - loc;
                        } else if ( right_type & REG_UTR3 ) {
                            // all the cds region enclosed in one exon, closed with UTRs
                            int loc1 = right_loc - right_offset;
                            if ( loc1 <= 0 ) {
                                loc = 1 - loc;
                                type = REG_CODING;
                            } else {
                                loc = loc1;
                                type = right_type;
                            }                        
                        }
                    }
                } else if ( left_type & REG_CODING ) {

                    if ( right_type & REG_UTR3 ) {
                        // plus strand
                        loc = right_loc - right_offset;
                        if ( loc > 0 ) {
                            type = right_type;
                        } else {
                            loc = left_loc + left_offset;
                        }

                    } else if ( right_type & REG_UTR5 ) {
                        // minus strand
                        loc = left_loc - left_offset;
                        if ( loc <= 0 ) {
                            loc = 1 - loc;
                            type = right_type;
                        }
                    }
                    
                } else if ( left_type & REG_UTR3 ) {
                    // should only happened in minus strand
                    loc = left_loc - left_offset;
                    if ( loc <= 0 ) {
                        if ( right_type & REG_CODING) {
                            loc = 1 - loc;
                            type = right_type;
                        } else if (right_type & REG_UTR5) {
                            int loc1 = right_loc - right_offset;
                            if ( loc1 <= 0 ) {
                                loc = 1 - loc;
                                type = REG_CODING;
                            } else {
                                loc = loc1;
                                type = right_type;
                            }
                        }
                    }
                } // end left type is UTR3
            } // end differnt types of start and end of block      
            PRINT_LOC(type, &temp, loc, func);
        }

        if ( i == 0 ) {
            c->type.func = func;
            if ( describe_vartype(line, l1, l2, des, c) )
                return 1; // no such transcript in refMrna.fa
        }

        if ( des->type == var_type_snp || (des->type == var_type_dels && des->ref_length == 1 ))
            i = 2; // skip if the variant is snv        
    }
    
#undef PRINT_LOC
#undef POS_LOCATE        

    kstring_t *string = &c->str;
    // generate coding and noncoding flag    
    if ( line->cdsstart == line->cdsend)
        kputs("n.", string);
    else
        kputs("c.", string);
    
    c->l_type = string->l;
    
    // generate location
    kputs(temp.s, string);
    free(temp.s);
    return 0;   
}

// generate_hgvs_core1 - convert genome coordinate to hgvs position based on the aligned transcripts regions.
// here will check the start position of the variant first, in most of cases, it is no necessary to check
// the end position, because most of the variants only affect one position. however, in some complex variants
// like deletions and delins, check the end position will help the research understand this variant better
// To describe the function regions of the variants, here only check the start position for simplity. For some
// huge deletions may affect more than one exon or function region, please refer to bedanno() to get all the
// function regions related.
static int generate_hgvs_changes(struct genepred *line, struct hgvs_des *des, struct hgvs_core *c)
{
    kstring_t *string = &c->str;

    // if no-change
    if (c->type.vartype == var_is_reference) {
        kputc('=',string);
        return 0;
    }
    // generate variant changes
    char *ref = strndup(des->ref, des->ref_length);
    char *alt = strndup(des->alt, des->alt_length);
    
    if ( des->type == var_type_snp) {
        ksprintf(string, "%s>%s", ref, alt);
    } else if ( des->type == var_type_dels) {
        // for one base deletion, format like NM_0001:c.123del
        if ( des->ref_length == 1 ) 
            ksprintf(string, "del%s", des->ref);
        else
            kputs("del", string);
    } else if ( des->type == var_type_ins ) {
        ksprintf(string, "ins%s", alt);
    } else if ( des->type == var_type_delins ) {
        // delins format should be NM_0001:c.123_125delinsXX
        ksprintf(string, "delins%s", alt);
    } else {
        error("Unknown type : %d.", des->type);
    }
    // memory release
    if ( des->ref_length )
        free(ref);
    if ( des->alt_length )
        free(alt);
    return 0;
}
static int set_empty_hgvs(struct hgvs *name )
{
    if ( name->m == 0 ) {
        name->m = 1;
        name->a = (struct hgvs_core*)calloc(1, sizeof(struct hgvs_core));
    }
    kstring_t *string = &name->a[0].str;
    if ( name->i == 0 ) {    
        string->l = string->m = 0;
        string->s = NULL;
        name->i = 1;
    }
    name->l = 1;
    hgvs_core_clear(&name->a[0]);
    kputc('.', string);
    return 0;
}
// return 0 on success, including ref pos, 1 on failure
static int generate_hgvsname1(struct genepred *line, struct hgvs_des *des, struct hgvs_core *c)
{
    hgvs_core_clear(c);
    kstring_t *string = &c->str;
  
    // out of range
    if (des->start > line->txend || des->start <= line->txstart)
	return 1;

    // put names into string, usually transcript 
    // if no transcript name, use gene name instead
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
	c->l_name2 = string->l;
    }
    if (c->l_name == 0) 
	error ("Failed to parse genepred database. %s\t%u\t%u", line->chrom, line->txstart, line->txend);
    kputc(':', string);

   
    if ( generate_hgvs_location(line, des, c) ) {
        warnings("Failed to parse hgvs pos in genepred record. %s vs %s", describe_description(des), line->name1);
        return 1;
    }
    
    generate_hgvs_changes(line, des, c);    
    
    return 0;
}

// -- assuming type of hgvs.name tag is string and number is R. this is mandontary for VCF file 
// genepred_memory        cached memory pool for gene predictions records
// line        input variant record in bcf1_t format
// hgvs_cache will init or update in this function
static void generate_hgvs(struct refgene_options *opts, bcf1_t *line)
{
    struct gp_mempool *buffer = &opts->buffer;
    if ( buffer->l == 0 )
      return;
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
        if ( des->type == var_type_nonref || des->type == var_type_ref) {
            set_empty_hgvs(name);
            hgvs_des_destory(des);
            continue;
        }
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
	    struct genepred *gl = &buffer->a[j];
	    // check strand of transcript and reverse the description if it is minus strand
	    hgvs_des_reverse(des, gl->strand);
	    // generate hgvs name core
            if ( des->type != var_type_snp )
                indel2repeat(des, gl, opts->refseq_fai, opts->fp);
            
	    if ( generate_hgvsname1(gl, des, core) != 0)
		continue;
            
	    name->l++;
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
static void setter1_exin_id(bcf_hdr_t *hdr, bcf1_t *line, char *string)
{
    bcf_update_info_string(hdr, line, "ExIn_id", string);
}
static void setter1_vartype(bcf_hdr_t *hdr, bcf1_t *line, char *string)
{
    bcf_update_info_string(hdr, line, "VarType", string);
}
//static void setter1_flankseq(bcf_hdr_t *hdr, bcf1_t *line, char *string)
//{
//    bcf_update_info_string(hdr, line, "FLKSEQ", string);
//}
static void setter_hgvs_string(bcf_hdr_t *hdr, bcf1_t *line, const char *key, char *string)
{
    if (key == NULL)
        return;
    if (strcmp(key, "Gene") == 0) {
	setter1_gene_names(hdr, line, string);	
    } else if (strcmp(key, "HGVSDNA") == 0) {
	setter1_hgvs(hdr, line, string);		
    } else if (strcmp(key, "Transcript") == 0) {
	setter1_transcripts_names(hdr, line, string);
    } else if (strcmp(key, "ExIn_id") == 0 ) {
        setter1_exin_id(hdr, line, string);
    } else if (strcmp(key, "VarType") == 0 ) {
        setter1_vartype(hdr, line, string);
        //  } else if (strcmp(key, "FLKSEQ") == 0 ) {
//        setter1_flankseq(hdr, line, string);
    } else {
	error("Unrecongnized tag %s, only Gene, HGVSDNA, and Transcript are supported for now.", key);
    }
}
static char *generate_transcript_string(struct hgvs_cache *cache)
{
    kstring_t str = KSTRING_INIT;
    int i, j;    
    for ( i = 0; i < cache->l; ++i ) {
	/* if ( i ) */
	/*     kputc(',', &str); */
        // foreach allele
        struct hgvs *hgvs = &cache->a[i];
        // if hgvs is empty, skip
        if ( hgvs->l == 0 || (hgvs->l == 1 && hgvs->a[0].str.s[0] == '.')) 
            continue;	
	for ( j = 0; j < hgvs->l; ++j ) {
            if ( j )
                kputc('|', &str);
	    struct hgvs_core *core = &hgvs->a[j];
	    if ( core->l_name )
		kputsn(core->str.s, core->l_name, &str);
	    else
		kputc('.', &str);
	}
        break;
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
            if ( j )
                kputc('|', &str);
	    struct hgvs_core *core = &name->a[j];
	    if (core->l_name)
		kputs(core->str.s, &str);
	    else
		kputc('.', &str);
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
        // if hgvs is empty, skip
        if ( hgvs->l == 0 || (hgvs->l == 1 && hgvs->a[0].str.s[0] == '.')) 
            continue;        
	for (j = 0; j < hgvs->l; ++j) {
            
	    struct hgvs_core *core = &hgvs->a[j];
	    if (core->l_name  && core->l_name2 > 0) {
		kstring_t *string = &core->str;
		assert(string->s[core->l_name]  == '(');
		kputsn(string->s + core->l_name + 1, core->l_name2 - core->l_name -2, &str);
		return str.s;
	    }
	}
    }
    return NULL;
}
static char *generate_vartype_string(struct hgvs_cache *cache)
{
    int i, j;
    kstring_t str = KSTRING_INIT;
    for ( i = 0; i < cache->l; ++i ) {
        if ( i )
            kputc(',', &str);
        struct hgvs *hgvs = &cache->a[i];
        for ( j = 0; j < hgvs->l; ++j ) {
            if ( j )
                kputc('|', &str);
            struct var_func_type *type = &hgvs->a[j].type;
            kputs(var_type_string(type->vartype), &str);
        }
    }
    return str.s;
}
static char *generate_exin_ids(struct hgvs_cache *cache)
{
    int i, j;
    kstring_t str = KSTRING_INIT;
    for (i = 0; i < cache->l; ++i) {
        struct hgvs *hgvs = &cache->a[i];
        // if hgvs is empty, skip
        if ( hgvs->l == 0 || (hgvs->l == 1 && hgvs->a[0].str.s[0] == '.')) 
            continue;
        for (j=0; j < hgvs->l; ++j) {
            if ( j )
                kputc('|', &str);
            struct var_func_type *type = &hgvs->a[j].type;
            if ( type->func == func_region_intron ) {
                assert(type->start_flag&1);
                kputs("Intron", &str);
                kputw(type->start_flag/2 + 1, &str);
            } else if (type->func == func_region_intergenic ) {
                kputs("Intergenic", &str);
            } else {
                kputs("Exon", &str);
                kputw(type->start_flag/2 + 1, &str);
                if ( type->func == func_region_utr5)
                    kputs("(UTR5)", &str);
                else if ( type->func == func_region_utr3)
                    kputs("(UTR3)", &str);
            }
        }
        break;
    }
    return str.s;
}

int hgvs_names_init(struct refgene_options *opts, bcf1_t *line)
{
    struct hgvs_cache *cache = &opts->cache;
    // check the hgvs names inited already ?
    hgvs_cache_clear(cache);
    // fill memory pool, only if no transcriptc cover this variant return 1, else return 0
    if ( hgvs_fill_memory(opts, line) != 0 )
	return 1;    

    generate_hgvs(opts, line);
    
    return 0;
}

int anno_refgene_core(struct refgene_options *opts, bcf1_t *line)
{
    // if refgene options not inited, skip all the next steps
    if (opts->refgene_is_inited == 0)
        return 1;
    // retrieve the regions this variants located first, check the memory pool and update the pool if the regions is out of
    // cached positions. just skip if the variant type of line is a ref.    
    if ( line->pos < 0 )
	return 1;
    // skip reference positions, usually bcfanno will check the type of position in the first step, get error if see ref here
    if ( bcf_get_variant_types(line) == VCF_REF )
	error("This is a reference position. Should not come here. %s : %d", bcf_seqname(opts->hdr_out, line), line->pos+1);

    // fill memory pool and init hgvs names structure, return 1 if no transcript found
    if ( hgvs_names_init(opts, line) != 0 )
	return 1;
    // setter each column into vcf line
    int i;
    for ( i = 0; i < opts->n_cols; ++i ) {
	struct anno_col *col = &opts->cols[i];
	col->setter.hgvs(opts, line, col);
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
    } else if ( strcmp(key, "ExIn_id") == 0 ) {
        return generate_exin_ids(cache);
    } else if ( strcmp(key, "VarType") == 0 ) {
        return generate_vartype_string(cache);
    } else {
	error("Failed to recongnise key %s", key);
    }
}
int setter_hgvs_names(struct refgene_options *opts, bcf1_t *line, struct anno_col *col)
{
    char *string = hgvs_get_names_string(&opts->cache, col->hdr_key);

    if ( string == NULL )
      return 1;
    setter_hgvs_string(opts->hdr_out, line, col->hdr_key, string);
    if ( string )
      free(string);
    return 0;
}
int refgene_columns_parse(struct refgene_options *opts, char *column)
{
    if ( column == NULL)
	return 1;
    kstring_t string = KSTRING_INIT;
    kputs(column, &string);
    if (opts->hdr_out == NULL)
	error("Usage error. opts::hdr_out header structure is not inited yet.");
    int nfields = 0;
    int *splits = ksplit(&string, ',', &nfields);
    if (nfields == 0)
	return 1;
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
	    col->icol = hgvs_bcf_header_add_hgvsdna(opts->hdr_out);
        } else if ( strcmp(col->hdr_key, "ExIn_id") == 0 ) {
            col->icol = hgvs_bcf_header_add_exid(opts->hdr_out);
        } else if ( strcmp(col->hdr_key, "VarType") == 0 ) {
            col->icol = hgvs_bcf_header_add_vartype(opts->hdr_out);
	} else {
	    error("Failed to parse column tag %s", col->hdr_key);
	}
	col->number =  bcf_hdr_id2length(opts->hdr_out, BCF_HL_INFO, col->icol);
	col->setter.hgvs = setter_hgvs_names;
    }
    if ( opts->refgene_is_inited == 0 )
	opts->refgene_is_inited = 1;
    free(string.s);
    free(splits);
    return 0;
}
int refgene_set_refseq_fname(struct refgene_options *opts, const char *fname)
{
    const char **var = &opts->refseq_fname;
    *var = fname;
    opts->refseq_fai = load_refseq_file(fname);
    refseq_cache.idx = opts->refseq_fai;
    return 0;
}
int refgene_set_refgene_fname(struct refgene_options *opts, const char *fname)
{
    if ( opts->genepred_fname == NULL ) {
	const char **var = &opts->genepred_fname;
	*var = fname;
    }
    opts->fp = hts_open(opts->genepred_fname, "r");
    if (opts->fp == NULL)
	error("Failed to open %s", opts->genepred_fname);
        
    opts->genepred_idx = load_genepred_file(opts->genepred_fname);
    if (opts->genepred_idx == NULL)
	error("Failed to load index of %s", opts->genepred_fname);
    return 0;
}
int refgene_set_trans_fname(struct refgene_options *opts, const char *fname)
{
    const char **var = &opts->trans_list_fname;
    *var = fname;
    // todo: set transcripts list here
    return 0;
}
int refgene_set_genes_fname(struct refgene_options *opts, const char *fname)
{
    const char **var = &opts->genes_list_fname;
    *var = fname;
    // todo: set gene list here
    return 0;
}
int refgene_options_init(struct refgene_options *opts)
{
    memset(opts, 0, sizeof(struct refgene_options));    
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
static const char *hts_bcf_wmode(int file_type)
{
    if ( file_type == FT_BCF )
	return "wbu";    // uncompressed BCF
    if ( file_type & FT_BCF )
	return "wb";      // compressed BCF
    if ( file_type & FT_GZ )
	return "wz";       // compressed VCF
    return "w";                                 // uncompressed VCF
}

int main(int argc, char **argv)
{
    int i;
    const char *out_type_string = 0;
    const char *out_fname = 0;
    const char *input_fname = 0;
    const char *columns = "Gene,Transcript,HGVSDNA,ExIn_id,VarType";
    for ( i = 1; i< argc; ) {
	const char *a = argv[i++];
	if (strcmp(a, "-h") == 0) 
	    return usage(argc, argv);
	const char **arg_var = 0;
	if (strcmp(a, "-o") == 0 && out_fname == 0)
	    arg_var = &out_fname;
	else if (strcmp(a, "-O") == 0 && out_type_string == 0)
	    arg_var = &out_type_string;
	else if (strcmp(a, "-data") == 0 && opts.genepred_fname == 0)
	    arg_var = &opts.genepred_fname;
	else if (strcmp(a, "-refseq") == 0 && opts.refseq_fname == 0)
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
 
    if (opts.genepred_fname == 0)
	error("genepredPlus database is required.");
    if ( opts.refseq_fname == 0 )
        error("Refseq databases is required.");
    
    refgene_set_refgene_fname(&opts, opts.genepred_fname);
    refgene_set_refseq_fname(&opts, opts.refseq_fname);

    htsFile *fp = hts_open(input_fname, "r");
    if (fp == NULL)
	error("Failed to open %s.", input_fname);

    // check input type is VCF/BCF or not
    htsFormat type = *hts_get_format(fp);
    if (type.format  != vcf && type.format != bcf)
	error("Unsupported input format! %s", input_fname);

    int out_type = FT_VCF;
    if (out_type_string != 0) {
	switch (out_type_string[0]) {
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
    LOG_print("Parse header.");
    // duplicate header file and sync new tag in the output header .
    // assuming output is stdout in Vcf format in default. 
    htsFile *fout = out_fname == 0 ? hts_open("-", hts_bcf_wmode(out_type)) : hts_open(out_fname, hts_bcf_wmode(out_type));
    opts.hdr_out = bcf_hdr_dup(hdr);
    set_format_genepredPlus();
    // parse configure columns
    refgene_columns_parse(&opts, (char*)columns);   
    bcf_hdr_write(fout, opts.hdr_out);    
    // init gene_pred or refgene database, hold tabix index cache and faidx cache in memory 
    bcf1_t *line = bcf_init();
    while ( bcf_read(fp, hdr, line) == 0 ) {     
	anno_refgene_core(&opts, line);
	bcf_write1(fout, opts.hdr_out, line);
    }
    LOG_print("Clear ...");
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    bcf_hdr_destroy(opts.hdr_out);
    refgene_options_destroy(&opts);
    hts_close(fp);
    hts_close(fout);
    double c1 = get_time();
    fprintf(stderr, "Run time: %.2fs\n", c1 -c0);
    return 0;
}
#endif
