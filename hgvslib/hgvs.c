#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include <htslib/khash.h>
#include <htslib/khash_str2int.h>
#include <htslib/kstring.h>
#include "hgvs.h"

// keep gene list and transcript list into two hash structures 
KHASH_MAP_INIT_STR(list, char*)


// files open and keep in mempool
struct refgene_mempools {
    htsFile *fp;
    const char *fn;
    const tbx_t *idx;
    const char *fa; // transcript fasta file
    const faidx_t *fai;
    /* todo: badaln */
    int tid;
    int n, m;
    struct refgene_entry *entries;
    int begin;
    int end;
    int lastbegin;
    int lastend;
    int has_genelist;
    int has_translist;
    khash_t(list) *ghash;
    khash_t(list) *thash;
} pool = {
    .fp=NULL,
    .fn=NULL,
    .idx=NULL,
    .n=0,
    .m=0,
    .entries = NULL,
    .begin = -1,
    .end =-1,
    .tid=-1,
    .lastbegin = -1,
    .lastend =-1,
    .has_genelist = 0,
    .has_translist = 0,
};

static void empty_mempool();

int generate_gene_hash(char *fn)
{
    int i, n=0;
    char **names = hts_readlist(fn, 1, &n);
    if (n == 0)
	return 0;
    pool.ghash = kh_init(list);
    int ignore, k;
    for (i=0; i<n; ++i) {
	k = kh_put(list, ghash, names[i], &ignore);
    }
    return n;
}
int generate_trans_hash(char *fn)
{
    int i, n=0;
    char **names = hts_readlist(fn, 1, &n);
    if (n == 0)
	return 0;
    pool.thash = kh_init(list);
    int ignore, k;
    for (i=0; i<n; ++i) {
	k = kh_put(list, thash, names[i], &ignore);
    }
    return n;
}

int generate_gene_hash(char *fn)
{
    int i, n=0;
    char **names = hts_readlist(fn, 1, &n);
    pool.ghash = kh_init(list);
    int ignore, k;
    for (i=0; i<n; ++i) {
	k = kh_put(list, ghash, names[i], &ignore);
    }
}

void destroy_mempool()
{
    empty_mempool();
    if (pool.fp) hts_close(pool.fp);
    if (pool.idx) tbx_destroy(pool.idx);
    if (pool.m) free(pool.entries);
}
int cal_end(bcf1_t *line)
{
    int end = line->pos+line->rlen;
    return end;
}
static void init_entry(struct refgene_entry *entry)
{
    entry->nfields = -1;
    entry->splits = 0;
    entry->name1 = entry->name2 = NULL;
    entry->flag = 0;
    entry->bin = -1;
    entry->tid = -1;
    entry->cdsStart = entry->cdsEnd = -1;
    entry->exonStarts = entry->exonEnds = entry->exonFrames = NULL;
    entry->score = -1;    
    entry->txStart = entry->txEnd = -1;
    entry->flag = 0;
    entry->buffer.l = entry->buffer.m = 0;
    entry->buffer.s = 0;
}
static void clean_entry1(struct refgene_entry *entry)
{
    assert(entry);
    if (entry->empty == 1) return;
    
    if (entry->flag & REFGENE_PRASE_EXONS) {
	free(entry->exonStarts);
	free(entry->exonEnds);
	free(entry->exonFrames);
    }
    if (entry->buffer.m > 0) {
	free(entry->buffer.s);
	entry->buffer.l = entry->buffer.m = 0;
	entry->buffer.s = 0;
    }
    if (entry->nfields) {
	free(entry->splits);
    }
    entry->flag = 0;
    entry->empty = 1;
}

static void empty_mempool()
{
    if (pool.n == 0) return;
    int i;
    for (i=0; i< pool.n; ++i) {
	clean_entry1(&pool.entries[i]);
    }
    pool.n=0;
    
    pool.tid = -1;
    pool.lastbegin = pool.begin;
    pool.lastend = pool.end;
    pool.begin = -1;
    pool.end = -1;
    return;
}

void refgene_entry_prase1(const bcf_hdr_t *h, struct refgene_entry *entry, int flag)
{
    assert(entry);
    assert(entry->empty == 0);
    // init buffer
    if (entry->nfields == -1) {
	entry->splits = ksplit(&entry->buffer, 0, &entry->nfields);
    }
    
    assert(entry->nfields == 16);
    
    flag &= ~entry->flag;

    if (flag & REFGENE_PRASE_BIN) {
	entry->bin = atoi(entry->buffer.s + entry->splits[0]);
	entry->flag |= REFGENE_PRASE_BIN;
    }
    
    if (flag & REFGENE_PRASE_NAME1) {
	assert(entry->name1 == NULL);
	//entry->name1 = strdup(entry->buffer.s + entry->splits[1]);
	entry->name1 = entry->buffer.s + entry->splits[1];
//	fprintf(stderr, "name1: %s\n", entry->name1);
	entry->flag |= REFGENE_PRASE_NAME1;
    }
    
    if (flag & REFGENE_PRASE_REG) {
	const char *name = entry->buffer.s + entry->splits[2];
	entry->tid = bcf_hdr_name2id(h, name);
	entry->strand = memcmp(entry->buffer.s + entry->splits[3], "+", 1) ?
	    GENOME_STRAND_PLUS : GENOME_STRAND_MINUS;	    
	entry->txStart = atoi(entry->buffer.s + entry->splits[4]);
	entry->txEnd = atoi(entry->buffer.s + entry->splits[5]);
	entry->cdsStart = atoi(entry->buffer.s + entry->splits[6]);
	entry->cdsEnd = atoi(entry->buffer.s + entry->splits[7]);
	entry->flag |= REFGENE_PRASE_REG;
//	fprintf(stderr, "txstart: %d\ttxend: %d\tcdstart: %d\tcdsend: %d\n", entry->txStart, entry->txEnd, entry->cdsStart, entry->cdsEnd);
    }
    
    if (flag & REFGENE_PRASE_EXONS) {
	entry->exonCount = atoi(entry->buffer.s + entry->splits[8]);
	entry->exonStarts = (int*)malloc(entry->exonCount*sizeof(int));
	entry->exonEnds = (int*)malloc(entry->exonCount*sizeof(int));
	entry->exonFrames = (int*)malloc(entry->exonCount*sizeof(int));
	int i;
	char *ss = entry->buffer.s + entry->splits[9];
	char *se = entry->buffer.s + entry->splits[10];
	char *sa = entry->buffer.s + entry->splits[15];
	for (i=0; i<entry->exonCount; ++i) {
	    entry->exonStarts[i] = atoi(strchr(ss, ','));
	    entry->exonEnds[i] = atoi(strchr(se,','));
	    entry->exonFrames[i] = atoi(strchr(sa, ','));
	}
	entry->flag |= REFGENE_PRASE_EXONS;
    }
    
    if (flag & REFGENE_PRASE_NAME2) {
	//entry->name2 = strdup(entry->buffer.s + entry->splits[12]);
	entry->name2 = entry->buffer.s + entry->splits[12];
	entry->flag |= REFGENE_PRASE_NAME2;
//	fprintf(stderr, "name2: %s\n", entry->name2);
    }
    
    if (flag & REFGENE_PRASE_SUFFIX) {
	entry->score = atoi(entry->buffer.s + entry->splits[11]);
	if (!strcmp(entry->buffer.s+entry->splits[13], "none")) {
	    entry->cdsStartStat = CDSSTAT_NONE;	    
	} else if (!strcmp(entry->buffer.s+entry->splits[13], "unk")) {
	    entry->cdsStartStat = CDSSTAT_UNKN;
	} else if (!strcmp(entry->buffer.s+entry->splits[13], "incmpl")) {
	    entry->cdsStartStat = CDSSTAT_INCMPL;
	} else if (!strcmp(entry->buffer.s+entry->splits[13], "cmpl")) {
	    entry->cdsStartStat = CDSSTAT_CMPL;
	} else {
	    fprintf(stderr, "[hgvslib] %d unknown type : %s\n", __LINE__, entry->buffer.s + entry->splits[13]);
	    exit(1);
	}
	if (!strcmp(entry->buffer.s+entry->splits[14], "none")) {
	    entry->cdsEndStat = CDSSTAT_NONE;	    
	} else if (!strcmp(entry->buffer.s+entry->splits[14], "unk")) {
	    entry->cdsEndStat = CDSSTAT_UNKN;
	} else if (!strcmp(entry->buffer.s+entry->splits[14], "incmpl")) {
	    entry->cdsEndStat = CDSSTAT_INCMPL;
	} else if (!strcmp(entry->buffer.s+entry->splits[14], "cmpl")) {
	    entry->cdsEndStat = CDSSTAT_CMPL;
	} else {
	    fprintf(stderr, "[hgvslib] %d unknown type : %s\n", __LINE__, entry->buffer.s + entry->splits[14]);
	    exit(1);
	}
	flag |= REFGENE_PRASE_SUFFIX;
    }
}

const char *construct_name(const bcf_hdr_t *h, int tid, int start, int end)
{
    const char *name = bcf_hdr_id2name(h, tid);
    kstring_t str = {0, 0, 0};
    kputs(name, &str);
    kputc(':', &str);
    kputw(start, &str);
    kputc('-', &str);
    kputw(end, &str);
    return str.s;
}

void push_mempool(const bcf_hdr_t *h, kstring_t *str, int flag)
{

    if (pool.n == pool.m) {
	pool.m += 10;
	pool.entries = (struct refgene_entry *)realloc(pool.entries, pool.m*sizeof(struct refgene_entry));
    }
    struct refgene_entry *entry = &pool.entries[pool.n];
    pool.n++;
    init_entry(entry);
    entry->buffer.s = strndup(str->s, str->l);
    entry->buffer.m = str->l;
    entry->buffer.l = str->l;
    entry->empty = 0;
    refgene_entry_prase1(h, entry, flag);
    //entry->empty = 0;
}
void update_mempool(const bcf_hdr_t *h)
{
    if (pool.n == 0) return;
    pool.begin = pool.end = -1;
    pool.tid = -1;
    int i;
    for (i = 0; i <pool.n; ++i) {
	struct refgene_entry *entry = &pool.entries[i];
	if (entry->empty == 1) continue;
	if (!(entry->flag & REFGENE_PRASE_REG))
	    refgene_entry_prase1(h, entry, REFGENE_PRASE_REG);
	if (pool.begin == -1 || pool.begin > entry->txStart)
	    pool.begin = entry->txStart;
	if (pool.end < entry->txEnd)
	    pool.end = entry->txEnd;
    }
    pool.tid = pool.entries[i].tid;
}
void filter_mempool(const bcf_hdr_t *h)
{
    khash_t(list) *h1, *h2;
    int i;
    h1 = pool.ghash;
    h2 = pool.thash;
    
    for (i=0; i<pool.n; ++i) {
	struct refgene_entry *entry = &pool.entries[i];
	if (entry->empty == 1) continue;
	int del_it = 0;
	khiter_t iter1, iter2;
	if (pool.has_genelist == 1) {
	    del_it = 1;
	    if (!(entry->flag & REFGENE_PRASE_NAME1))
		refgene_entry_prase1(h, entry, REFGENE_PRASE_NAME1);
	    iter1 = kh_get(list, h1, entry->name1);
	    if (iter1 != kh_end(h1))
		del_it = 0;
	}
	
	if (pool.has_translist == 1) {
	    del_it = 1;
	    if (!(entry->flag & REFGENE_PRASE_NAME2))
		refgene_entry_prase1(h, entry, REFGENE_PRASE_NAME2);
	    iter2 = kh_get(list, h2, entry->name2);
	    if (iter2 != kh_end(h1))
		del_it = 0;
	}
	if (del_it == 1)
	    clean_entry1(entry);
    }
}

void fill_mempool(const bcf_hdr_t *h, int tid, int start, int end, int flag)
{
    if (tid == -1) return;
    assert (pool.fn != NULL);
    if (pool.fp == NULL)
	pool.fp = hts_open(pool.fn, "r");
    assert (pool.fp != NULL);
    if (pool.idx == NULL) {
	pool.idx = tbx_index_load(pool.fn);
	if ( !pool.idx ) {
	    fprintf(stderr, "Could not load .tbi/.csi index of %s\n", pool.fn);
	    exit(1);
	}
    }
    assert (pool.idx != NULL);
    update_mempool(h);

    if (pool.tid != -1 && pool.end != -1 &&
	tid == pool.tid && start <= pool.end+UTR3_REG) {
	return;
    } else {
	empty_mempool();
    }
    
    const char *reg = construct_name(h, tid, start, end);
    hts_itr_t *itr = tbx_itr_querys(pool.idx, reg);
    free((void*)reg);
    if ( !itr) return;
    
    kstring_t str = {0, 0, 0};
    while (tbx_itr_next(pool.fp, pool.idx, itr, &str) >=0) {
	push_mempool(h, &str, flag);
	str.l = 0;
    }
    free(str.s);
    tbx_itr_destroy(itr);
    update_mempool(h);
}


struct hgvs_record *(bcf_hdr_t *h, bcf1_t *line)
{
    int i,n=0;
    int end = cal_end(line);
    int flag = REFGENE_PRASE_REG | REFGENE_PRASE_NAME1 | REFGENE_PRASE_NAME2;
    fill_mempool(h, line->rid, line->pos, end, flag);
    filter_mempool(h); // filter records if set transcript or gene list

    struct hgvs_record *rec;
    rec = (struct hgvs_record *)calloc(1, sizeof(struct hgvs_record));
    rec->ntrans = rec->mtrans = 0;
    
    int x = bcf_get_variant_type(line);

    // skip ref record in the first place, it is rudentency, remove this part later
    if (x == VCF_REF|| line->n_allele < 2) {
	free(rec);
	return NULL;
    }

    // loop the mempool
    for (i=0; i<pool.n; ++i) {

	if (pool.entries[i].empty == 1)
	    continue;

	struct refgene_entry *entr = &pool.entries[i];
	int start = entr->txStart;
	int stop = entr->txEnd;

	if (line->pos >= start && line->pos <= stop) {
	    if (rec->ntrans == rec->mtrans) {
		rec->mtrans = rec->mtrans == 0 ? 2 : rec->mtrans + 2;
		rec->trans = (struct hgvs1*)realloc(rec->trans, rec->mtrans*sizeof(struct hgvs1));
	    }
	    struct hgvs1 * htr = rec->trans[rec->ntrans++];
	    htr->empty = 1;
	    htrr->n_allele = line->n_allele;
	    // trans convertor name come here
	    // check coding or noncoding transcript
	    assert(htr->name1);
	    if (entr->name1[1] == 'M') htr->type = TYPE_CODING;
	    else htr->type = TYPE_NONCODING;
	    // calculate position
	    int j;
	    int loc = 0; // locations in exon or cds
	    int off = 0;  // offset in noncoding region

/* 
   UTR_REG5    .  UTR_REG3   CODING     UTR_REG5.UTR_REG3
   ============ATG=================================
*/
		
	    // tx region
	    int pos = line->pos; // variant position

	    if (pos < )
	    if (entr->type != TYPE_CODING) {
		// just count tx position
		htr->type = TYPE_NONCODING;
		if ()
		for (j=0; j<entr->exonCount; ++j) {
		    if (line->pos > entr->exonEnds[j]) {
			pos += entr->exonEnds[j] - entr->exonStarts[j] + 1;
		    } else {
			if (pos > entr->exonStarts[j]) {
			    int l = line->pos - entr->exonStarts[j] + 1;
				pos += l;
				if (l <= UTR3_REG) {
				    htr->func = FUNC_SPLITSITE5; 
				} else {
				    htr->func = FUNC_UNKNOWN;
				}
			    } else if (pos == entr->exonStarts[j]) {
				pos +=1;
				htr->func = FUNC_SPLITSITE;
			    } else {
				//pos < entr->exonStarts[j]
				// TODO: check downstream or upstream in a intron
				offset = line->pos - entr->exonStarts[j];
				if (j==0) htr->func = FUNC_UPSTREAM;
				else htr->func = FUNC_INTRON;
			    }
			}
		    } // exon loop
		    htr->exId = j;
		    htr->cdsId = -1;
		    if (j == entr->exonCount && line->pos > entr->exonEnds[j-1]) {
			offset = line->pos - entr->exonEnds[j-1];
			htr->func = FUNC_DOWNSTREAM;
		    }
			
		} else {
		    htr->type = TYPE_CODING;
		    int cc = 0;
		    int n=0;
		    /* if (line->pos < entr->cdsStart) { */
		    /* 	htr->func = FUNC_UPSTREAM; */
		    /* 	if (line->pos >= entr->txStart) { */
		    /* 	    htr->func = FUNC_5UTR; */
		    /* 	    for (j=0; j<entr->exonCount; ++j) */
		    /* 		if (line->pos */
		    /* 	} */
		    /* 	if (line->pos >= entr->cdsStart - UTR5_REG) { */
		    /* 	    htr->func = FUNC_SPLITSITE5; */
		    /* 	} */
		    /* 	pos = 0; */
		    /* 	offset = line->pos - entr->cdsStart; */
		    /* 	htr->cdsId = 0; // cdsId == 0 stand for UTR5 or Intergenic */
		    /* }  */

		    /*
		        *       *         *    *       *           *
		                            =======  ============
		           -----------   ----------  ------------
		     */


		    if (line->pos < entr->txStart) {
			if (line->pos >= entr->txStart + UTR5_REG) {
			    htr->func = FUNC_SPLITSITE5;
			    htr
		    }
		    for (j=0; j<entr->exonCount; ++j) {
			if (line->pos >= entr->exonStarts[j] && line->pos <= entr->exonEnds[j]) {
			    htr->exId = j;
			    
			    
			    
			}
			if (entr->exonEnds[j] > entr->cdsStart) {
			    // cds region come here
			    cc++;
			    if (entr->exonStarts[j] <= entr->cdsStart) {
				pos += exonEnds[j] - entr->cdsStart + 1;
			    } else {
				pos += exonEnds[j] - entr->exonStarts[j] + 1;
			    }
			    
			} else {
			    continue;
			}
		    }
		    
		}
	    }
	    // function predition
	    
	    // generate name for each allele
	    entr->als = (struct hgvs_ale*)malloc(entr->n_allele*sizeof(struct hgvs_ale));

	    for (j=1; j< entr->n_allele; ++j) {
		
	    }
	    
	    
	}
    }


    
}

#ifdef HGVS_MAIN
#include <getopt.h>
#include <htslib/synced_bcf_reader.h>
int usage() {
    fprintf(stderr, "hgvs refgene.txt.gz in.vcf\n");
    return 1;
}
int main(int argc, char **argv)
{
    if (argc - optind != 2)
	return usage();
    bcf_srs_t *file = bcf_sr_init();

    if (!bcf_sr_add_reader(file, argv[2])) {
	fprintf(stderr, "Failed to open %s : %s\n", argv[2], bcf_sr_strerror(file->errnum));
	exit(1);
    }
    pool.fn = argv[1];

    int flag = REFGENE_PRASE_REG | REFGENE_PRASE_NAME1 | REFGENE_PRASE_NAME2;

    while (bcf_sr_next_line(file)) {
	bcf1_t *line = file->readers[0].buffer[0];
	int end = cal_end(line);
	fill_mempool(file->readers[0].header, line->rid, line->pos, end, flag);
    }

    bcf_sr_destroy(file);
    destroy_mempool();
    return 0;
}

#endif
