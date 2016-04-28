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

struct refgene_mempools pool = {
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
    .lastend =-1
};

int cal_end(bcf1_t *line)
{
    int end = line->pos+line->rlen;
    return end;
}
static void clean_entry(struct refgene_entry * entry)
{
    if (entry->flag & REFGENE_PRASE_NAME1)
	free(entry->name1);
    if (entry->flag & REFGENE_PRASE_NAME2)
	free(entry->name2);
    if (entry->flag & REFGENE_PRASE_EXONS) {
	free(entry->exonStarts);
	free(entry->exonEnds);
	free(entry->exonFrames);
    }
    entry->txStart = entry->txEnd = -1;
    entry->flag = 0;
}
static void empty_mempool()
{
    if (pool.m == 0) return;
    int i;
    for (i=0; i< pool.n; ++i)
	clean_entry(&pool.entries[i]);
    free(pool.entries);
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
    if (entry->nfields != -1) {
	entry->splits = ksplit(&entry->buffer, '\t', &entry->nfields);
    }
    assert(entry->nfields == 16);

    flag &= ~entry->flag;
    if (flag & REFGENE_PRASE_BIN) {
	entry->bin = atol(entry->buffer.s + entry->splits[0]);
	entry->flag |= REFGENE_PRASE_BIN;
    }
    if (flag & REFGENE_PRASE_NAME1) {
	entry->name1 = strdup(entry->buffer.s + entry->splits[1]);
	fprintf(stderr, "name1: %s\n", entry->name1);
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
	fprintf(stderr, "txstart: %d\ttxend: %d\tcdstart: %d\tcdsend: %d\n", entry->txStart, entry->txEnd, entry->cdsStart, entry->cdsEnd);
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
	entry->name2 = strdup(entry->buffer.s + entry->splits[12]);
	entry->flag |= REFGENE_PRASE_NAME2;
	fprintf(stderr, "name2: %s\n", entry->name2);
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
    entry->buffer.s = str->s;
    entry->buffer.m = str->m;
    entry->buffer.l = str->l;
    refgene_entry_prase1(h, entry, flag);    
}
void update_mempool(const bcf_hdr_t *h)
{
    int i;
    for (i = 0; i <pool.n; ++i) {
	struct refgene_entry *entry = &pool.entries[i];
	if (entry->flag ^ REFGENE_PRASE_REG)
	    refgene_entry_prase1(h, entry, REFGENE_PRASE_REG);
	if (pool.begin == -1 || pool.begin > entry->txStart)
	    pool.begin = entry->txStart;
	if (pool.end < entry->txEnd)
	    pool.end = entry->txEnd;
    }
}
void fill_pool(const bcf_hdr_t *h, int tid, int start, int end, int flag)
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
    if (tid == pool.tid && start <= pool.end+UTR3_REG) 
	return;
    else
	empty_mempool();
    const char *reg = construct_name(h, tid, start, end);
    hts_itr_t *itr = tbx_itr_querys(pool.idx, reg);
    if ( !itr) return;
    kstring_t str = {0, 0, 0};
    while (tbx_itr_next(pool.fp, pool.idx, itr, &str) >=0) 
	push_mempool(h, &str, flag);
    update_mempool(h);
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
	fprintf(stderr, "%d\t%d\t%d\t%d\n", line->pos+1, end, pool.begin, pool.end);
	fill_pool(file->readers[0].header, line->rid, line->pos, end, flag);
    }
    bcf_sr_destroy(file);
    empty_mempool();
    return 0;
}

#endif
