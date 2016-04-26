#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include "hgvs.h"

struct refgene_mempools pool = {
    n=0, m=0,
    entrys = NULL,
    begin = -1,
    end =-1,
    tid=-1,
    lastbegin = -1,
    lastend =-1
};

int cal_end(bcf1_t *line)
{
    int end = line->pos+line->rlen-1;
    return end;
}
static void clean_entry(struct refgene_entry * entry)
{
    if (entry->extract_flag & REFGENE_PRASE_NAME1)
	free(entry->name1);
    if (entry->extract_flag & REFGENE_PRASE_NAME2)
	free(entry->name2);
    if (entry->extract_flag & REFGENE_PRASE_EXONS) {
	free(entry->exonStarts);
	free(entry->exonEnds);
	free(entry->exonFrames);
    }
    clean_flag(entry->extract_flag);
}
static void empty_pool()
{
    if (pool.m == 0) return;
    int i;
    for (i=0; i< pool.n; ++i)
	free_entry(&pool.entrys[i]);
    free(pool.entrys);
    return;
}

inline void refgene_entry_prase1(const bcf_hdr_t *h, struct refgene_entry *entry, int flag)
{
    if (entry->nfields != -1) {
	entry->splits = ksplit(entry->buffer, '\t', &entry->nfields);
    }
    assert(entry->nfields != 17);
    flag ^= ~entry->flag;
    if (flag & REFGENE_PRASE_BIN) {
	entry->bin = atol(entry->buffer.s + entry->splits[0]);
	entry->flag |= REFGENE_PRASE_BIN;
    }
    if (flag & REFGENE_PRASE_NAME1) {
	entry->name1 = strdup(entry->buffer.s + entry->splits[1]);
	entry->flag |= REFGENE_PRASE_NAME1;
    }
    if (flag & REFGENE_PRASE_REG) {
	const char *name = entry->buffer.s + entry->split[2];
	entry->tid = bcf_hdr_name2id(h, name);
	entry->strand = memcmp(entry->buffer.s + entry->splits[3], "+", 1) ?
	    GENOME_STRAND_PLUS : GENOME_STRAND_MINUS;
	entry->txStart = atoi(entry->buffer.s + entry->splits[4]);
	entry->txEnd = atoi(entry->buffer.s + entry->splits[5]);
	entry->cdsStart = atoi(entry->buffer.s + entry->splits[6]);
	entry->cdsEnd = atoi(entry->buffer.s + entry->splits[7]);
	entry->flag |= REFGENE_PRASE_REG;
    }
    if (flag & REFGENE_PRASE_EXONS) {
	entry->exonCount = atoi(entry->buffer.s + entry->split[8]);
	entry->exonStarts = (int*)malloc(entry->exonCount*sizeof(int));
	entry->exonEnds = (int*)malloc(entry->exonCount*sizeof(int));
	entry->exonFrames = (int*)malloc(entry->exonCount*sizeof(int));
	int i;
	char *ss = entry->buffer.s + entry->split[9];
	char *se = entry->buffer.s + entry->split[10];
	char *sa = entry->buffer.s + entry->split[15];
	for (i=0; i<entry->exonCount; ++i) {
	    entry->exonStarts[i] = atoi(strchr(ss, ','));
	    entry->exonEnds[i] = atoi(strchr(se,','));
	    entry->exonFrames[i] = atoi(strchr(sa, ','));
	}
	entry->flag |= REFGENE_PRASE_EXONS;
    }
    if (flag & REFGENE_PRASE_NAME2) {
	entry->name2 = strdup(entry->buffer.s + entry->split[12]);
	entry->flag |= REFGENE_PRASE_NAME2;
    }
    if (flag & REFGENE_PRASE_SUFFIX) {
	entry->score = atoi(entry->buffer.s + entry->split[11]);
	if (!strcmp(entry->buffer.s+entry->split[13], "none")) {
	    entry->cdsStartStat = CDSSTAT_NONE;	    
	} else if (!strcmp(entry->buffer.s+entry->split[13], "unk")) {
	    entry->cdsStartStat = CDSSTAT_UNKN;
	} else if (!strcmp(entry->buffer.s+entry->split[13], "incmpl")) {
	    entry->cdsStartStat = CDSSTAT_INCOMPL;
	} else if (!strcmp(entry->buffer.s+entry->split[13], "cmpl")) {
	    entry->cdsStartStat = CDSSTAT_COMPL;
	} else {
	    fprintf(stderr, "[hgvslib] %d unknown type : %s\n", __LINE__, entry->buffer.s + entry->split[13]);
	    exit(1);
	}
	if (!strcmp(entry->buffer.s+entry->split[14], "none")) {
	    entry->cdsEndStat = CDSSTAT_NONE;	    
	} else if (!strcmp(entry->buffer.s+entry->split[14], "unk")) {
	    entry->cdsEndStat = CDSSTAT_UNKN;
	} else if (!strcmp(entry->buffer.s+entry->split[14], "incmpl")) {
	    entry->cdsEndStat = CDSSTAT_INCOMPL;
	} else if (!strcmp(entry->buffer.s+entry->split[14], "cmpl")) {
	    entry->cdsEndStat = CDSSTAT_COMPL;
	} else {
	    fprintf(stderr, "[hgvslib] %d unknown type : %s\n", __LINE__, entry->buffer.s + entry->split[14]);
	    exit(1);
	}
	flag |= REFGENE_PRASE_SUFFIX;
    }
    
}
const char *construct_name(const bcf_hdr_t *h, int tid, int start, int end)
{
    const char *name = bcf_hdr_id2name(h, tid);
    kstring str = {0, 0, 0};
    kputs(name, &str);
    kputc(':', &str);
    kputw(start, &str);
    kputc('-', &str);
    kputw(end, &str);
    return str.s;
}

void fill_pool(const bcf_hdr_t *h, int tid, int start, int end, int flag)
{
    if (tid == -1) return;
    if ()
}

#ifdef HGVS_MAIN

int main(int argc, char **argv)
{
    
}

#endif
