#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>
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

struct hgvs *pos_anno_hgvs(bcf1_t *line)
{
    int end = line->pos+line->rlen-1;
    
}

void fill_pool (int tid, int start, int end)
{
    
}
