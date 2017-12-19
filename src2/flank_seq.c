/*  
    Copyright (C) 2016,2017  BGI Research

    Author: Shi Quan (shiquan@genomics.cn)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE. 
*/

// Add the FLANKSEQ tag for each variant in the INFO
#include "utils.h"
#include "htslib/faidx.h"
#include "htslib/vcf.h"
#include "anno_flank.h"

// export flank sequence arount target variant
static int flank_size = 10;

void set_flksize(int size)
{
    assert(flank_size > 0);
    flank_size = size;
}

struct seqidx* load_sequence_index(const char *file)
{
    struct seqidx *idx = malloc(sizeof(*idx));
    idx->file = file;
    idx->idx = fai_load(file);
    if ( idx->idx == NULL ) {
        warnings("failed to load index of %s.", file);
        free(idx);
        return NULL;
    }
    return idx;
}
struct seqidx *sequence_index_duplicate(struct seqidx *idx)
{
    struct seqidx *d = malloc(sizeof(*d));
    d->file = idx->file;
    d->idx = fai_load(d->file);
    assert(d->idx);
    return d;
}

void sequence_index_destroy(struct seqidx *idx)
{
    fai_destroy(idx->idx);
    free(idx);
}

int bcf_header_add_flankseq(bcf_hdr_t *hdr)
{
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, "FLKSEQ");
    if (id == -1) {
	bcf_hdr_append(hdr, "##INFO=<ID=FLKSEQ,Number=1,Type=String,Description=\"Upstream and downstream bases of current position on DNA sequence.\">");
	bcf_hdr_sync(hdr);
	id = bcf_hdr_id2int(hdr, BCF_DT_ID, "FLKSEQ");
	assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));
    }
    return id;
}
int bcf_add_flankseq(struct seqidx *idx, bcf_hdr_t *hdr, bcf1_t *line)
{
    const char *name = bcf_hdr_id2name(hdr, line->rid);
    int end = line->pos + line->rlen + flank_size;
    int start = line->pos + 1 - flank_size;
    int l_seq = 0;
    char *seq = faidx_fetch_seq(idx->idx, name, start-1, end-1, &l_seq);
    if ( seq == NULL || end - start + 1 != l_seq ) {
        if ( seq ) free(seq);
        return 1;
    }
    kstring_t str = { 0, 0, 0,};
    kputsn(seq, flank_size, &str);
    kputc('.', &str);
    kputsn(seq + (l_seq - flank_size), flank_size, &str);
    bcf_update_info_string(hdr, line, "FLKSEQ", str.s);
    free(seq);
    free(str.s);
    return 0;
}
