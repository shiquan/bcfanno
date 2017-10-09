// Add the FLANKSEQ tag for each variant in the INFO
#include "utils.h"
#include "htslib/faidx.h"
#include "htslib/vcf.h"

// export flank sequence arount target variant
static int flank_size = 10;

void set_flksize(int size)
{
    assert(flank_size > 0);
    flank_size = size;
}

struct seqidx {
    const char *file;
    faidx_t *idx;
} idx = {
    .file = NULL,
    .idx = NULL,
};

int load_sequnce_index(const char *file)
{
    idx.file = file;
    idx.idx = fai_load(file);
    if ( idx.idx == NULL )
        return 1;
    return 0;
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
int bcf_add_flankseq(bcf_hdr_t *hdr, bcf1_t *line)
{
    const char *name = bcf_hdr_id2name(hdr, line->rid);
    int end = line->pos + line->rlen + flank_size;
    int start = line->pos + 1 - flank_size;
    int l_seq = 0;
    char *seq = faidx_fetch_seq(idx.idx, name, start-1, end-1, &l_seq);
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

void seqidx_destroy()
{
    fai_destroy(idx.idx);
}
