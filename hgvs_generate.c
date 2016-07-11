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


KHASH_MAP_INIT_STR(list, char*)
typedef khash_t(list) glist_t;

/* for each transcript, the region span several exon partial regions, use exon_pair[] stand each exons
 * and exon_offset_pair[] is the offset / coding coordiante consider of cdsstart 
 */
struct exon_pair {
    uint32_t start;
    uint32_t end;
};

/*
 * The dna reference offset for DNA coding and noncoding transcript consist of two parts. the offset value
 * and the function region construct a 32-bits value.
 *
 *   coding/nocoding coordinate
 *   32                        4321
 *   |_________________________||||
 *
 * ONLY one-bit below is accept per value   
 * offset 1:  is_noncoding
 * offset 2:  is_coding
 * offset 3:  is_utr5
 * offset 4:  is_utr3
 *  
 */
struct exon_offset_pair {
    uint32_t start;
    uint32_t end;
};

struct gene_predictions {
    char *id;
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

struct gene_predictions_mempool {
    int rid;
    uint32_t start;
    uint32_t end;
    int l, m;
    struct gene_predictions **a;
};

glist_t *init_gene_list(const char *list)
{
    int i, n = 0;
    char **names = hts_readlist(list, 1, &n);
    if (n== 0) return NULL;
    glist_t *glist = kh_init(list);
    int ig, k;
    for (i=0; i<n; ++n) {
	k = kh_put(list, glist, names[i]; &ig);
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
struct gene_predictions *genepred_prase_line(kstring_t *str)
{
    struct gene_predictions * gp = (struct gene_predictions*)malloc(sizeof(struct gene_predictions));
    int nfields = 0;
    int *splits = ksplit(str, 0, &nfields);

    /* accept any genepred-like format, like refGene, ensGene, gene_predictions*/
    assert(nfields == 16);
    gp->id = strdup(str->s + splits[2]);
    gp->name1 = atoi(str->s + splits[1]);
    gp->strand = memcmp(str->s + splits[3], "+", 1) ? '-' : '+';
    gp->txstart = atoi(str->s + splits[4]);
    gp->txend = atoi(str->s + splits[5]);
    gp->cdsstart = atoi(str->s + splits[6]);
    gp->cdsend = atoi(str->s + splits[7]);
    gp->exoncount = atoi(str->s + splits[8]);
    gp->name1 = atoi(str->s + splits[12]);
    /* the exons region in genepred file is like  1,2,3,4,5. try to init the exon_pair[] and exon_offset_pair[] by exoncount */
    char *ss = str->s + splits[9];
    char *se = str->s + splits[10];
    gp->exons = (struct exon_pair*)calloc(gp->exoncount, sizeof(struct exon_pair));
    gp->dna_ref_offsets = (struct exon_offset_pair*)calloc(gp->exoncount, sizeof(struct exon_offset_pair));

    int i;
    char *ss1, *se1;
    int dna_refernece_length = 0;
    int dna_coding_reference_length = 0;
    int dna_utr5_reference_length = 0;
    int dna_utr3_reference_length = 0;
    int is_coding_reference = gp->cdsstart == gp->cdsend ? 0 : 1;
    
    for (i=0; i<gp->exoncount; ++i) {
	ss1 = ss;
	while(*ss1 && *ss1 != ',') ss1++;
	ss1[0] = '\0';
	gp->exons[i].start = atoi(ss);
	ss = ss1++;
	se1 = se;
	while(*se && *se1 != ',') se1++;
	se1[0] = '\0';
	gp->exons[i].end = atoi(se);
	se = se1++;
	dna_reference_length += gp->exons[i].end - gp->exons[i].start + 1;
	if (is_coding_reference == 0) continue;
	if (gp->exons[i].end < gp->cdsstart) {
	    dna_utr5_reference_length += gp->exons[i].end - gp->exons[i].start +1;
	} else {
	    /* first cds */
	    if (gp->cdsstart > gp->exons[i].start) {
		dna_utr5_reference_length += gp->cdsstart - gp->exons[i].start;
	    } else {
		if (gp->cdsend < gp->exons[i].start) {
		    dna_utr3_reference_length += gp->cdsend - gp->cdsstart + 1;
		} else {
		    
		}
	    }
	}		

	    
    }
    /*  Coding DNA reference:
     *                   -3                                            *3
     *                    -2                                          *2
     * -93   -45 -44       -1 1                187 188           351 *1   *96 *97    *223
     *  |_______|____________|====================|=================|________|_______|
     *         / \            ATG                / \             TGA        / \
     *        /   \                             /   \                      /   \
     *       /gtga,,,g                         /gta,,,act                 /gtag,,,cc
     *   -45+1        -44-1                187+1        188-1         *96+1        *97-1
     *    -45+2      -44-2                  187+2      188-2           *96+2      *97-2
     *     -45+3    -44-3                    187+3    188-3             *96+3    *97-3
     *
     *
     * Noncoding DNA reference:
     *
     *  1        49 50             280 281                 540 541             667
     *  |__________|__________________|_______________________|_________________|
     * 
     */
    for (){
    }

}

void push_mempool(struct gene_predictions_mempool *pool, kstring_t *str)
{
    if (pool->m == pool->l) {
	pool->m = pool->m == 0 ? 2 : pool->m << 1;
	pool->a = (struct gene_predictions *)realloc(pool->a, sizeof(struct gene_predictions)*pool->m);
    }
    /* the prase func must return a point or go abort */
    pool->a[pool->l++] = gene_predictions_prase_line(str);
}
void update_mempool(struct gene_predictions_mempool *pool)
{    
    if (pool->l == 0) return;
    /* try to loop this pool and find the edges*/
    pool->start = pool->a[0]->txStart;
    pool->end = pool->a[0]->txEnd;
    int i;    
    for (i=1; i<pool->l; ++i) {
	//if (pool->a[i] == NULL) continue;
	if (pool->a[i]->txStart < pool->start) pool->start = pool->a[i]->txStart;
	if (pool->a[i]->txEnd > pool->end) pool->end = pool->a[i]->txEnd;
    }
}
void fill_mempool(struct gene_predictions_mempool *pool, htsFile *fp, tbx_t *tbx, int tid, uint32_t start, uint32_t end)
{
    hts_itr_t *itr = tbx_itr_queryi(tbx, tid, start, end);
    kstring_t str = KSTRING_INIT;
    while (tbx_itr_next(fp, idx, itr, &str) >= 0) {
	push_mempool(pool, &str);
	str.l = 0;
    }
    if (str.m) free(str.s);
    tbx_itr_destroy(itr);
    update_mempool(pool);    
}
void anno_hgvs_core(bcf_hdr_t *hdr, bcf1_t *line)
{
    
}



#include <sys/time.h>
static double get_time() {
    struct timeval tv;
    if ( gettimeofday(&tv, 0) != 0 ) 
	error("Failed to get time of day.");
    return (double)tv.tv_sec + 1.0e-6 * (double)tv.tv_usec;
}

struct args {
    char *gene_predictions_file;
    tbx_t *gpidx;
    faidx_t *faidx;
    glist_t *glist;
    struct gene_predictions_mempool pool;
};

int usage(int argc, char **argv)
{
    fprintf(stderr, "Usage: %s [-h] -refseq refseq.fa.gz -data gene_predictions.txt.gz|refgene.txt.gz -O [u|v|z|b] -o output_file in_file.vcf.gz \n", argv[0]);
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
    const char *in_file = 0; // stdin in default
    const char *out_file = 0; // stdout in default
    const char *out_type = 0; // output format, u in default
    const char *gp_file = 0; // error abort if not set gene_predictions file

    for (int i = 1; i< argc;) {
	const char *a = argv[i++];
	if (strcmp(a, "-h") == 0) 
	    return usage(argc, argv);

	const char **arg_var = 0;
	if (strcmp(a, "-o") == 0 && out_file == 0)
	    arg_var = &out_file;
	else if (strcmp(a, "-O") == 0 && out_type == 0)
	    arg_var = &out_type;
	else if (strcmp(a, "-data") == 0 && gp_file == 0)
	    arg_var = &gp_file;
	
	if (arg_var != 0) {
	    if (i == argc)
		error("Missing arg after %s", a);
	    *arg_var = argv[i++];
	    continue;
	}

	if (in_file == 0) {
	    in_file = a;
	    continue;
	}
	error("Unknown arg : %s", a);
    }
    /* assuming input file is stdin, only vcf, gzipped vcf, bcf file is accept,
       err msg will output if file type unrecongnized*/
    if (in_file == 0 && (!isatty(fileno(stdin)))) in_file = "-";

    htsFile *fp = hts_open(in_file, "r");
    if (fp == NULL)
	error("Failed to open %s.", in_file);
    htsFormat type = *hts_get_format(fp);
    if (type->format  != vcf && type->format != bcf)
	error("Unsupported input format! %s", input_fname);
    int type = FT_VCF;
    if (out_type != 0) {
	switch (out_type[0]) {
	    case 'b':
		type = FT_BCF_GZ; break;
	    case 'u':
		type = FT_BCF; break;
	    case 'z':
		type = FT_VCF_GZ; break;
	    case 'v':
		type = FT_VCF; break;
	    default :
		error("The output type \"%s\" not recognised\n", out_type);
	};
    }

    double c0 = get_time();    
    const bcf_hdr_t *hdr = bcf_hdr_read(fp);    
    /* duplicate header file and sync new tag in the output header .
     * assuming output is stdout in Vcf format in default. */
    bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);	
    htsFile *fout = out_file == 0 ? hts_open("-", hts_bcf_wmode(type)) : hts_open(out_file, hts_bcf_wmode(type));
    
    /* init gene_predictions or refgene database, hold tabix index cache and faidx cache in memory */
    bcf1_t *line = bcf_init();
    while ( bcf_read(fp, hdr, line) ) {
	anno_hgvs_core(line);
	bcf_write1(fout, hdr_out, line);
    }
    
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    bcf_hdr_destroy(hdr_out);
    hts_close(fout);
    hts_close(fp);
    double c1 = get_time();
    fprintf(stderr, "Run time: %ds\n", c1 -c0);
    return 0;
}
