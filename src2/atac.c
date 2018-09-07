// atac.c
// Transcription factors bind the open chromatin region and may further influence the histone
// states and gene regulation (McVicker et al. Science 2013). A genetic variant happened at binding 
// regions may contribute to weakness of motifs binding, this variant 'disrupt' the motif. We will
// observe mutated alleles disappeared or reduced at genome accessible regions, we call this as allele
// imbalance or allele specific peak. This program is trying to calculate the allele counts from ATAC
// reads for each genomic position at ATAC peaks.


#include "utils.h"
#include "bed_utils.h"
#include "wrap_pileup.h"
#include "htslib/vcf.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"


enum abt {
    allele_imbalance,
    allele_balance,
    allele_lowcover,
};

static const char *hts_bcf_wmode(int file_type)
{
    if ( file_type == FT_BCF ) return "wbu";    // uncompressed BCF
    if ( file_type & FT_BCF ) return "wb";      // compressed BCF
    if ( file_type & FT_GZ ) return "wz";       // compressed VCF
    return "w";                                 // uncompressed VCF
}

struct args {
    const char *vcf_fname;
    const char *bam_fname;
    const char *bed_fname;
    const char *fasta_fname;
    const char *out_fname;
    
    struct bedaux *bed;
    
    // vcf input/output
    htsFile *fp_input;
    htsFile *fp_out;

    bcf_hdr_t *bcf_hdr;
    int output_type;
    int qual_thres;
    
} args = {
    .vcf_fname = NULL,
    .bam_fname = NULL,
    .bed_fname = NULL,
    .fasta_fname = NULL,
    .out_fname = NULL,
    
    .bed = NULL,
    .fp_input = NULL,
    .fp_out = NULL,
    .bcf_hdr = NULL,
    .output_type = 0,
    .qual_thres = 20,
};

// Chromatin Accessible Auxiliary, CAA
// For each peak we cache all covered reads and piled for each position.
struct CAA {
    int tid;
    int start, end; // current block
    samFile *fp; // for each CAA, hold a file handler, thread safe
    hts_idx_t *sam_idx; // index of alignment file
    hts_itr_t *sam_itr; // update for each position?
    bam_hdr_t *sam_hdr; // 

    bcf_hdr_t *bcf_hdr; // point to args::bcf_hdr
    
    //struct plp_ref *ref;
    bam_plp_t plp_iter; // returned by bam_plp_init
    int n_plp;
    //faidx_t *fai;
    //char *ref;
    //int ref_len;

    struct bedaux *bed; // point to args::bed, do NOT free it.
    struct args *args; // point to args, for accessing some parameters

    int id; // GT id
    int sample_id;

    int last_pos;
    const bam_pileup1_t *last_plp;
    
    char *name; // sample name in vcf
};

#define CAA_seqname(_C,_T) _C->sam_hdr->target_name[_T]

void CAA_destroy(struct CAA *CAA)
{
    //fai_destroy(CAA->fai);
    bam_hdr_destroy(CAA->sam_hdr);
    sam_close(CAA->fp);
    hts_idx_destroy(CAA->sam_idx);

    //if ( CAA->->ref[0] ) free(CAA->ref->ref[0]);
    //if ( CAA->ref->ref[1] ) free(CAA->ref->ref[1]);
    if ( CAA->sam_itr ) free(CAA->sam_itr);

    free(CAA);
}

struct CAA *CAA_init(const char *alignment_fname)//, const char *reference_fname)
{
    struct CAA *CAA = malloc(sizeof(*CAA));
    memset(CAA, 0, sizeof(*CAA));

    CAA->tid = -1;
    htsFile *fp = hts_open(alignment_fname, "r");
    if (fp == NULL) error("Failed to open %s: %s.", alignment_fname, strerror(errno));
    htsFormat type = *hts_get_format(fp);
    hts_close(fp);

    CAA->fp = sam_open_format(alignment_fname, "rb", &type);
    CAA->sam_hdr = sam_hdr_read(CAA->fp);
    CAA->sam_idx = sam_index_load(CAA->fp, alignment_fname);
    //CAA->fai = fai_load(reference_fname);

    if ( CAA->fp == NULL ) error("%s: %s.", alignment_fname, strerror(errno));
    if ( CAA->sam_hdr == NULL ) error("Failed to read BAM header of %s.", alignment_fname);
    if ( CAA->sam_idx == NULL ) error("Failed to load index of %s.", alignment_fname);
    //if ( CAA->fai == NULL ) error("Failed to load index %s.", reference_fname);

    //static plp_ref_t ref = {{NULL,NULL}, {-1,-1}, {0,0}};
    //CAA->ref = &ref;
    return CAA;
}

const bam_pileup1_t *CAA_variant_pos(struct CAA *CAA, int _tid, int _pos);

static int plp_func(void *data, bam1_t *b)
{
    int ret;

    struct CAA *CAA = (struct CAA*)data;
    struct args *args = CAA->args;

    for ( ;; ) {
        ret = CAA->sam_itr == NULL ? sam_read1(CAA->fp, CAA->sam_hdr, b) : sam_itr_next(CAA->fp, CAA->sam_itr, b);
        if ( ret < 0 ) break;
        if ( b->core.tid < 0  || (b->core.flag & BAM_FUNMAP)) continue;
        if ( b->core.qual < args->qual_thres || (b->core.flag & BAM_FQCFAIL)) continue;
        if ( b->core.flag & BAM_FDUP ) continue;
        break;
    }

    return ret;
}

const bam_pileup1_t *CAA_new_region_init(struct CAA *CAA, int tid, int pos)
{
    int start, end;
    int ret, has_ref;
    struct bedaux *bed = CAA->bed;
    ret = bed_position_covered(bed, CAA_seqname(CAA,tid), pos, &start, &end);
    
    if ( ret == 0 ) return NULL; // Not located in open accessible regions
   
    if ( CAA->sam_itr ) hts_itr_destroy(CAA->sam_itr);

    CAA->tid = tid;
    CAA->start = start;
    CAA->end = end;   

    CAA->sam_itr = bam_itr_queryi(CAA->sam_idx, tid, start, end);
    CAA->plp_iter = bam_plp_init(plp_func, CAA);
    
    /*
    has_ref = plp_get_ref(CAA->ref, CAA->sam_hdr, CAA->fai, tid, &CAA->ref, &CAA->ref_len);
    
    if (has_ref == 0 ) {
        warning("No such chromosome %s at reference.", CAA_seqname(CAA,tid));
        return NULL;
    }
    */
    return CAA_variant_pos(CAA, tid, pos);
}

// 
const bam_pileup1_t *CAA_variant_pos(struct CAA *CAA, int _tid, int _pos)
{
    // new chromosome or new region
    if ( CAA->tid != _tid || CAA->end < _pos )
        return CAA_new_region_init(CAA, _tid, _pos);

    // last position cached
    if ( CAA->last_pos == _pos ) {
        debug_print("Return last cached plp. No test.");
        return CAA->last_plp;
    }
    // if ( CAA->has_ref == 0 ) return NULL;

    int n_plp;
    int tid;
    int pos;
    for ( ;; ) {
        
        const bam_pileup1_t *plp = bam_plp_auto(CAA->plp_iter, &tid, &pos, &n_plp);
        //debug_print("1: %s\t%d", CAA_seqname(CAA,_tid), _pos);
        if ( plp == NULL ) break; // ends
        //debug_print("2: %s\t%d", CAA_seqname(CAA,tid), pos);
        pos = pos +1; // 0 based to 1 based
        if ( pos < _pos ) 
            continue;
        else if ( pos == _pos ) {
            CAA->n_plp = n_plp;
            return plp; // that's it
        }
        else {
            // there is no reads cover this location, plp_auto will skip uncover locs
            CAA->n_plp = n_plp;
            CAA->last_pos = pos;
            CAA->last_plp = plp;
            return NULL;
        }
        // error("Input record is not sorted ? %s:%d,%d", CAA_seqname(CAA,tid), _pos,pos); 
    }
    return NULL;
}

int anno_vcf_atac(struct CAA *CAA, bcf1_t *l, const bam_pileup1_t *plp)
{
    int i;
    // char *ref
    // chr pos ref alt ref_depth alt_depth
    bcf_unpack(l, BCF_UN_ALL);
    uint32_t *d = calloc(l->n_allele, sizeof(*d));
    float *f = calloc(l->n_allele, sizeof(*f));
    uint32_t sum = 0;
    for ( i = 0; i < CAA->n_plp; i++ ) {
        const bam_pileup1_t *p = plp+i;
        if ( !p->is_del ) {
            if ( p->qpos >= p->b->core.l_qseq ) continue;
            char c = "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(p->b), p->qpos)];
            if (c == 'N' ) continue;
            //debug_print("%d\t%d\t%c",p->b->core.pos,p->qpos,c);
            int j;
            for ( j = 0; j < l->n_allele; ++j ) {
                if ( c== *l->d.allele[j] ) { d[j]++; sum++; break; }               
            }
        }
        else {

            // for indels
        }
    }
    if ( sum > 0 ) {
        for ( i = 0; i < l->n_allele; ++i ) f[i] = (float)d[i]/sum;
        bcf_update_format_float(CAA->bcf_hdr, l, "PeakAF", f, l->n_allele);
    }
    // printf("%s\t%d\t%s\t%s\t%d\t%d\n", CAA->bcf_hdr->id[BCF_DT_CTG][l->rid].key, l->pos+1, l->d.allele[0],
    // l->n_allele > 1? l->d.allele[1] : ".", d[0], l->n_allele>1? d[1] : 0);
    bcf_update_format_int32(CAA->bcf_hdr, l, "PeakAC", d, l->n_allele);

    return 0;
}
int bcf2bam_rid(bcf_hdr_t *bcf_hdr, bam_hdr_t *sam_hdr, int tid)
{
    const char *n = bcf_hdr->id[BCF_DT_CTG][tid].key;
    int i;
    for ( i = 0; i < sam_hdr->n_targets; ++i ) {
        if ( strcmp(n, sam_hdr->target_name[i]) == 0 ) 
            return i;
    }
    return -1;
}
int anno_vcf_atac_main(struct CAA *CAA)
{
    //int tid, pos;
    //char *ref = NULL;
    //int ref_len = 0;

    bcf1_t *line = bcf_init();

    for (;;) {
        if ( bcf_read(args.fp_input, args.bcf_hdr, line)!= 0 ) break;
        if ( line->rid == -1 ) goto out;
        
        // if ( bcf_get_variant_types(line) == VCF_REF ) goto output_line;
        // if ( bed_region_covered(args.bed, args.bcf_hdr->target_name[tid], line->pos+1, line->pos+1) == 0 ) goto out;
        int tid = bcf2bam_rid(args.bcf_hdr, CAA->sam_hdr, line->rid);
        const bam_pileup1_t *plp = CAA_variant_pos(CAA, tid, line->pos+1);
        if ( plp != NULL && CAA->n_plp != 0 )
            anno_vcf_atac(CAA, line, plp);

      out:
        bcf_write1(args.fp_out, args.bcf_hdr, line);
        continue;
    }

    return 0;
}

int usage()
{
    // SAMPLE_TREATMENT_AF
    // SAMPLE_TREATMENT_INTERUPT
    // SAMPLE_
    fprintf(stderr, "bcfanno_atac\n");
    fprintf(stderr, "  -bed    region.bed\n");
    fprintf(stderr, "  -bam    aln.bam\n");
    fprintf(stderr, "  -vcf    wgs.vcf\n");
    //fprintf(stderr, "  -fasta  ref.fa\n");
    fprintf(stderr, "  -s      sample name, if not set annotated to first sample in the VCF\n");
    return 1;
}
int parse_args(int argc, char **argv)
{
    int i;
    const char *output_type = 0;
    // if ( argc == 1 ) return usage();
    
    for ( i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        if ( strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0 ) return usage();

        if ( strcmp(a, "-bed") == 0 )
            var = &args.bed_fname;
        else if ( strcmp(a, "-bam") == 0 )
            var = &args.bam_fname;
        else if ( strcmp(a, "-vcf") == 0 )
            var = &args.vcf_fname;
        else if ( strcmp(a, "-fasta") == 0 )
            var = &args.fasta_fname;
        else if ( strcmp(a, "-O") == 0 )
            var = &output_type;
        else if ( strcmp(a, "-o") == 0 )
            var = &args.out_fname;
        
        if ( var != 0 ) {
            if ( i == argc ) error("Missing an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        error("Unknown argument: %s, use -h see help information.", a);
    }

    if ( args.vcf_fname == NULL ) error("Parameter -vcf is required. Use -h for more information");
    if ( args.bam_fname == NULL ) error("Parameter -bam is required. Use -h for more information");
    if ( args.bed_fname == NULL ) error("Parameter -bed is required. Use -h for more information");
    //if ( args.fasta_fname == NULL ) error("-fasta is required.");

    args.bed = bedaux_init();
    bed_read(args.bed, args.bed_fname);
    if ( args.bed->flag & bed_bit_empty) error("Could not load BED file. %s.", args.bed_fname);
    bed_merge(args.bed);

    args.fp_input = hts_open(args.vcf_fname, "r");
    if ( args.fp_input == NULL ) error("%s : %s.", args.vcf_fname, strerror(errno));
    
    int out_type = FT_VCF;
    if ( output_type != 0 ) {
	switch (output_type[0]) {
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
    args.fp_out = args.out_fname == NULL ? hts_open("-", hts_bcf_wmode(out_type)) :
        hts_open(args.out_fname, hts_bcf_wmode(out_type));
    
    args.bcf_hdr = bcf_hdr_read(args.fp_input);

#define BRANCH(_key, _description) do {                                 \
        int id;                                                         \
        id = bcf_hdr_id2int(args.bcf_hdr, BCF_DT_ID, _key);             \
        if (id == -1) {                                                 \
            bcf_hdr_append(args.bcf_hdr, _description);                 \
            bcf_hdr_sync(args.bcf_hdr);                                 \
            id = bcf_hdr_id2int(args.bcf_hdr, BCF_DT_ID, _key);         \
            assert(bcf_hdr_idinfo_exists(args.bcf_hdr, BCF_HL_FMT, id)); \
        }                                                               \
    } while(0)

    BRANCH("PeakAC", "##FORMAT=<ID=PeakAC,Number=R,Type=Integer,Description=\"Counts in peak region respect to each allele.\">");
    BRANCH("PeakAF", "##FORMAT=<ID=PeakAF,Number=R,Type=Float,Description=\"Count frequency in peak region respect to each allele.\">");
#undef BRANCH

    bcf_hdr_write(args.fp_out, args.bcf_hdr);
    
    return 0;
}
void memory_release()
{
    hts_close(args.fp_input);
    hts_close(args.fp_out);
}
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) ) return 1;
    
    struct CAA *CAA = CAA_init(args.bam_fname);
    CAA->args = &args;
    CAA->bed = args.bed;
    CAA->bcf_hdr = args.bcf_hdr;
    CAA->id = bcf_hdr_id2int(args.bcf_hdr, BCF_DT_ID, "GT");    
    if ( CAA->id <0 ) error("No GT tag found at input VCF.");
    
    anno_vcf_atac_main(CAA);

    // release memory
    CAA_destroy(CAA);
    memory_release();
    
    return 0;
}
