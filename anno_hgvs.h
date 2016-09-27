#ifndef VCFANNO_HGVS_HEADER
#define VCFANNO_HGVS_HEADER

#include "anno.h"

// for each transcript, the region span several exon partial regions, use exon_pair[] stand each exons
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

#define OFFSET_BITWISE   4
#define REG_NONCODING  1
#define REG_CODING     2
#define REG_UTR5       4
#define REG_UTR3       8
#define REG_MASK       0xF

struct exon_offset_pair {
    uint32_t start;
    uint32_t end;
};

struct genepred_line {
    // mark the memory is allocated or not inside this struct,
    // clear == 0 for empty, clear == 1 recall clear_genepred_line to free it.
    int clear;
    // chromosome names, usually names is start with "chr", like chr1, chr2, ..., chrM, which is UCSC style.
    // however, for other insitutions they have different preference, Ensemble like to use ensemble id, NCBI like
    // use RefSeq accession to represent the chromosome and updated version. Please make sure all the databases have
    // same nomenclature before annotation.
    char *chrom;
    
    uint32_t txstart;
    uint32_t txend;
    // '-' for minus, '+' for plus strand
    char strand;
    // usually gene name, or ensemble gene id
    char *name1;
    // transcript name or null for some gene_pred file
    char *name2; 
    uint32_t cdsstart;
    uint32_t cdsend;
    uint32_t exoncount;
    struct exon_pair *exons;
    struct exon_offset_pair *dna_ref_offsets;
};

enum func_region_type {
    func_region_unknown,
//    func_region_split_sites,
    func_region_cds,
    func_region_intron,
    func_region_utr5,
    func_region_utr3,
    func_region_intergenic,
};

enum var_type {
    _var_type_promoter_to_int = -1,
    var_is_reference,
    var_is_synonymous,
    var_is_inframe_insertion,
    var_is_inframe_deletion,
    var_is_stop_gained,
    var_is_stop_lost,
    var_is_stop_retained,
    var_is_splice_donor,
    var_is_splice_acceptor,
};

static inline char *var_type_string(enum var_type type)
{
    static const char vartypes[10] = {
        "reference",
        "synonymous",
        "inframe insertion",
        "inframe deletion",
        "stop gained",
        "stop lost",
        "stop retained",
        "splice donor",
        "splice acceptor",
        NULL,
    };
    return vartypes[type];
}

struct var_func_type {
    enum func_region_type func;
    enum var_type vartype;
    int count; // for cds or intron count
};

//  hgvs_core keeps transcript/protein name, the format of hgvs name is construct by
//                                    l_name   l_name2  l_type
//                                    |        |        |
//  [transcript|protein|gene|ensemble](gene_id):"prefix".postion...
//
//  l_name          l_type
//  |               |
//  [chrom]:"prefix".postion...
struct hgvs_core {    
    uint16_t l_name; // name offset in the data cache, for chrom l_name == 0;
    uint16_t l_name2; // name2 offset in data cache, if no name2, l_name2 should be 0
    uint16_t l_type; // the byte after type offset, usually offset of '.'
    kstring_t str;
    struct var_func_type type;   
};

// Dynamic allocate and init technology.
// For struct { int l, m, i; void *a }, a is the cached array of predefined type. l for used length, m for max length,
// i for inited length. And i always >= l, and m always >= i. If l == m, the cache array should be reallocated by a new
// memory size. This structure used to reuse the cache complex structures other than points.
struct hgvs {
    int l, m, i; 
    struct hgvs_core *a;
};
struct hgvs_cache {
    int l, m, i; // l == n_allele -1
    struct hgvs *a;
};
#define HGVS_CACHE_INIT { 0, 0, 0, 0 }
// HGVS nomenclature : 
// DNA recommandations 
// - substitution variant, 
// Format: “prefix”“position_substituted”“reference_nucleoride””>”new_nucleoride”, e.g. g.123A>G
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
// hgvs_variant_type is different from var_type, for var_type is more concern about function, but
// hgvs_variant_type is just describe the change on DNA level
enum hgvs_variant_type {
    var_type_nonref = -1, // for gatk <NONREF> allele
    var_type_ref = 0,
    var_type_snp,
    var_type_dels,
    var_type_ins,
    var_type_delins,
    var_type_copy, // dup for ins
    var_type_complex, 
    var_type_unknow,
};
// descriptions struct of hgvs name
struct hgvs_des {
    // if type ==  var_type_ref, hgvs name generater will skip to construct a name, only type will be inited, 
    // DONOT use any other value in this struct then 
    enum hgvs_variant_type type;
    uint8_t strand; // for plus strand start <= end, for minus start >= end
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

// Refgene gene names, transcript names, HGVS names generate
struct genepred_memory {
    int rid;
    uint32_t start;
    uint32_t end;
    // l for used length, m for max length, i for inited length
    int l, m, i;
    struct genepred_line *a;
};

#define GENEPRED_MEMORY_INIT {.rid = -1, .start = 0, .end = 0, .l = 0, .m = 0, .i = 0, .a = 0 }

// map each column right for each format.
// in default refgene have one more `bin` column in the first column than genepred format
struct genepred_format {
    int chrom;
    int name1;
    int name2; 
    int strand;
    int txstart;
    int txend;
    int cdsstart;
    int cdsend;
    int exoncount;
    int exonstarts;
    int exonends;
};

// See UCSC website for more details about refgene format, there is no 
// version information in refgene file in default, however, we can retrieve
// this version number from other file in the same directrary. But remeber 
// to check the names accordly in the refrna file.
struct refgene_options {
    int refgene_is_inited;
    // gene prediction file, include exons and cds regions, get it from UCSC, see our manual for details
    const char *genepred_fname;
    // transcripts sequences in fasta format, notice the name of transcripts should be exists in gene pred file
    const char *refseq_fname;
    // file handler of genepred
    htsFile *fp;
    // refgene should be sorted and indexed by tabix
    tbx_t *genepred_idx;
    // if refseq_fname is provided, check_refseq == 1, else check_refseq == 0
    int check_refseq;
    // faidx of refseq
    faidx_t *refseq_fai;
    // alias to hdr point
    bcf_hdr_t *hdr_out;
    // transcripts list for screening datasets    
    const char *trans_list_fname;
    // gene list for screening datasets
    const char *genes_list_fname;
    int screen_by_genes;
    int screen_by_transcripts;
    void *genehash;
    void *transhash;
    char *columns;
    int n_cols;
    struct anno_col *cols;
    // memory pool
    struct genepred_memory buffer;
    struct hgvs_cache cache;
};

// generate hgvs name
// for col keys, must be one of HGVSDNA, Gene, Transcript
extern int setter_hgvs_names(struct refgene_options *opts, bcf1_t *line, struct anno_col *col);
extern int anno_refgene_core(struct refgene_options *opts, bcf1_t *line);
extern int hgvs_bcf_header_add_gene(bcf_hdr_t *hdr);
extern int hgvs_bcf_header_add_hgvsdna(bcf_hdr_t *hdr);
extern int hgvs_bcf_header_add_trans(bcf_hdr_t *hdr);
extern void set_format_refgene();
extern void set_format_genepred();
//extern int refgene_opts_destroy(struct refgene_options *opts);
extern int refgene_columns_parse(struct refgene_options *opts, char *columns);
extern int refgene_set_refseq_fname(struct refgene_options *opts, char *fname);
extern int refgene_set_refgene_fname(struct refgene_options *opts, char *fname);
extern int refgene_set_trans_fname(struct refgene_options *opts, char *fname);
extern int refgene_set_genes_fname(struct refgene_options *opts, char *fname);
extern int refgene_options_init(struct refgene_options *opts);
extern int refgene_options_destroy(struct refgene_options *opts);

#endif
