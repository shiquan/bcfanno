#ifndef GEA_HEADER
#define GEA_HEADER

#include "utils.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"

// The design of GEA header structure borrow a lot of definitions from VCF/BCF format and UCSC genepred format,
// this design makes GEA flexible and extensible. Meanwhile, some definitions in VCF/BCF format are not suitable
// for GEA, like FILTER/REF/ALT/GT .. have been discarded.

/*
  demo 1 : stores sample peaks in GEA format.

  ##fileformat=GenomeElementAnnotation v1.0
  ##contig=<ID=chr1,length=199999999>
  ##bioType=<ID="peak",Description="Peaks.">
  ##cellType=<ID="CD38",Description="">
  ##individual=<ID="P1",Phenotype="Lung cancer.",Description="Patient.">
  ##individual=<ID="A1",Phenotype="None.",Description="Health control.">
  ##experiment=<ID="ChIPseq+H3K23me3",Description="">
  ##experiment=<ID="ATACseq",Description="">
  ##sample=<ID="P1_CD38_H3K23me3",individual="P1",cellType="CD38",experiment="ChIPseq+H3K23me3",Description="CD38 cell from patient P1.">
  ##sample=<ID="P1_CD4_H3K23me3",individual="P1",cellType="CD4",experiment="ChIPseq+H3K23me3",Description="CD4 cell from patient P1.">
  ##sample=<ID="A1_fibroblast_H3K23me3",individual="A1",cellType="fibroblast",featureType="H3K23me3",experiment="ChIPseq",Description="Fibroblast cell.">
  ##comments
  ##FORMAT=<ID="expLevel",Number=1,Type=Float,Range=0:100,Description="Expression level.">
  ##FORMAT=<ID="peaksLevel",Number=block,Type=Float,Range=0:100,Description="Normalized peaks level.">
  #chrom chromStart chromEnd name bioType geneName strand cdsStart cdsEnd blockCount blockStarts blockEnds alignment_state INFO FORMAT P1_CD38_H3K23me3 P1_CD4_H3K23me3 ...

  *******************
  demo 2 : stores annotations in GEA database
  
  ##fileformat=GenomeElementAnnotation v1.0
  ##contig=<ID=chr1,length=199999999>
  ##bioType=<ID="open_chromatin",Description="Open chromatin.">
  ##bioType=<ID="TFBS",Description="Transcript factor binding sites or CTCF binding sites.">
  ##bioType=<ID="enhancer",Description="Enhancer.">
  ##bioType=<ID="promoter_flanking_region",Description="Promoter flanking region.">
  ##bioType=<ID="promoter",Description="Promoter.">
  ##bioType=<ID="mRNA",Description="Coding RNA.">
  ##bioType=<ID="ncRNA",Description="Noncoding RNA.">
  ##bioType=<ID="gene",Description="Gene.">
  ##comments
  ##INFO=<ID="PMT",Number=1,Type=Int,Description="Pathogenic mutations in this transcript.">
  ##INFO=<ID="PME",Number=block,Type=Int,Description="Pathogenic mutations in each exon.">
  ##INFO=<ID="CVT",Number=1,Type=Int,Description="Common variants in this transcript.">
  ##INFO=<ID="CVE",Number=block,Type=Int,Description="Common variants in each exon.">
  #chrom chromStart chromEnd name bioType geneName strand cdsStart cdsEnd blockCount blockStarts blockEnds alignment_state INFO
*/
#define GEA_SAMPLE_COL  14

// variant length for each tag
#define GEA_VL_FIXED     0 // fixed length
#define GEA_VL_VAR       1 // varied
#define GEA_VL_BLOCK     2 // value consist of records seperated by comma in which one record for each block

// type of value for each tag
#define GEA_TYPE_STR    0
#define GEA_TYPE_REAL   1
#define GEA_TYPE_INT    2
#define GEA_TYPE_FLAG   3

// Dictionary id
// This design is simliar to VCF structure. The header keeps three dictionaries. First dict keeps IDs
// in "cellType/bioType/INFO/FORMAT" lines, the second keeps the sequence names and lenths in the contig lines and the
// last keeps the tissue/cell type/sample names.

// Dictionary type
#define GEA_DICT_ALL    8

#define GEA_DT_ID       0
#define GEA_DT_BIOTYPE  1
#define GEA_DT_CELLTYPE 2
#define GEA_DT_INDVI    3
#define GEA_DT_EXPR     4
#define GEA_DT_CTG      5
#define GEA_DT_SAMPLE   6
#define GEA_DT_COMMENT  7

struct gea_hrec {
    int type; // one of GEA_HL_* type
    char *key; // The part before '='. Should be contig/cellType/bioType/INFO/FORMAR/SAMPLE or genic key.
    char *value; // Set only for generic lines. NULL for contig/cellType/bioType/INFO/FORMAT/SAMPLE.
    int n_key;
    char **keys, **vals; // The key=value pairs.
};

// header line
#define GEA_HL_BIOTYPE   0
#define GEA_HL_INFO      1
#define GEA_HL_FMT       2
#define GEA_HL_CTG       3
#define GEA_HL_STR       4 // structured header line TAG=<A=..,B=..>
#define GEA_HL_GEN       5 // generic header line TAG=.
#define GEA_HL_CELL      6
#define GEA_HL_INDIVI    7
#define GEA_HL_EXPR      8
#define GEA_HL_SAMPLE    9

enum gea_biotype_predefined {
    biotype_not_support = -1,
    biotype_genome = 0,
    biotype_gene,
    biotype_coding_transcript,
    biotype_noncoding_transcript,
    biotype_peak,
    biotype_regulatory_region,
    biotype_tfbs,
    biotype_enhancer,
    biotype_promotor,
    biotype_promotor_flanking_region,
    biotype_open_chromatin,
    //biotype_topological_domain,
    //biotype_predicted_functional_region,
};

static const char *gea_biotype_predefined_string[] = {
    "genome",
    "gene",
    "mRNA",
    "ncRNA",
    "peak",
    "regulatory_region",
    "TFBS",
    "enhancer",
    "promotor",
    "promotor_flanking_region",
    "open_chromatin",
    NULL,
};

struct gea_id_info {
    int id;
    // Predefined types.
    // INFO/FOMAT, Number:20, var:4, Type:4, ColType:4, consistance with VCF/INFO
    // Not like var of VCF, GEA do not need genotype information, but blocks type is support
    // for CONTIG, info is the length
    uint32_t info;
    // bioType
    int biotype;
    // SAMPLE
    int indiv, celltype, exper;

    // other key=value infomation, i.e. Description
    struct gea_hrec *hrec;
};

struct gea_id_pair {
    const char *key;
    const struct gea_id_info *val;
};

struct gea_hdr {
    int32_t n[GEA_DICT_ALL]; // n : the size of dict block in use
    int32_t m[GEA_DICT_ALL]; // m : allocated size of each dict block
    char *version; // version of gea format, prserved for further update
    struct gea_id_pair *id[GEA_DICT_ALL];
    void *dict[GEA_DICT_ALL];
    char **samples;
    // comments
    struct gea_hrec **hrec;
    int n_hrec, dirty;
    kstring_t mem;
};

#define gea_hdr_id2type(hdr, id) ((hdr)->id[GEA_DT_ID][id].val->info>>4 & 0xf)
#define gea_hdr_id2length(hdr, id) ((hdr)->id[GEA_DT_ID][id].val->info>>8 & 0xf)
#define gea_hdr_id2number(hdr, id) ((hdr)->id[GEA_DT_ID][id].val->info>>12)
#define gea_hdr_id2coltype(hdr, id) ((hdr)->id[GEA_DT_ID][id].val->info&0xf)
#define gea_hdr_idinfo_exists(hdr, id) ((id < 0 || gea_hdr_id2coltype(hdr, id) == 0xf) ? 0 : 1)

// GEA record
struct gea_format {
    int id; // id : numeric tag id, the corresponding string is gea_hdr::id[GEA_DT_ID][$id].key
    int n, size, type; // n : number of values per-sample; size : number of bytes per-sample; type : one of GEA_BT_* types
    uint8_t *p;        // 
    uint32_t p_len;
    uint32_t p_off:31, p_free:1;
};

struct gea_dec {
    int m_info, m_fmt;
    bcf_info_t *info; // INFO
    struct gea_format *fmt; // FORMAT and individual sample
    int shared_dirty; // if set, share.s must be recreated on GEA output
    int indiv_dirty;  // if set, indiv.s must be recreated on GEA output
};

#define GEA_ERR_CTG_UNDEF   1
#define GEA_ERR_TAG_UNDEF   2
#define GEA_ERR_NCOLS       4
#define GEA_ERR_LIMITS      8
#define GEA_ERR_CHAR        16
#define GEA_ERR_CTG_INVALID 32
#define GEA_ERR_TAG_INVALID 64

// CIGARS
#define CIGAR_UNKNOWN_BASE '*'
#define CIGAR_MATCH_BASE   'M'
#define CIGAR_DELETE_BASE  'D'
#define CIGAR_INSERT_BASE  'I'

// The fields cigar packed counts and cigar type, as follows:
// uint32_t type:4,count:28
//   offset   bits  value
//   0        4     cigar_type
//   5        32    cigar_counts
//
#define CIGAR_PACKED_FIELD 4
#define CIGAR_UNKNOWN_TYPE 0
#define CIGAR_MATCH_TYPE   1
#define CIGAR_DELETE_TYPE  2
#define CIGAR_INSERT_TYPE  3
#define CIGAR_MASK_TYPE    0xf

//
// The gea_record structure corresponds to one GEA line. Reading from GEA fule is slower because the string
// is first to be parsed, packed into GEA line (done in gea_parse), then unpacked into internal gea_record
// structure. If it is known in advance that some of the fields will not be required (notably the sample
// columns), parsing of these can be skipped by setting max_unpack appropriately.
// Similarly, it is fast to output a GEA line because the columns (kept in shared.s, indiv.s, etc.) are
// written directly by gea_write, whereas a GEA line musy be formatted in gea_format.
struct gea_coding_transcript {
    int *loc[2];
    // UTR5 length of coding transcript, for both strands utr5 should be always on the upstream of coding region
    int utr5_length;
    // CDS length should be length of UTR5 plus CDS here
    int cds_length;
    // Length of this transcript, intron emitted.
    int reference_length;
};

enum strand {
    strand_is_unknown = -1, // "."
    strand_is_plus, // "+"
    strand_is_minus, // "-"
    strand_is_both, // "+-"
};
struct gea_record {
    int32_t rid; // chrom
    int chromStart; // chromStart
    int chromEnd; // chromEnd
    enum strand strand;
    char *name; // name of this record, for Gene could be gene name, for transcript could be [XN][MR] id
    //enum gea_biotype_predefined bioType; // biological type of this record, this value should be well defined
    // if bioType == biotype_not_support, init type string
    // For now GEA format only support severl biological types, new type (not support) will be set here.
    //char *biotype_not_support;
    int biotype;
    // Set to GEA_UN_BLOCK, GEA_UN_INFO,GEA_UN_CIGAR,GEA_UN_FMT to boost performance of gea_parse when some
    // of the fields won't be needed
    //int max_unpack; 

    // GEA_UN_BLOCK
    char *geneName;
    int cStart; // coding start for mRNA, core start for functional region
    int cEnd; // coding end for mRNA, core end for functional region        
    int blockCount; 
    int *blockPair[2];

    // GEA_UN_CIGAR
    int n_cigar;
    int *cigars;

    uint32_t n_info:16, n_fmt:16;
    uint32_t n_sample;
    kstring_t shared, indiv;
    
    struct gea_dec d;
    struct gea_coding_transcript c; // only set for coding transcript and GEA_UN_BLOCK unpacked
    
    int unpacked; // remember what has been unpaked to allow calling gea_unpack() repeatedly without redoing the work

    int errcode; // one of GEA_ERR_* codes
};

#define gea_hdr_nsamples(hdr) (hdr)->n[GEA_DT_SAMPLE]

// gea_unpack() - unpack/decode a GEA record
//#define GEA_UN_BLOCK 1
#define GEA_UN_CIGAR 1
#define GEA_UN_TRANS 2
#define GEA_UN_FMT   4
#define GEA_UN_INFO  8
#define GEA_UN_ALL   (GEA_UN_CIGAR|GEA_UN_FMT|GEA_UN_INFO)

int gea_unpack(const struct gea_hdr *hdr, struct gea_record *b, int which);
int gea_parse(kstring_t *s, const struct gea_hdr *h, struct gea_record *v);
struct gea_record *gea_init();
void gea_destroy(struct gea_record *rec);

// return 0 for gea format, -1 for unknown
int gea_check_format(const char *fn);

struct gea_hdr *gea_hdr_read(htsFile *fp);
int gea_hdr_write(htsFile *fp, const struct gea_hdr *hdr);

// Append new GEA header line, return 0 on success
int gea_hdr_append(struct gea_hdr *hdr, const char *line);

// Read or write GEA record 
int gea_read(htsFile *fp,const struct gea_hdr *hdr, struct gea_record *rec);
int gea_write(htsFile *fp, const struct gea_hdr *hdr, struct gea_record *rec);

struct gea_format *gea_get_fmt(const struct gea_hdr *hdr, struct gea_record *rec, const char *key);
bcf_info_t   *gea_get_info(const struct gea_hdr *hdr, struct gea_record *rec, const char *key);

struct gea_format *gea_get_fmt_id(struct gea_record *rec, const int id);
bcf_info_t   *gea_get_info_id(const struct gea_hdr *hdr, struct gea_record *rec, const int id);

// predefined compact structure for fast access
struct gea_reader {
    htsFile *file;
    tbx_t *idx;
    struct gea_hdr *hdr;
    struct gea_record **buffer;
    int n, m; 
};

struct gea_reader *gea_read_file(const char *fn);
int gea_reader_destroy(struct gea_reader *reader);
int gea_reader_next_line(struct gea_reader *reader, struct gea_record *rec);


#endif
