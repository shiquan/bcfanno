#ifndef HGVSLIB_HEADER
#define HGVSLIB_HEADER
#include <stdio.h>
#include <stdlib.h>
#include "sequence.h"

#define SPLIT_RANGE 3

enum func_region_type {
    _func_region_promote_to_int = -1,
    func_region_unknown,
    //func_region_split_sites,
    func_region_cds,
    func_region_noncoding,
    func_region_intron,
    func_region_utr5,
    func_region_utr3,
    func_region_intergenic,
    // If variants across more than one exon, define large fragment.
    func_region_large,
};

struct var_func_type {
    enum func_region_type func;
    enum var_type vartype;
    // Exon or intron count. Start from 1.
    int count;
    // Location of amino first influenced. Sometime a indel may change the amino sequence after severl condons.
    int loc_amino;
    // Original amino acid and mutated amino acid. Check amino type only if variants happed in cds region.
    int ori_amino;
    int mut_amino;
    // Terminal sequence after loc_amino.
    int fs;
};

struct hgvs_name {
    char *name1; // transcripts name, locus name
    char *name2; // gene name or null 

    // The gene position format is same with genepred dna reference offset consist of two parts: the offset value
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
    int pos;
    int end_pos;
    // location for function region.
    int loc;
    int end_loc;
    
    int offset;
    int end_offset;
    // If strand is '-', convert sequence to complement strand. Default is '+'.
    char strand;
};
struct hgvs_core {
    struct hgvs_name name;
    struct var_func_type type;
};

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
    var_type_del,
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
    char *chrom;
    int32_t start; // position on genome
    // for var_type_snp end == start, for var_type_dels end > start, for var_type_ins end = start +1, 
    // for delvar_type_ins end > start 
    int32_t end;
    // Assume all the ref and alt sequence located in the '+' strand, for some transcripts located in the '-'
    // strand, convert the sequence to complement strand therefore.
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
    int tandam_repeat_number;
    int tandam_repeat_number_ori;

    int l, m;
    struct hgvs_core *a;
};

extern void hgvs_des_clear(struct hgvs_des *des);
extern int init_hgvs_spec(const char *fname, const char *fasta);
extern int set_transcripts_list(const char *fname);
extern int set_genes_list(const char *fname);
extern void hgvs_spec_destroy();

// Parse the description name. 
extern int setter_description(const char *name, int _pos, char *ref, char *alt);
extern int parse_hgvs_name(const char *name);
extern struct hgvs_des *fill_hgvs_name();

#endif
