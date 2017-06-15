#include "utils.h"
#include "hgvs.h"
#include "genepred.h"
#include "number.h"
#include "sort_list.h"
#include "htslib/kstring.h"
#include <string.h>
#include <regex.h>
#include <ctype.h>

// Avoid overhang point.
char str_init[2];

struct hgvs_spec {
    struct genepred_spec *data;
    struct hgvs_des des;
    // To parse hgvs nomen.
    regex_t *exp;
} spec;

struct splice_sites_definition {
    int splice_sites5_upstream;
    int splice_sites5_downstream;
    int splice_sites3_upstream;
    int splice_sites3_downstream;
} ss = {
    .splice_sites5_downstream = 3,
    .splice_sites5_upstream = 3,
    .splice_sites3_upstream = 3,
    .splice_sites3_downstream = 3,
};

void splice_sites_definition_update(int downstream5, int upstream5, int downstream3, int upstream3)
{
    ss.splice_sites5_downstream = downstream5;
    ss.splice_sites5_upstream = upstream5;
    ss.splice_sites3_downstream = downstream3;
    ss.splice_sites3_upstream = upstream3;
}

void hgvs_core_clear(struct hgvs_core *core)
{
    if ( core->name.name1 != NULL)
        free(core->name.name1);
    if ( core->name.name2 != NULL)
        free(core->name.name2);
    memset(core, 0, sizeof(struct hgvs_core));
}
void hgvs_des_clear(struct hgvs_des *des)
{
    int i;
    for ( i = 0; i < des->l; ++i) {
        hgvs_core_clear(&des->a[i]);
    }
    free(des->a);
    if ( des->chrom != NULL )
        free(des->chrom);
    if ( des->ref != NULL && des->ref_length > 0)
        free(des->ref);
    if ( des->alt != NULL && des->alt_length > 0) 
        free(des->alt);    
    memset(des, 0, sizeof(struct hgvs_des));    
}
void hgvs_spec_destroy()
{
    genepred_spec_destroy(spec.data);
    hgvs_des_clear(&spec.des);
    if ( spec.exp != NULL )
        free(spec.exp);
}
int init_hgvs_spec(const char *fname, const char *fasta)
{
    memset(&spec, 0, sizeof(struct hgvs_spec));
    spec.data = genepred_spec_init();
    if ( genepred_load_data(spec.data, fname)  )
        return 1;
    if ( genepred_load_fasta(spec.data, fasta) )
        return 1;
    spec.exp = NULL;
    /* spec.exp = (regex_t*)malloc(sizeof(regex_t)); */
    /* int rv; */
    /* rv = regcomp(spec.exp, "(.*):[gcn][.](.*)", REG_EXTENDED); */
    /* if ( rv != 0 ) { */
    /*     error("regcomp failed with %d.", rv); */
    /* }     */
    return 0;
}

int set_transcripts_list(const char *fname)
{
    if ( genepred_load_trans(spec.data, fname) == 1)
        return 1;
    return 0;
}

int set_genes_list(const char *fname)
{
    if ( genepred_load_genes(spec.data, fname) == 1 )
        return 1;
    return 0;
}

enum nametype {
    name_is_unknown = -1,
    name_is_chromosome,
    name_is_coding_rna,
    name_is_noncoding_rna,
};
static int convert_loc2position(struct genepred_line *line, int location, int offset)
{
    int i;
    if ( line->strand == '+') {
        for ( i = 0; i < line->exon_count; ++i ) {
            if ( location >= line->loc[BLOCK_START][i] && location <= line->loc[BLOCK_END][i])
                break;            
        }
        
        if ( location > line->loc[BLOCK_END][i])
            return 0;
        if ( location == line->loc[BLOCK_START][i] )
            return line->exons[BLOCK_START][i] + offset;
        if ( location == line->loc[BLOCK_END][i] )
            return line->exons[BLOCK_END][i] + offset;
        
        if ( offset != 0 && line->exons[BLOCK_START][i] != location && line->exons[BLOCK_END][i] != location ) {
            error("Inconsist between offset and exon edge. exon location : %d, offset: %d.", location, offset);
        }
        return line->exons[BLOCK_START][i] + location - line->loc[BLOCK_START][i];        
    } else {
        for ( i = 0; i < line->exon_count; ++i ) {
            if ( location >= line->loc[BLOCK_END][i] && location <= line->loc[BLOCK_START][i])
                break;            
        }
        if ( location == line->loc[BLOCK_START][i] )
            return line->exons[BLOCK_START][i] + offset;
        if ( location == line->loc[BLOCK_END][i] )
            return line->exons[BLOCK_END][i] + offset;

        if ( offset != 0 && line->exons[BLOCK_START][i] != location && line->exons[BLOCK_END][i] != location ) {             
            error("Inconsist between offset and exon edge. exon location : %d, offset: %d.", location, offset);           
        }
        return line->exons[BLOCK_START][i] + line->loc[BLOCK_START][i] - location;
    }
    return 0;
}
// ss point the start of string, se point the end of the convert string.
static int parse_position(char *ss, char *se, struct genepred_line *line)
{
    // struct hgvs_des *des = &spec.des;
    
    // Parse the position string. The possible position may be look like 123, *123, -123, 123+12, 123-12, etc.        
    int position = 0;
    int offset = 0;
    // Assume all the position located in the cds region first. For intron region, offset should NOT be 0, and
    // pos_type is the type of nearby region.        
    enum func_region_type pos_type = func_region_cds;
    char *s1 = ss;
    char *s2;
    if ( *s1 == '*' ) {
        pos_type = func_region_utr3;
        ++s1;
    } else if ( *s1 == '-') {
        pos_type = func_region_utr5;
        ++s1;
    }
    
    int i;
    for ( s2 = s1, i = 0; *s2 >= 48 && *s2 <= 57; ++s2, ++i );
    position = str2int_l(s1, i);        
    if ( s2 != se ) {
        for ( i = 0, s1 = s2; s2 != se; ++s2, ++i);
        if ( check_num_likely_l(s1, i) ) {
            offset = str2int_l(s1, i);
        } else {
            error("Failed to parse offset. %s.", s1);
        }
    }

    // Convert functional location to gene location.
    if ( pos_type == func_region_cds ) {
        position += line->utr5_length;
    } else if ( pos_type == func_region_utr5 ) {
        position = line->utr5_length - position + 1;
    } else if ( pos_type == func_region_utr3 ) {
        position = line->reference_length - line->utr3_length + position;    
    } else {
        error("Impossible situation.");
    }
    return convert_loc2position(line, position, offset);
}
static int is_nucletide(char *ss)
{
    if ( ss == NULL )
        return 1;

    if ( *ss == 'A' || *ss == 'a' || *ss == 'C' || *ss == 'c' || *ss == 'G' || *ss == 'g' || *ss == 'T' || *ss == 't' )
        return 0;

    return 1;
}
// Here only support several universal HGVS styles.
// simple convert style:  [ACGT]n>[ACGT]n
// deletion style:       del[ACGT]n
// Insertion Style:      Ins[Acgt]N
// duplicate style:      dup[ACGT]n
// Do NOT support other complex styles, like inv, con, etc.
static int parse_var(char *ss, char *se, int strand, int *ref_length, char **ref, int *alt_length, char **alt, enum hgvs_variant_type *type)
{
    char *s1;
    char *s2;

    s1 = ss;
    s2 = ss;
    
    if ( ss == NULL || strlen(ss) < 3 )
        return 1;

    func_dup_seq func = strand == '+' ? strndup : rev_seqs;

    *type = var_type_ref;
    
    int i;
    if ( s1[0] == 'd' && s1[1] == 'e' && s1[2] == 'l' ) {
        s1 += 3;
        s2 = s1;
        for ( i = 0 ; s1 <= se && is_nucletide(s1) == 0; ++s1, ++i );
        *alt_length = 0;
        *alt = NULL;
        *ref_length = i;
        *ref = func(s2, i);
        *type = var_type_del;
    } else if ( (s1[0] == 'i' && s1[1] == 'n' && s1[2] == 's') || (s1[0] == 'd' && s1[1] == 'u' && s1[2] == 'p') ) {
        s1 += 3;
        s2 = s1;
        for ( i = 0 ; s1 <= se && is_nucletide(s1) == 0; ++s1, ++i );
        *ref_length = 0;
        *ref = NULL;
        *alt_length = i;
        *alt = func(s2, i);
        *type = var_type_ins;
    } else {
        s2 = s1;
        for ( i = 0 ; s1 <= se && is_nucletide(s1) == 0; ++s1, ++i );        
        *ref_length = i;
        *ref = func(s2, i);
        if ( *s1 == '>' ) {
            ++s1;
            s2 = s1;
            for ( i = 0 ; s1 <= se && is_nucletide(s1) == 0; ++s1, ++i );
            *alt_length = i;
            *alt = func(s2, i);
        } else {
            return 1;
        }        

        if ( *alt_length == 1 && *ref_length == 1 ) {
            *type = var_type_snp;
        } else {
            *type = var_type_delins;
        }
    }
    return 0;
}
// Standard HGVS name should be NM_0001.2:c.123A>G; tolerant format could be NM_0001:c.123A>G (no version number);
// Check string format could be parsed and convert cds position to genome position.
// This code is fragile, improve me.
int check_hgvs_name(const char *name)
{
    if ( spec.exp == NULL ) {
        spec.exp = (regex_t*)malloc(sizeof(regex_t));
        int rv;
        rv = regcomp(spec.exp, "(.*):[gcn][.](.*)", REG_EXTENDED);
        if ( rv != 0 ) {
            error("regcomp failed with %d.", rv);
        }
    }
    // Only allow match once.
    regmatch_t matches[1]; 
    if ( regexec(spec.exp, name, 1, matches, 0) ) {
        int i, l;
        l = strlen(name);
        fprintf(stderr, "%s\n", name);
        for (i = 0; i < matches[0].rm_so; ++i ) {
            fprintf(stderr, ANSI_COLOR_RED "%c" ANSI_COLOR_RESET, '^');
        }
        for (i = matches[0].rm_so; i < matches[0].rm_eo; ++i ) {
            fprintf(stderr, ANSI_COLOR_GREEN "%c" ANSI_COLOR_RESET, '~');
        }
        for (i = matches[0].rm_eo; i < l; ++i) {
            fprintf(stderr, ANSI_COLOR_RED "%c" ANSI_COLOR_RESET, '^');
        }
        fprintf(stderr, "\n");
        return 1;
    }
    return 0;
}
int setter_description(const char *name, int _pos, char *ref, char *alt)
{
    struct hgvs_des *des = &spec.des;
    hgvs_des_clear(des);

    //func_dup_seq func = strand == '+' ? strndup : rev_seqs;
    
    const char *a = alt;
    const char *r = ref;
    // pos may differ from _pos, but _pos will not change in this function
    int pos = _pos;
    des->chrom = strdup(name);
    
    // for GATK users, <NON_REF> allele will be treat as ref
    if ( strcmp(alt, "<NON_REF>") == 0 ) {
      des->type = var_type_nonref;
      return 0;
    }
    // skip same string from start
    while (*a && *r && toupper(*a) == toupper(*r)) { a++; r++; pos++; }
    
    if ( !a[0] && !r[0] ) {
	des->type = var_type_ref;
	return 0;
    }
    
    // if ref and alternative allele are 1 base, take it as snp
    if ( *a && *r && !a[1] && !r[1] ) {
	// mpileip may output X allele, treat as ref
	if ( *a == '.' || *a == 'X' || *a == '*') {
	    des->type = var_type_ref;
	    return 0;
	}
	// for most cases, it is a snp
	des->type = var_type_snp;
	des->start = des->end = pos;
	des->ref_length = des->alt_length = 1;
	des->ref = strndup(r, 1);
	des->alt = strndup(a, 1);
	return 0;
    }
    // if alternate allele longer than ref, take it as insertion
    if ( *a && !*r ) {
	des->type = var_type_ins;
	while ( *a ) a++;
	des->start = pos;
	des->ref_length = 0;
	des->ref = str_init;
	des->end = pos +1;
	des->alt_length = (a-alt) - (r-ref);
	des->alt = strndup(a-des->alt_length, des->alt_length);
	return 0;
    }
    // if ref allele longer than alt, should be deletion
    if ( !*a && *r) {	
	des->type = var_type_del;
	while ( *r ) r++;
	des->start = pos;
	des->ref_length = (r-ref) - (a-alt);
	des->ref = strndup(r-des->ref_length, des->ref_length);
	des->end = pos + des->ref_length -1;
	des->alt_length = 0;
	des->alt = str_init;
        return 0;
    }

    // trim tails if ends are same
    const char *ae = a;
    const char *re = r;
    while ( ae[1] ) ae++;
    while ( re[1] ) re++;
    while ( re > r && ae > a && toupper(*re) == toupper(*ae) ) {
	re--;
	ae--;
    }
    if (ae == a && re == r) {
        des->type = var_type_snp;
        des->start = pos;
        des->end = pos;
        des->ref_length =  1;
        des->ref = strndup(r, 1);
        des->alt_length = 1;
        des->alt = strndup(a, 1);
        return 0;
    }
    // check a and e in first step, so "re==r && ae == a" would not happen here 
    if ( ae == a) {
	des->type = var_type_del;
	des->start = pos;
	des->ref_length = re -r + 1;
	des->ref = strndup(r, des->ref_length);
	des->end = pos + des->ref_length -1;
	des->alt_length = 0;
	des->alt = str_init;
        return 0;
    }
    if (re == r) {
	des->type = var_type_ins;
	des->start = pos;
	des->end = pos+1;
	des->ref_length = 0;
	des->ref = str_init;
	des->alt_length = ae - a+1;
	des->alt = strndup(a, des->alt_length);
        return 0;
    }
    // delins
    des->type = var_type_delins;
    des->start = pos;
    des->ref_length = re-r+1;
    des->ref = strndup(r, des->ref_length);
    des->alt_length = ae-a+1;
    des->alt = strndup(a, des->alt_length);
    des->end = pos + des->ref_length -1;
    return 0;    
}

int parse_hgvs_name(const char *name)
{
    if ( check_hgvs_name(name) )
        return 1;

    // If format test success.    
    struct hgvs_des *des = &spec.des;
    struct genepred_spec *data = spec.data;
    // Check the chromosome or transcript name.
    enum nametype type = name_is_unknown;
    char *ss = (char*)name;
    char *safe_lock = ss;
    char *se = ss;
    char *se1 = ss;
    int i;
    while ( safe_lock && *safe_lock )
        safe_lock++;

    // Keep the safe_lock always point to the last unit.
    --safe_lock;

    // Eliminate the quotes.
    while ( *ss == '"' && *safe_lock == '"' ) {
        ++ss;
        --safe_lock;
    }
    while ( *ss == '\'' && *safe_lock == '\'') {
        ++ss;
        --safe_lock;
    }
    
    for ( ; se != safe_lock && *se != ':'; ++se );    
    for ( i = 0; *se1 != '(' && se1 != se; ++se1, ++i);
    
    kstring_t string = { 0, 0, 0};
    kputsn(ss, i, &string);   
    se++;
    if ( *se == 'g' ) {
        // Assume genome locations.
        type = name_is_chromosome;
    } else if ( *se == 'c' ) {
        type = name_is_coding_rna;
    } else if ( *se == 'n' ) {
        type = name_is_noncoding_rna;
    } else {
        error("Undefined variant type. %c.", se[1]);
    }
    // Check the start position.
    for ( ; se != safe_lock && *se != '.'; ++se);
    se1 = ++se;

    // The possible character of position should be numberic and :
    // -  alias UTR5 region
    // *  alias UTR3 or intron region
    // +  alias intron
    for ( i = 0; se != safe_lock && (se[0] == '*' || se[0] == '-' || se[0] == '+' || (se[0] >= 48 && se[0] <= 57)); ++se, ++i);
    if ( se == se1 )
        error("No position found.");

    char strand = '+';
    if ( type == name_is_chromosome ) {
        des->chrom = strdup(string.s);
        // For genome locations, only number is accept.
        if ( check_num_likely_l(se1, i)) {
            des->start = str2int_l(se1, i);
        } else {
            error("Genome position not readable. %s, %d.", se1, i);
        }

        // If end position is specified.
        if ( se[0] == '_') {
            se1 = ++se;
            for ( i = 0; se[0] == '?' || se[0] == '*' || se[0] == '-' || se[0] == '+' || (se[0] >= 48 && se[0] <= 57); ++se, ++i);
            if ( se == se1 )
                error("No position found.");            
            if ( check_num_likely_l(se1, i) == 0 ) {
                des->end = str2int_l(se1, i);
            } else {
                error("Genome position not readable. %s.", se1);
            }            
        } else {
            des->end = des->start;
        }        
    } else {
        // Here only retrieve one record for each transcript. For some reasons, like alternative locus, one transcript
        // may align to several genome regions, that usually mislead researchers one gene will be expressed by different
        // genome locus in one cell which it is not true. For these multi records, I will give a warning message.
        struct genepred_line *line = genepred_retrieve_trans(data, string.s);
        if ( line == NULL )
            error("No transcript found in data. %s.", string.s);        
        parse_line_locs(line);        
        des->start = parse_position(se1, se, line);
        des->chrom = strdup(line->chrom);
        strand = line->strand;
        // Check the end position.    
        if ( se[0] == '_') {
            se1 = ++se;
            for ( i = 0; se[0] == '?' || se[0] == '*' || se[0] == '-' || se[0] == '+' || (se[0] >= 48 && se[0] <= 57); ++se, ++i);
            if ( se == se1 )
                error("No position found.");
            des->end = parse_position(se1, se, line);
            list_lite_delete(&line, genepred_line_destroy);
        } else {
            des->end = des->start;
        }

        if ( line->next != NULL ) {
            warnings("Multiple transcript hits. Only random pick one record. %s.", line->name1);
        }        
        list_lite_delete(&line, genepred_line_destroy);
    }

    for ( se1 = se; se != safe_lock && *se != '('; ++se);
    
    if ( parse_var(se1, se, (int)strand, &des->ref_length, &des->ref, &des->alt_length, &des->alt, &des->type) )
        error("Failed to parse variants. %s.", se1);

    return 0;
}
// 
static int find_the_block(struct genepred_line *line, int *blk_start, int *blk_end, int pos)
{
    *blk_start = 0;
    *blk_end = 0;
    int i;
    for ( i = 0; i < line->exon_count; ++i ) {
        int start = line->exons[BLOCK_START][i];
        int end = line->exons[BLOCK_END][i];
        if ( pos < start) {
            *blk_end = i;
            break;
        } else {
            *blk_start = i;
            if ( pos <= end ) {
                *blk_end = i;
                break;
            }           
        }
    }
    // Always return 0 ? if out of range return 1??
    return 0;
}
static int find_locate(struct genepred_line *line, int *pos, int *offset, int start, int *exon_id)
{
    int block1 = 0;
    int block2 = 0;
    if ( find_the_block(line, &block1, &block2, start ) )
        return 1;
    if ( block1 == block2 ) {
        if ( line->strand == '+' ) {
            *pos = line->loc[BLOCK_START][block1] + start - line->exons[BLOCK_START][block1];
        } else {
            *pos = line->loc[BLOCK_START][block1] - (start - line->exons[BLOCK_START][block1]);
        }
        *offset = 0;
    } else {
        int upstream =  start - line->exons[BLOCK_END][block1];
        int downstream = line->exons[BLOCK_START][block2];
        if ( upstream > downstream ) {
            *pos = line->loc[BLOCK_START][block2];
            *offset = line->strand == '+' ? -downstream : downstream;
        } else {
            *pos = line->loc[BLOCK_END][block1];
            *offset = line->strand == '+' ? upstream : -upstream;
        }        
    }
    *exon_id = block1;
    return 0;
}
static int check_func_vartype(struct genepred_line *line, int pos, int offset, int ref_length, char *ref, int alt_length, char *alt, struct var_func_type *type)
{    
    if ( line->cdsstart == line->cdsend ) {
        type->func = func_region_noncoding;
    } else {
        if ( pos <= line->utr5_length ) {
            type->func = func_region_utr5;
            // pos = line->utr5_length - pos + 1;
        } else if ( line->reference_length - pos <= line->utr3_length ) {
            type->func = func_region_utr3;
            // pos = line->utr3_length - (line->reference_length - pos);
        } else {
            type->func = func_region_cds;
            // pos = pos - line->utr5_length;
        }
    }

    struct hgvs_des *des = &spec.des;
    // Exon or intron id.
    int i;
    // CDS id.
    int h;
    int cds_pos;
    type->vartype = var_is_unknown;

#define BRANCH(_type) do {                              \
        if ( type->vartype == var_is_unknown ) {        \
            type->vartype = _type;                      \
        }                                               \
} while(0)

    if ( line->strand == '+' ) {
        for ( i = 0, h = 0; i < line->exon_count; ) {
            if ( line->loc[BLOCK_END][i] > line->utr5_length &&
                 pos > line->utr5_length &&
                 pos <= line->reference_length - line->utr3_length)
                h++;

            if ( pos >= line->loc[BLOCK_START][i] && pos <= line->loc[BLOCK_END][i] ) {
                if ( pos <= line->loc[BLOCK_START][i] + SPLIT_RANGE || pos >= line->loc[BLOCK_END][i] - SPLIT_RANGE) {
                    type->vartype = var_is_splice_site;
                }

                break;
            }
            ++i;
        }

        if ( offset == 0 ) {
            type->count = i+1;
            type->count2 = h;
        } else {
            type->count = offset < 0 ? i : i + 1;
            type->count2 = 0;
        }
    } else {
        for ( i = 0, h = 0; i < line->exon_count; ) {
            int j = line->exon_count - i - 1;
            if ( line->loc[BLOCK_START][j] > line->utr5_length &&
                 pos > line->utr5_length &&
                 pos <= line->reference_length - line->utr3_length)
                h++;

            if ( pos >= line->loc[BLOCK_END][j] && pos <= line->loc[BLOCK_START][j] ) {
                if ( pos <= line->loc[BLOCK_END][j] + SPLIT_RANGE || pos >= line->loc[BLOCK_START][j] - SPLIT_RANGE) {
                    type->vartype = var_is_splice_site;
                }                
                break;
            }
            ++i;
        }
            
        if ( offset == 0 ) {
            type->count = i+1;
            type->count2 = h;
        } else {
            type->count = offset < 0 ? i : i + 1;
            type->count2 = 0;
        }
    }
    // check if variantion located in the splice sites around the edge of utr and cds regions
    
    if ( offset != 0 ) {        
        if ( offset < SPLIT_RANGE && offset > -SPLIT_RANGE ) {
            type->vartype = var_is_splice_site;
        } else {
            type->vartype = var_is_intron;
        }
        goto no_amino_code;
    } else if ( pos > line->utr5_length - SPLIT_RANGE && pos < line->utr5_length + SPLIT_RANGE ) {
        type->vartype = var_is_splice_site;
    } else if ( line->reference_length - pos > line->utr3_length - SPLIT_RANGE && line->reference_length -pos < line->utr3_length + SPLIT_RANGE ) {
        type->vartype = var_is_splice_site;
    }
    
    
    // Check if noncoding transcript.
    if ( line->cdsstart == line->cdsend ) {
        // no cds count for noncoding transcript
        type->count2 = 0;
        if ( line->strand == '+') {
            for ( i = 0; i < line->exon_count; ++i ) {
                if ( pos >= line->loc[BLOCK_START][i] && pos <= line->loc[BLOCK_END][i]) {
                    if ( pos <= line->loc[BLOCK_START][i] + SPLIT_RANGE || pos >= line->loc[BLOCK_END][i] - SPLIT_RANGE ) {
                        type->vartype = var_is_splice_site;
                    } else {
                        type->vartype = var_is_noncoding;
                    }
                    goto no_amino_code;
                }
            }
        } else {
            for ( i = 0; i < line->exon_count; ++i ) {
                if ( pos <= line->loc[BLOCK_START][i] && pos >= line->loc[BLOCK_END][i]) {
                    if ( pos <= line->loc[BLOCK_END][i] + SPLIT_RANGE || pos >= line->loc[BLOCK_START][i] -SPLIT_RANGE ) {
                        type->vartype = var_is_splice_site;
                    } else {
                        type->vartype = var_is_noncoding;
                    }
                    goto no_amino_code;
                }
            }
        }
        type->vartype = var_is_noncoding;
        goto no_amino_code;
    }

    // Check if utr regions.
    if ( pos <= line->utr5_length -SPLIT_RANGE ) {
        BRANCH(var_is_utr5);
        goto no_amino_code;
    }

    if ( pos <= line->utr5_length ) {
        type->vartype = var_is_splice_site;
        goto no_amino_code;
    }

    cds_pos = line->reference_length - line->utr3_length;
    
    if ( pos >= cds_pos + SPLIT_RANGE ) {
        BRANCH(var_is_utr3);
        goto no_amino_code;
    }

    if ( pos >= cds_pos ) {
        type->vartype = var_is_splice_site;
        goto no_amino_code;
    }

    // amino_code:
    // Check if cds regions.
    cds_pos = pos - line->utr5_length;
    // 0 based.
    int start = (cds_pos-1)/3*3 + line->utr5_length;
    type->loc_amino = (cds_pos-1)/3 + 1;
    int cod = (cds_pos-1) % 3;
    
    int l = 0;
    char *name = line->name1;
    char *ori_seq = faidx_fetch_seq(spec.data->fai, name, start, start + 1000, &l);
    
    if ( ori_seq == NULL || l == 0 )
        goto failed_check;
    
    char *ref_seq = ref == NULL ? NULL : strdup(ref);
    char *alt_seq = alt == NULL ? NULL : strdup(alt);

    if ( line->strand == '-' ) {        
        compl_seq(ref_seq, strlen(ref_seq));
        compl_seq(alt_seq, strlen(alt_seq));
    }
    // Check the ref sequence consistant with database or not.        
    if ( ref_length > 0 && ref != NULL ) {
        int i;    
        for ( i = 0; i < ref_length; ++i ) {
            // debug_print("%d, %d, %d, %s.", line->utr5_length, start, cod, ori_seq);
            if ( seq2code4(ori_seq[cod+i]) != seq2code4(ref_seq[i]) ) {
                if ( seq2code4(ori_seq[cod+i]) == seq2code4(alt_seq[i]) ) {
                    type->vartype = var_is_no_call;
                } else {
                    warnings("Inconsistance nucletide. %s:%d, %c vs %c.", des->chrom, des->start+i,ori_seq[cod+i], ref_seq[i]);
                    break;
                }
            }
        }
    }

    if ( ref_length == alt_length ) {
        char codon[4];
        memcpy(codon, ori_seq, 3);
        codon[3] = '\0';
        type->ori_amino = codon2aminoid(codon);
        if ( ref_length == 1 ) {
            codon[cod] = *alt_seq;
        } else {
            int i;
            for ( i = cod; i < 3; ++i ) {
                codon[i] = alt_seq[i-cod];
            }                           
        }
        
        type->mut_amino = codon2aminoid(codon);
        if ( type->ori_amino == type->mut_amino ) {
            if ( type->ori_amino == 0 ) {
                BRANCH(var_is_stop_retained);
            } else {
                BRANCH(var_is_synonymous);
            }
        } else {
            if ( type->ori_amino == 0 ) {
                //type->vartype = var_is_stop_lost;
                BRANCH(var_is_stop_lost);
            } else if ( type->mut_amino == 0 ) {
                //type->vartype = var_is_nonsense;
                BRANCH(var_is_nonsense);
            } else {
                //type->vartype = var_is_missense;
                BRANCH(var_is_missense);
            }
        }                
    } else {
        // Insertion
        
        // Deletion

        // Delins
    }

    if (ref_seq != NULL ) {
        free(ref_seq);        
    }

    if ( alt_seq != NULL ) {
        free(alt_seq);
    }
        
    return 0;
    
#undef BRANCH
    
  failed_check:
    type->vartype = var_is_unknown;
    
  no_amino_code:
    type->loc_amino = 0;
    type->ori_amino = 0;
    type->mut_amino = 0;
    type->fs = 0;
    return 0;
}
int print_hgvs_summary();

// return 0 on success, 1 on out of range.
int generate_hgvs_core(struct genepred_line *line, struct hgvs_core *core, int start, int end, int ref_length, char *ref, int alt_length, char *alt)
{
    // int i;
    // int blk_start = 0, blk_end = 0;
    struct hgvs_name *name = &core->name;
    struct var_func_type *type = &core->type;
    
    if ( line->loc_parsed == 0 )
        parse_line_locs(line);

    // generate amino acid length
    name->aa_length = line->cdsstart == line->cdsend ? 0 : (line->reference_length - line->utr5_length - line->utr3_length)/3;
        
    // in case dual strands transcript RNAs, record strand for each transcript
    name->strand = line->strand;
    int exon_id1 = 0;
    int exon_id2 = 0;
    if ( find_locate(line, &name->pos, &name->offset, start, &exon_id1) )
        return 1;
    
    // Locate end. For most cases variants are snps, start == end.
    if ( end != start ) {
        find_locate(line, &name->end_pos, &name->end_offset, end, &exon_id2);
        if ( exon_id1 != exon_id2 ) {
            type->func = func_region_large;
        }
        if ( line->strand == '-') {
            int temp = name->pos;
            name->pos = name->end_pos;
            name->end_pos = temp;
            temp = name->offset;
            name->offset = name->end_offset;
            name->end_offset = temp;
        }        
    }

    name->name1 = strdup(line->name1);
    name->name2 = strdup(line->name2);

    // Check the vartype.
    if ( check_func_vartype(line, name->pos, name->offset, ref_length, ref, alt_length, alt, type) )
        return 1;

    // Update to function location.
    if ( type->func == func_region_utr5 ) {
        name->loc = line->utr5_length - name->pos + 1;
        name->end_loc = line->utr5_length - name->end_pos + 1;
    } else if ( type->func == func_region_utr3 ) {
        name->loc = line->utr3_length - ( line->reference_length - name->pos );
        name->end_loc = line->utr3_length - ( line->reference_length - name->end_pos );            
    } else if ( type->func == func_region_cds ) {
        name->loc = name->pos - line->utr5_length;
        name->end_loc = name->end_pos - line->utr5_length;
    } else {
        name->loc = name->pos;
        name->end_loc = name->end_pos;
    }
    
    return 0; 
}
// Fill all possible HGVS names for this variant.
struct hgvs_des *fill_hgvs_name()
{
    // Check the position inited.
    struct hgvs_des *des = &spec.des;
    if ( des->start == 0 )
        error("Variant position is not inited.");

    struct genepred_line *line = genepred_retrieve_region(spec.data, des->chrom, des->start-1, des->end);

    for ( ;; ) {
        if ( line == NULL )
            break;
        
        if ( des->l == des->m ) {
            des->m += 2;
            des->a = (struct hgvs_core*)realloc(des->a, des->m*sizeof(struct hgvs_core));
            memset(&des->a[des->l], 0, sizeof(struct genepred_line));
        }

        hgvs_core_clear(&des->a[des->l]);
        if ( generate_hgvs_core(line, &des->a[des->l], des->start, des->end, des->ref_length, des->ref, des->alt_length, des->alt) == 0 )
            des->l++;
        
        struct genepred_line * temp = line;
        line = line->next;
        genepred_line_destroy(temp);
    }
    // print_hgvs_summary();
    return des;
}

int print_hgvs_summary()
{
    int i;
    struct hgvs_des *des = &spec.des;
    kstring_t string = { 0, 0, 0 };
    ksprintf(&string, "\n%12s\t%10s\t%10s\t%5s\t%5s\n","#Chromosome", "Start", "End", "Ref", "Alt");
    ksprintf(&string, "%12s\t%10d\t%10d\t%5s\t%5s\n\n", des->chrom, des->start, des->end, des->ref, des->alt);
    ksprintf(&string, "%8s\t%15s\t%15s\t%15s\t%10s\t%15s\n", "#Gene","Transcript","cHGVS","pHGVS","ExInt", "VarType");
             
    for ( i = 0; i < des->l; ++i ) {
        struct hgvs_core *core = &des->a[i];
        struct var_func_type *type = &core->type;
        if ( core->name.name1 == NULL )
            error("No transcript name found.");
        ksprintf(&string, "%8s\t%15s\t", core->name.name2, core->name.name1);

        kstring_t temp = { 0, 0, 0 };
        if ( type->func == func_region_noncoding ) {
            kputs("n.", &temp);
        } else if ( type->func == func_region_cds ) {
            kputs("c.", &temp);
        } else if ( type->func == func_region_utr5 ) {
            kputs("c.-", &temp);
        } else if ( type->func == func_region_utr3 ) {
            kputs("c.*", &temp);
        }
        ksprintf(&temp,"%d", core->name.loc);
        if ( core->name.offset  > 0 ) {
            ksprintf(&temp, "+%d", core->name.offset);
        } else if ( core->name.offset < 0 ) {
            ksprintf(&temp, "%d", core->name.offset);
        }

        if ( des->start !=  des->end ) {
            kputc('_', &temp);
            if ( type->func == func_region_utr5 ) {
                kputc('-',&temp);
            } else if ( type->func == func_region_utr3 ) {
                kputc('*',&temp);
            }
            ksprintf(&temp, "%d", core->name.end_loc);
            if ( core->name.offset  > 0 ) {
                ksprintf(&temp, "+%d", core->name.offset);
            } else if ( core->name.offset < 0 ) {
                ksprintf(&temp, "%d", core->name.offset);
            }
        }
        
        if ( des->ref_length == 0 ) {
            if ( des->alt != NULL ) {
                ksprintf(&temp, "ins%s", des->alt);
            } else {
                ksprintf(&temp, "ins%d", des->alt_length);
            }
        } else if ( des->alt_length == 0 ) {
            if ( des->ref != NULL ) {
                ksprintf(&temp, "del%s", des->ref);
            } else {
                ksprintf(&temp, "del%d", des->ref_length);
            }                
        } else {
            ksprintf(&temp, "%s>%s", des->ref, des->alt);                
        }
        // cHGVS
        ksprintf(&string, "%15s\t", temp.s);
        
        // pHGVS
        temp.l = 0;
        if ( type->loc_amino > 0 ) {
            ksprintf(&temp, "p.%s%d%s", codon_names[type->ori_amino], type->loc_amino, codon_names[type->mut_amino]);
            ksprintf(&string, "%15s\t", temp.s);
        } else {
            ksprintf(&string, "%15s\t", "-");
        }
        // Exon id
        ksprintf(&string, "%10d\t", type->count);
        // VarType
        ksprintf(&string, "%15s", var_type_string(type->vartype));
        kputc('\n', &string);
        free(temp.s);
    }
    puts(string.s);
    free(string.s);
    return 0;
}

#ifdef HGVS_MAIN
int usage()
{
    fprintf(stderr,
            "Usage: hgvs_converter NM001:c.123A>T\n"
        );
    return 1;
}
int parse_args(int ac, char **av)
{
    if ( ac == 0 )
        return usage();
    const char *data_fname = 0;
    const char *fasta = 0;
    const char *name = 0;
    int i;
    for ( i =0; i < ac; ) {
        const char *a = av[i++];
        const char **var = 0;
        if ( strcmp(a, "-data") == 0  && data_fname == 0 )
            var = &data_fname;
        else if ( strcmp(a, "-fasta") == 0 && fasta == 0 )
            var = &fasta;

        if ( var != 0) {
            if (i == ac) {
                fprintf(stderr, "Missing an argument after %s.", a);
                return 1;
            }
            *var = av[i++];
            continue;
        }
        if ( name == 0 ) {
            name = a;
            continue;
        }
        fprintf(stderr, "Unknown argument : %s", a);
        return 1;
    }

    if ( data_fname == NULL )
        error("-data genepred databases is required.");

    if ( fasta == NULL )
        error("-fasta transcripts fasta file is required.");
    if ( init_hgvs_spec(data_fname, fasta) )
        return 1;
    if ( parse_hgvs_name(name) ) {
        error("Failed to parse name string.");
        return 1;
    }        
    return 0;
}

void convert_hgvs()
{
    fill_hgvs_name();
    print_hgvs_summary();
}
void release_memory()
{
}
int main(int argc, char **argv)
{
    if ( parse_args(--argc, ++argv) )
        return 1;
    convert_hgvs();
    release_memory();
    return 0;
}
#endif
