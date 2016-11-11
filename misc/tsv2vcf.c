/*  tsv2vcf.c -- convert from tab-separated text file to VCF/BCF with extra INFO field
 * 
 *  Author:
 *      Shi Quan  shiquan.cn@gmail.com
 */

#include "utils.h"
#include <stdlib.h>
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/vcf.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"

#ifndef KSTRING_INIT
#define KSTRING_INIT { 0, 0, 0 }
#endif

static int pos_is_set = 0;
static int start_is_set = 0;
static int end_is_set = 0;

void set_start ()
{
    start_is_set = 1;
}
void set_pos ()
{
    pos_is_set = 1;
}
void set_end ()
{
    end_is_set = 1;
}


#define check_start (start_is_set == 1)
#define check_pos (pos_is_set == 1)
#define check_end (end_is_set == 1)

struct line {
    int n;
    int *splits;
    kstring_t string;
};

void clear_line(struct line *line)
{
    if (line->n > 0)
        free(line->splits);
    line->n = -1;
    line->string.l = 0;
}

struct tsv_col;
typedef int (*tsv_setter_t)( bcf_hdr_t *hdr, struct tsv_col *col, bcf1_t *, struct line *);

struct tsv_col {
    char *key; // tags name
    int col; // column number in databases
    int type;
    int hdr_id;
    tsv_setter_t setter;
};

char *get_col_string(struct line *line, int col)
{
    if (col == -1 )
        return NULL;
    assert(col>=0);
    if ( col >= line->n )
        return NULL;
    return line->string.s + line->splits[col];
}
struct ref_alt_spec {
    int chrom_col;
    int pos_col;
    int end_col; 
    int ref_col;
    int alt_col;
    kstring_t string;
};
struct args {
    const char *header_fname;
    const char *input_fname;
    const char *output_fname;
    const char *reference_fname;
    int pos_column;
    int start_column;
    int end_column;
    int chr_column;
    int force;
    // char *columns;
    int output_type;
    int n_cols, m_cols;
    struct tsv_col *cols;
    faidx_t *fai;
    kstring_t comment;
    struct ref_alt_spec alleles;
};

struct args args = {
    .header_fname = 0,
    .input_fname = 0,
    .output_fname = 0,
    .reference_fname = 0,
    .pos_column = -1,
    .start_column = -1,
    .end_column = -1,
    .chr_column = -1,
    .force = 0,
    .output_type = FT_VCF,
    .n_cols = 0,
    .m_cols = 0,
    .cols = 0,
    .fai = 0,
    .comment = KSTRING_INIT,
    .alleles = { -1, -1, -1, -1, -1, KSTRING_INIT },    
};

int seq2num(char c)
{
    int i;
    switch(c) {
        case 'a':
        case 'A':
            i = 0;
            break;
        case 'c':
        case 'C':
            i = 1;
            break;
        case 'g':
        case 'G':
            i = 2;
            break;
        case 't':
        case 'T':
            i = 3;
            break;

        default:
            i = 4;
            break;
    }
    return i;
}
// 1 on success, 0 on failure
int check_is_number(char *string)
{
    char *ss = string;

    for ( ;ss && *ss && *ss != '\0'; ss++ ) {
        if ( isspace(*ss) )
            continue;
        if ( !isdigit((int)*ss)) {
            warnings("%s is looks not like a number.", string);
            return 0;
        }
    }
    return 1;
}
void clear_spec(struct ref_alt_spec *spec)
{
    spec->string.l = 0;
}
// construct chrom, pos and ref, alt (optional) for bcf1_t *rec
// 1 on failure, 0 on success
int construct_basic_inf(faidx_t *fai, struct ref_alt_spec *spec, bcf_hdr_t *hdr, struct line *line, bcf1_t *rec)
{
    // todo : convert region only
    if ( spec->ref_col == -1 && spec->alt_col == -1 ) {
        warnings("no ref and alt columns specified.");
        return 1;
    }

    const char seqs[5] = "ACGTN";

    char *chrom = get_col_string(line, spec->chrom_col);
    if (chrom == NULL)
        error("no chrom column found.");
    int id = bcf_hdr_id2int(hdr, BCF_DT_CTG, chrom);
    if ( id == -1 ) {
        warnings("Chromosome %s not found in reference.", chrom);
        return 1;
    }
    rec->rid = id;
    clear_spec(spec);
    char *start_s = get_col_string(line, spec->pos_col);
    if ( check_is_number(start_s) == 0 ) {
        rec->pos = -1;
        return 1;
    }
    
    int start = atoi(start_s);
    if ( pos_is_set )
        start --;

    assert(start >= 0);
    int end = 0;    
    if ( end_is_set ) {
        assert(spec->end_col >= 0);
        char *end_s = get_col_string(line, spec->end_col);
        if ( check_is_number(end_s) ) {
            end = atoi(end_s);           
        }
    }    
    char *alt = NULL;
    char *ref = NULL;
    int n;
    char *seq = NULL;
    int is_capped = 1; // assume the indel variants are capped to left base
    int ref_length, alt_length;
    int i = 0;
    
    if ( spec->ref_col == -1 ) {
        // if end column is not set, assume one position, this is very risk, todo: force check the length of reference
        ref_length = end_is_set ? end - start : 1;
        
        if ( spec->alt_col == -1 ) { // assume only convert regions into vcf
            if ( end_is_set ) {
                int id = bcf_hdr_id2int(hdr, BCF_DT_ID, "END");
                if ( id == -1 ) {
                    bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position.\">");
                    bcf_hdr_sync(hdr);
                    id = bcf_hdr_id2int(hdr, BCF_DT_ID, "END");
                    assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));
                }
                bcf_update_info_int32(hdr, rec, "END", &end, 1);
            }
        } else {
            ref = faidx_fetch_seq(fai, chrom, start, start+ref_length-1, &n);
            alt = get_col_string(line, spec->alt_col);
            alt_length = strlen(alt);
            for ( i = 0; i < ref_length; i++)
                kputc(seqs[seq2num(ref[i])], &spec->string);
            free(ref);
            if (alt_length == 1 && *alt == '.')
                goto update_alleles;
            kputc(',', &spec->string);
            for ( i = 0; i < alt_length; i++)
                kputc(seqs[seq2num(alt[i])], &spec->string);
        }
        goto update_alleles;
    } 

    ref = get_col_string(line, spec->ref_col);
    if (ref == NULL || *ref == '.'  || *ref == '\0' || strcmp(ref, "(null)") == 0) {
        is_capped = 0;
        ref_length = 0;
    } else {
        ref_length = strlen(ref);                    
    }

    // if end position specified, check the length of reference consistent with region length
    if ( is_capped && end_is_set && ref_length != end - start )
        error("Unconsistent reference length : %s\t%d\t%d vs. %d", chrom, start, end, ref_length);
    // if nocapped, left offset 1 base in the genome corrdinate
    if ( is_capped == 0) {
        --start;
        ref_length = 1;
    }
    seq = faidx_fetch_seq(fai, chrom, start, start+ref_length-1, &n);    
    int strand = 0; // 0 for plus, 1 for minus        
    if ( ref_length == 1) {
        // assume plus strand
        if ( is_capped == 0 ) {
            kputc(seqs[seq2num(seq[0])], &spec->string);
        } else {
            if ( seq2num(ref[0]) == seq2num(seq[0]) ) {
                kputc(seqs[seq2num(ref[0])], &spec->string);
            } else if ( seq2num(ref[0]) + seq2num(seq[ref_length-i-1]) == 3 ) {
                kputc(seqs[3-seq2num(ref[0])], &spec->string);
                strand = 1;
            } else {
                if ( args.force ) {
                    warnings("bad seq at %s:%d %s vs %s", chrom, start+1, ref, seq);
                    kputc(seqs[seq2num(ref[0])], &spec->string);
                } else {
                    error("bad seq at %s:%d %s vs %s", chrom, start+1, ref, seq);
                }
            }
        }
    } else {
        for ( i = 0; i < ref_length;  i++) {
            if ( seq2num(ref[i]) == 4 )  {
                if ( args.force ) {
                    warnings("bad seq at %s:%d %s vs %s", chrom, start+1, ref, seq);                    
                } else {
                    error("bad seq at %s:%d %s vs %s", chrom, start+1, ref, seq);
                }
            }
            if ( seq2num(ref[i]) != seq2num(seq[i]) ) {
                if (seq2num(ref[ref_length-i-1]) + seq2num(seq[i]) == 3) {
                    strand = 1;
                } else {
                    if ( args.force ) {
                        warnings("bad seq at %s:%d %s vs %s", chrom, start+1, ref, seq);                        
                    } else {
                        error("bad seq at %s:%d %s vs %s", chrom, start+1, ref, seq);
                    }
                }
            }                 
        }
        if ( strand ) {
            for ( i = 0; i < ref_length; i++) 
                kputc(seqs[seq2num(ref[ref_length-i-1])], &spec->string);
        } else {
            for ( i = 0; i < ref_length; i++)
                kputc(seqs[seq2num(ref[i])], &spec->string);     
        }        
    }

    alt = get_col_string(line, spec->alt_col);
    if ( alt == NULL || (alt[0] == '.' && alt[1] == 0) || strcmp(alt, "(null)") == 0) {
        // kputs(",.", &spec->string);
        goto update_alleles;
    }
    alt_length = strlen(alt);
    kputc(',', &spec->string);
    if ( is_capped == 0)
        kputc(seqs[seq2num(seq[0])], &spec->string);

    if ( strand ) {
        for ( i = 0; i < alt_length; ++i )
            kputc(seqs[seq2num(alt[alt_length-i-1])], &spec->string);
    } else {
        for ( i = 0; i < alt_length; ++i )
            kputc(seqs[seq2num(alt[i])], &spec->string);
    }

  update_alleles:
    if (seq)
        free(seq);
    // debug_print("%s\n", spec->string.s);
    rec->pos = start;
    bcf_update_alleles_str(hdr, rec, spec->string.s);
    return 0;
}

int setter_chrom( bcf_hdr_t *hdr, struct tsv_col *col, bcf1_t *rec, struct line *line)
{
    char *name = get_col_string(line, col->col);
    if (name == NULL)
        return 0;
    if (name[0] == '.')
        return 0;

    return 0;
}
/* int setter_pos( bcf_hdr_t *hdr, struct tsv_col *col, bcf1_t *rec, struct line *line) */
/* { */
/*     char *name = get_col_string(line, col->col); */
/*     if (name == NULL || (int)name[0] == '.') */
/*         return 0; */
/*     if (name[0] == '.') */
/*         return 0; */
/*     int pos = atoi(name); */
/*     if (pos < 0) */
/*         return 0; */
/*     rec->pos = pos-1; */
/*     return 0; */
/* } */
/* int setter_start( bcf_hdr_t *hdr, struct tsv_col *col, bcf1_t *rec, struct line *line) */
/* { */
/*     char *name = get_col_string(line, col->col); */
/*     if (name == NULL || (int)name[0] == '.') */
/*         return 0; */
/*     int pos = atoi(name); */
/*     if (pos < 0) */
/*         return 0; */
/*     rec->pos = pos; */
/*     return 0; */
/* } */

int vcf_setter_alleles( bcf_hdr_t *hdr, bcf1_t *rec,  char *allele_string)
{
    bcf_update_alleles_str(hdr, rec, allele_string);
    return 0;    
}

void *split_string(char *string, int *n, int type)
{
    if ( string == NULL || (string[0] == '.' && string[1] == 0)) {
        *n = 0;
        return NULL;
    }
    
    kstring_t tmp = KSTRING_INIT;
    kputs(string, &tmp);
    int *splits = ksplit(&tmp, ',', n);       
    int i;
    if ( type == BCF_HT_INT ) {
        int *a = (int *)calloc(*n, sizeof(int));
        for ( i = 0; i < *n; ++i )
            a[i] = atoi(tmp.s+splits[i]);
        free(tmp.s);
        free(splits);
        return (void*)a;
    }
    
    if ( type == BCF_HT_REAL ) {
        float *a = (float*)calloc(*n, sizeof(float));
        for ( i = 0; i < *n; ++i )
            a[i] = atof(tmp.s+splits[i]);
        free(tmp.s);
        free(splits);
        return (void*)a;    
    }
    
    if ( type == BCF_HT_STR || type == BCF_HT_FLAG) {
        //char **s = (char**)calloc(*n, sizeof(char*));
        //for ( i = 0; i < *n; ++i )
        //   s[i] = strdup(tmp.s+splits[i]);
        char *s = strdup(tmp.s);
        char *ss = s;
        while (*ss) {
            if (*ss == ';') {*ss = '|';}
            ss++;
        }
        free(tmp.s);
        free(splits);
        return (void*)s;
    }
    free(tmp.s);
    if ( splits )
        free(splits);
    *n = 0;
    error("Unknown type %s : %d.", string, type);
    return NULL;
}
// for INFO, return 1 on success, return 0 on failure.
int setter_info( bcf_hdr_t *hdr, struct tsv_col *col, bcf1_t *rec, struct line *line)
{
    char *name = get_col_string(line, col->col);
    if (name == NULL || ((int)name[0] == '.' && (int)name[1] == 0 ))
        return 0;

    int n = 0;
    void *value = split_string(name, &n, col->type);

    if ( n == 0 )
        return 0;
    bcf_update_info(hdr, rec, col->key, value, n, col->type);
    free(value);
    return 1;
}
int tsv_register( bcf_hdr_t *hdr, char *name, struct tsv_col *col)
{
    col->type = -1;
    col->hdr_id = -1;
    while ( *name == '#' )
        name++;
    
    if ( strcasecmp("chr", name) == 0 || strcasecmp("chrom", name) == 0) {
        col->setter = setter_chrom;
        col->key = strdup("CHR"); // should allocate memory for name
        return 0;
    }

    /* if ( strcasecmp("pos", name ) == 0 ) { */
    /*     col->setter = setter_pos; */
    /*     col->key = strdup("POS"); */
    /*     return 0; */
    /* } */

    /* if ( strcasecmp("start", name ) == 0 ) { */
    /*     col->setter = setter_start; */
    /*     col->key = strdup("START"); */
    /*     return 0; */
    /* } */

    if ( strcasecmp("end", name) == 0 ) {
        // append END tag in the header, this is mandontary
        col->hdr_id = bcf_hdr_id2int(hdr, BCF_HL_INFO, "END");
        if ( col->hdr_id == -1 ) {
            bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">");
            bcf_hdr_sync(hdr);
            col->hdr_id = bcf_hdr_id2int(hdr, BCF_HL_INFO, "END");
        }
        col->type = bcf_hdr_id2type(hdr, BCF_HL_INFO, col->hdr_id);
        col->key = strdup("END");
        return 0;
    }
    
    col->hdr_id = bcf_hdr_id2int(hdr, BCF_DT_ID, name);
    if (col->hdr_id == -1)
        return 1;
    col->type = bcf_hdr_id2type(hdr, BCF_DT_ID, col->hdr_id);
    col->setter = setter_info;        
    col->key = strdup(name);
    return 0;
}
static const char *hts_bcf_wmode(int file_type)
{
    if ( file_type == FT_BCF ) return "wbu";    // uncompressed BCF
    if ( file_type & FT_BCF ) return "wb";      // compressed BCF
    if ( file_type & FT_GZ ) return "wz";       // compressed VCF
    return "w";                                 // uncompressed VCF
}
static void bcf_hdr_set_chrs(bcf_hdr_t *hdr, faidx_t *fai)
{
    int i, n = faidx_nseq(fai);
    for (i=0; i<n; i++)
    {
        const char *seq = faidx_iseq(fai,i);
        int len = faidx_seq_len(fai, seq);
        bcf_hdr_printf(hdr, "##contig=<ID=%s,length=%d>", seq,len);
    }
}

int parse_line(struct line *line)
{
    if ( line->n > 0 )
        free(line->splits);
    line->splits = ksplit(&line->string, '\t', &line->n);
    if ( line->n < 3 )
        return 0;
    return 1;
}

int usage(char *prog)
{
    fprintf(stderr, "Usage : %s -header|-h header.txt -r reference.fa [-force -pos column -O z -o out.vcf.gz] in.tsv.gz\n", prog);
    fprintf(stderr, "        -header, -h     header file\n");
    fprintf(stderr, "        -r              reference file\n");
    fprintf(stderr, "        -pos            position column, if set will skip pos,start,end column in the title\n");
    fprintf(stderr, "        -start          start position column, inconsistance with -pos\n");
    fprintf(stderr, "        -end            end position column, only set if need add a END tag in the INFO\n");
    fprintf(stderr, "        -chr            chr column, if set will skip first column in the title\n");
    fprintf(stderr, "        -force          if reference seq and fasta file are inconsistent, just give a warning\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Homepage: https://github.com/shiquan/vcfanno\n");
    return 1;
}

int parse_args(int argc, char **argv)
{
    if ( argc == 1 )
        return usage(argv[0]);

    int i;
    for (i = 0; i < argc; ++i)
        kputs(argv[i], &args.comment);
    
    const char *out_type = 0;
    const char *pos_column = 0;
    const char *chr_column = 0;
    const char *start_column = 0;
    const char *end_column = 0;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        if ((strcmp(a, "-header") == 0 || strcmp(a, "-h") == 0) && args.header_fname == 0) 
            var = &args.header_fname;
        else if ( strcmp(a, "-O") == 0 && out_type == 0 )
            var = &out_type;
        else if ( strcmp(a, "-o") == 0 && args.output_fname == 0 )
            var = &args.output_fname;
        else if ( strcmp(a, "-r") == 0 && args.reference_fname == 0)
            var = &args.reference_fname;
        else if ( strcmp(a, "-pos") == 0 && args.pos_column == -1)
            var = &pos_column;
        else if ( strcmp(a, "-chr") == 0 && args.chr_column == -1)
            var = &chr_column;
        else if ( strcmp(a, "-start") == 0 && args.start_column == -1)
            var = &start_column;
        else if ( strcmp(a, "-end") == 0 && args.end_column == -1)
            var = &end_column;
       
        if ( var != 0 ) {
            if ( i == argc)
                error("Missing argument after %s", a);
            *var = argv[i++];
            continue;
        }

        if ( strcmp(a, "-force") == 0 ) {
            args.force = 1;
            continue;
        }
        
        if ( args.input_fname == 0) {
            args.input_fname = a;
            continue;
        }

        error("Unknown argument %s.", a);        
    }
    if ( args.header_fname == 0 )
        error("No header file; use -header / -h to specify a header.txt file.");
    
    if ( args.input_fname == 0 && !isatty(fileno(stdin)) )
        args.input_fname = "-";
    
    if ( args.input_fname == 0 )
        return usage(argv[0]);
    if ( args.reference_fname == 0)
        return usage(argv[0]);
    args.fai = fai_load(args.reference_fname);
    if ( out_type != 0 ) {
        switch (out_type[0]) {
            case 'b':
                args.output_type = FT_BCF_GZ;
                break;
            case 'u':
                args.output_type = FT_BCF;
                break;
            case 'z':
                args.output_type = FT_VCF_GZ;
                break;
            case 'v':
                args.output_type = FT_VCF;
                break;
            default:
                error("The output type %d not recognised.", args.output_type);
        }
    }
    if ( pos_column ) {
        args.pos_column = atoi(pos_column);
        set_pos();
    }
    if ( chr_column ) 
        args.chr_column = atoi(chr_column);
        
    if ( start_column ) {
        args.start_column = atoi(start_column);
        set_start();
    }

    if ( end_column ) {
        args.end_column = atoi(end_column);
        set_end();
    }

    if ( args.chr_column == -1)
        error("No chrom column set, please use -chr specify it.");
    
    if ( check_end )
        args.end_column = atoi(end_column);
    if ( pos_is_set == 0 && start_is_set == 0)
        error("No position column is set.");
    /* if ( start_is_set == 1 && end_is_set == 0) */
    /*     warnings("End column is not set, treat start position zero based."); */
    if ( pos_is_set == 1 && start_is_set == 1) {
        warnings("Redundancy columns, pos column and start column both set. Skip pos column..");
        pos_is_set = 0;        
    }
    if ( pos_is_set == 1 && end_is_set == 1)  {
        warnings("error Format; pos and end position both set; assuming pos is 1 based, but standard bed format required 0 based start position. Skip end position ..");
        end_is_set = 0;
    }
   
    return 0;
}

int init_columns(bcf_hdr_t *hdr)
{    
    // init header.txt
    htsFile *fp = hts_open(args.header_fname, "rb");
    if ( fp == 0)
        error("Failed to open %s.",args.header_fname);
    kstring_t str = KSTRING_INIT;
    while ( hts_getline(fp, KS_SEP_LINE, &str) > 0 ) {
        if ( strncmp(str.s,"##INFO=", 7) ) {
            warnings("Unknown format, check the line %s : %s.", args.header_fname, str.s);
            continue;
        }
        if ( bcf_hdr_append(hdr, str.s) )
            error("Could not parse %s : %s.", args.header_fname, str.s);        
        str.l = 0;
    } 
    bcf_hdr_sync(hdr);
    hts_close(fp);
    
    htsFile *fp_input = hts_open(args.input_fname, "r");

    kstring_t head = KSTRING_INIT;
    do {
        if ( hts_getline(fp_input, KS_SEP_LINE, &str) == 0 )
            break;
        if (str.s[0] == '#') {
            head.l = 0;
            kputs(str.s, &head);            
        } else {
            break;
        }
        str.l = 0;
    } while(1);
    free(str.s);
    hts_close(fp_input);
    if ( head.l == 0 )
        error("No title found in file, %s.", args.input_fname);
    int i;
    int n;
    int *splits = ksplit(&head, '\t', &n);
    if ( strcasecmp("#chr", head.s + splits[0]) && strcasecmp("#chrom", head.s + splits[0]) )  {
        error("Failed to parse title of %s, %s.", args.input_fname, head.s);        
    }

    for ( i = 0; i < n; ++i ) {
        if ( strcasecmp("ref", head.s+splits[i]) == 0 ) {
            args.alleles.ref_col = i;
            continue;
        }

        if ( strcasecmp("alt", head.s+splits[i]) == 0 ) {
            args.alleles.alt_col = i;
            continue;
        }
        if ( args.n_cols == args.m_cols ) {
            args.m_cols += 8;
            args.cols = (struct tsv_col *)realloc(args.cols, args.m_cols *sizeof(struct tsv_col));
        }
        if (args.chr_column == i + 1) {
            args.alleles.chrom_col = i;
            tsv_register(hdr, "CHR", &args.cols[args.n_cols]);
            goto col_increase;
        } 
        
        if (args.pos_column == i + 1) {
            args.alleles.pos_col = i;
            // tsv_register(hdr, "POS", &args.cols[args.n_cols]);
            continue;
        }

        if (args.start_column == i+1) {
            args.alleles.pos_col = i;
            // tsv_register(hdr, "start", &args.cols[args.n_cols]);
            continue;
        }
        
        if (args.end_column == i+1) {
            args.alleles.end_col = i;
            tsv_register(hdr, "end", &args.cols[args.n_cols]);
            goto col_increase;            
        }

        if ( tsv_register(hdr, head.s+splits[i], &args.cols[args.n_cols]) )
            continue;
      col_increase:
        args.cols[args.n_cols].col = i;
        args.n_cols++;        
    }

    if ( args.alleles.ref_col == -1 ) {
        warnings("No ref column.");
    }

    if ( args.alleles.alt_col == -1 ) {
        warnings("No alt column.");        
    }
            
    free(head.s);
    free(splits);
    return 0;
}
int convert_tsv_vcf()
{
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    bcf_hdr_set_chrs(hdr, args.fai);
    if (init_columns(hdr) )
        error("Empty columns.");

    htsFile *fp_output = args.output_fname == 0 ? hts_open("-", hts_bcf_wmode(args.output_type)) :
        hts_open(args.output_fname, hts_bcf_wmode(args.output_type));

    bcf_hdr_write(fp_output, hdr);
    bcf1_t *rec = bcf_init();
    bcf_float_set_missing(rec->qual);
    struct line line = {
        .n = -1,
        .splits = 0,
        .string = KSTRING_INIT,
    };
    htsFile *fp_input = hts_open(args.input_fname, "r");
    int i;
    int n;
    while ( hts_getline(fp_input, KS_SEP_LINE, &line.string) > 0 ) {
        if ( line.string.s[0] == '#')
            continue;
        if ( parse_line(&line) == 0)
            continue;
        bcf_clear(rec);
        construct_basic_inf(args.fai, &args.alleles, hdr, &line, rec);
        n = 0;

        for ( i = 0; i < args.n_cols; ++i )
            n += args.cols[i].setter(hdr, &args.cols[i], rec, &line);
        
        clear_line(&line);
        if (rec->rid == 0 && rec->pos == 0)
            continue;
        if ( n )
            bcf_write(fp_output, hdr, rec);
    }

    if ( hts_close(fp_output))
        error("Failed to close %s", args.output_fname);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp_input);
    return 0;    
}

int release_memory()
{
    int i;
    for (i = 0; i < args.n_cols; ++i )
        free(args.cols[i].key);
    free(args.alleles.string.s);
    free(args.cols);
    fai_destroy(args.fai);
    free(args.comment.s);
    return 0;
}

int main(int argc, char **argv)
{    
    if ( parse_args(argc, argv) )
        return 1;
    
    convert_tsv_vcf();

    release_memory();

    return 0;
}
