/* vcfconverter.c -- convert from tab-separated to VCF/BCF with extra INFO field
 * 
 *  Author:
 *      Shi Quan  shiquan.cn@gmail.com
 *
 * This program is reconstruct from Petr Danecek's tsv2vcf.c and vcfconvert.c. New functions 
 * added in vcfconverter, see man page for details.
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
    assert(col >= 0);
    if ( col >= line->n )
        return NULL;
    return line->string.s + line->splits[col];
}
struct ref_alt_spec {
    int ref_col;
    int alt_col;
    kstring_t string;
};

void contruct_alleles(faidx_t *fai, struct ref_alt_spec *spec, struct line *line, const char *chrom, int pos)
{
    static int seq2num[256] = {
	4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 
	4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 
	4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 
	0 , 1 , 2 , 3 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 
	4 , 0 , 4 , 1 , 4 , 4 , 4 , 2 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 
	4 , 4 , 4 , 4 , 3 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 
	4 , 0 , 4 , 1 , 4 , 4 , 4 , 2 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 
	4 , 4 , 4 , 4 , 3 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 
	4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 
	4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 
	4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 
	4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 
	4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 
	4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 
	4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 
	4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 4
    };
    const char seqs[5] = "ACGTN";
    
    spec->string.l = 0;
    
    char *name  = spec->ref_col == -1 ? "N" : get_col_string(line, spec->ref_col);
    int length = strlen(name);
    int n = 0;
    char *seq = faidx_fetch_seq(fai, chrom, pos, pos+length-1, &n);
    int i = 0;
    int strand = 0; // 0 for plus, 1 for minus
    if ( length == 1) {
        if (seq2num[(int)name[0]] == 4 ) {
            kputc('N', &spec->string); // assume plus strand
        } else {
            if ( seq2num[(int)name[0]] == seq2num[(int)seq[0]] ) {
                kputc(seqs[seq2num[(int)name[0]]], &spec->string);
            } else if ( seq2num[(int)name[0]] + seq2num[(int)seq[length-i-1]] == 3 ) {
                kputc(seqs[3-seq2num[(int)name[0]]], &spec->string);
                strand = 1;
            }
        }
    } else {
        for ( i = 0; i < length;  i++) {
            if ( seq2num[(int)name[i]] == 4 )
                error("bad seq : %s vs %s", name, seq);

            if ( seq2num[(int)name[i]] != seq2num[(int)seq[i]] ) {
                if (seq2num[(int)name[length-i-1]] + seq2num[(int)seq[i]] == 3) {
                    strand = 1;
                } else {
                    error("bad seq : %s vs %s", name, seq);
                }
            }                 
        }
        if ( strand ) {
            for ( i = 0; i < length; i++) 
                kputc((int)seq[seq2num[(int)name[length-i-1]]], &spec->string);
        } else {
            for ( i = 0; i < length; i++)
                kputc((int)seq[seq2num[(int)name[i]]], &spec->string);            
        }        
    }

    name = get_col_string(line, spec->alt_col);
    if ( name ) {
        length = strlen(name);
        kputc(',', &spec->string);
        if ( strand ) {
            for ( i = 0; i < length; ++i )
                kputc((int)seq[seq2num[(int)name[length-i-1]]], &spec->string);
        } else {
            for ( i = 0; i < length; ++i )
                kputc((int)seq[seq2num[(int)name[i]]], &spec->string);
        }
    }
}

int setter_chrom( bcf_hdr_t *hdr, struct tsv_col *col, bcf1_t *rec, struct line *line)
{
    char *name = get_col_string(line, col->col);
    if (name == NULL)
        return 0;
    int id = bcf_hdr_id2int(hdr, BCF_DT_CTG, name);
    if ( id == -1 )
        return 0;
    rec->rid = id;
    return 0;
}
int setter_pos( bcf_hdr_t *hdr, struct tsv_col *col, bcf1_t *rec, struct line *line)
{
    char *name = get_col_string(line, col->col);
    if (name == NULL || (int)name[0] == '.')
        return 0;
    int pos = atoi(name);
    if (pos < 0)
        return 0;
    rec->pos = pos-1;
    return 0;
}
int setter_start( bcf_hdr_t *hdr, struct tsv_col *col, bcf1_t *rec, struct line *line)
{
    char *name = get_col_string(line, col->col);
    if (name == NULL || (int)name[0] == '.')
        return 0;
    int pos = atoi(name);
    if (pos < 0)
        return 0;
    rec->pos = pos;    
    return 0;
}

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
        return (void*)a;
    }
    
    if ( type == BCF_HT_REAL ) {
        float *a = (float*)calloc(*n, sizeof(float));
        for ( i = 0; i < *n; ++i )
            a[i] = atoi(tmp.s+splits[i]);
        return (void*)a;    
    }
    
    if ( type == BCF_HT_STR || type == BCF_HT_FLAG) {
        //char **s = (char**)calloc(*n, sizeof(char*));
        //for ( i = 0; i < *n; ++i )
        //   s[i] = strdup(tmp.s+splits[i]);
        char *s = strdup(tmp.s);
        return (void*)s;
    }
    
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
    return 1;
}
int tsv_register( bcf_hdr_t *hdr, char *name, struct tsv_col *col)
{
    col->type = -1;
    col->hdr_id = -1;

    if ( strcasecmp("#chr", name) == 0 || strcasecmp("#chrom", name) == 0 ) {
        col->setter = setter_chrom;
        col->key = strdup("CHR"); // should allocate memory for name
        return 0;
    }

    if ( strcasecmp("pos", name ) == 0 ) {
        col->setter = setter_pos;
        col->key = strdup("POS"); 
        return 0;
    }

    if ( strcasecmp("start", name ) == 0 ) {
        col->setter = setter_start;
        col->key = strdup("START"); 
        return 0;
    }

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

struct args {
    const char *header_fname;
    const char *input_fname;
    const char *output_fname;
    const char *reference_fname;
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
    //.columns = 0,
    .output_type = FT_VCF,
    .n_cols = 0,
    .m_cols = 0,
    .cols = 0,
    .fai = 0,
    .comment = KSTRING_INIT,
    .alleles = { -1, -1, KSTRING_INIT },
};
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
    fprintf(stderr,"Usage : %s -header|-h header.txt -r reference.fa -O z -o out.vcf.gz in.tsv.gz\n", prog);
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
        
        if ( var != 0 ) {
            if ( i == argc)
                error("Missing argument after %s", a);
            *var = argv[i++];
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
    if ( (strcasecmp("#chr", head.s + splits[0]) && strcasecmp("#chrom", head.s + splits[0])) ||
         ( strcasecmp("pos", head.s + splits[1]) && strcasecmp("start", head.s + splits[1])) ) {
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
        
        if ( tsv_register(hdr, head.s+splits[i], &args.cols[args.n_cols]) )
            continue;
        args.cols[args.n_cols].col = i;
        args.n_cols++;        
    }

    if ( args.alleles.ref_col == -1 ) {
        warnings("No ref column.");
    }

    if ( args.alleles.alt_col == -1 ) {
        warnings("No alt column.");        
    }
    
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
        n = 0;
        for ( i = 0; i < args.n_cols; ++i ) {
            n += args.cols[i].setter(hdr, &args.cols[i], rec, &line);
        }
        contruct_alleles(args.fai, &args.alleles, &line, bcf_seqname(hdr, rec), rec->pos+1);
        vcf_setter_alleles(hdr, rec, args.alleles.string.s);
        if ( n )
            bcf_write(fp_output, hdr, rec);
        clear_line(&line);
    }
    if ( hts_close(fp_output))
        error("Failed to close %s", args.output_fname);

    free(line.string.s);
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
