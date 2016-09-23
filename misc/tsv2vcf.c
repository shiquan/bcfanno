/* vcfconverter.c -- convert from tab-separated to VCF/BCF with extra INFO field
 * 
 *  Author:
 *      Shi Quan  shiquan.cn@gmail.com
 *
 * This program is reconstruct from Petr Danecek's tsv2vcf.c and vcfconvert.c. New functions 
 * added in vcfconverter, see man page for details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/vcf.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"

struct line {
    int n;
    int *splits;
    kstring_t string;
};

void clear_line(struct line *line)
{
    if (line->n > 0)
        free(line->splits);
    line-n = -1;
    line->string.l = 0;
}

struct tsv_col;
typedef int (*tsv_setter_t)(bcf_hdr_t *hdr, struct tsv_col *col, bcf1_t *, struct line *);

struct tsv_col {
    char *name; // tags name
    int col; // column number in databases
    int type;
    int hdr_id;
    tsv_setter_t setter;
};

char *get_col_string(struct line *line, int col)
{
    assert(col->col >= 0);
    if ( col >= line->n )
        return NULL;
    return line->string.s + line->splits[col];
}
int setter_chrom(bcf_hdr_t *hdr, struct tsv_col *col, bcf1_t *line, struct line *line)
{
    char *name = get_col_string(line, col->col);
    if (name == NULL)
        return 0;
    int id = bcf_hdr_id2int(hdr, BCF_DT_CTG, name);
    if ( id == -1 )
        return 0;
    line->rid = id;
    return 0;
}
int setter_pos(bcf_hdr_t *hdr, struct tsv_col *col, bcf1_t *line, struct line *line)
{
    char *name = get_col_string(line, col->col);
    if (name == NULL || name[0] == '.')
        return 0;
    int pos = atoi(name);
    if (pos < 0)
        return 0;
    line->pos = pos-1;
    return 0;
}
int setter_start(bcf_hdr_t *hdr, struct tsv_col *col, bcf1_t *line, struct line *line)
{
    char *name = get_col_string(line, col->col);
    if (name == NULL || name[0] == '.')
        return 0;
    int pos = atoi(name);
    if (pos < 0)
        return 0;
    line->pos = pos;    
    return 0;
}
int setter_ref(bcf_hdr_t *hdr, struct tsv_col *col, bcf1_t *line, struct line *line)
{
    
    return 0;
}
int setter_alt(bcf_hdr_t *hdr, struct tsv_col *col, bcf1_t *line, struct line *line)
{
    return 0;
}
// for INFO, return 1 on success, return 0 on failure.
int setter_info(bcf_hdr_t *hdr, struct tsv_col *col, bcf1_t *line, struct line *line)
{
    return 1;
}
int tsv_register(bcf_hdr_t *hdr, char *name, struct tsv_col *col)
{
    col->type = -1;
    col->hdr_id = -1;

    if ( strcasecmp("#chr", name) == 0 || strcasecmp("#chrom", name) == 0 ) {
        col->setter = setter_chrom;
        col->name = strdup("CHR"); // should allocate memory for name
        return 0;
    }

    if ( strcasecmp("pos", name ) == 0 ) {
        col->setter = setter_pos;
        col->name = strdup("POS"); 
        return 0;
    }

    if ( strcasecmp("start", name ) == 0 ) {
        col->setter = setter_start;
        col->name = strdup("START"); 
        return 0;
    }

    if ( strcasecmp("ref", name) == 0 ) {
        col->setter = setter_ref;
        col->name = strdup("REF"); 
        return 0;
    }

    if ( strcasecmp("alt", name) == 0 ) {
        col->setter = setter_alt;
        col->name = strdup("ALT"); 
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
        col->name = strdup("END");
        return 0;
    }

    col->hdr_id = bcf_hdr_id2int(hdr, BCF_HL_INFO, name);
    if (col->hdr_id == -1)
        return 1;
    
    col->setter = setter_info;        
    col->name = strdup(name);
    return 0;
}
struct args {
    const char *header_fname;
    const char *input_fname;
    const char *output_fname;    
    //const char *columns;
    int output_type;
    int n_cols, m_cols;
    tsv_col_t *cols;
    faidx_t *fai;
    kstring_t comment;
};

struct args args = {
    .header_fname = 0,
    .input_fname = 0,
    .output_fname = 0,    
    //.columns = 0,
    .output_type = FT_VCF,
    .n_cols = 0,
    .m_cols = 0,
    .cols = 0,
    .fai = 0,
    .comment = KSTRING_INIT,
};
int parse_line(struct line *line)
{
    if ( line->n > 0 )
        free(line->splits);
    line->splits = ksplit(&line->string, '\t', &line->n);
    if ( line->n < 3 )
        return 0;    
}

int usage(char *prog)
{
    fprintf(stderr,"Usage : %s -header|-h header.txt -O z -o out.vcf.gz in.tsv.gz\n", prog);
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
    for (i = 0; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if ((strcmp(a, "-header") == 0 || strcmp(a, "-h") == 0) && args.header_fname == 0) 
            var = &header_fname;
        else if ( strcmp(a, "-O") == 0 && out_type == 0 )
            var = &out_type;
        else if ( strcmp(a, "-o") == 0 && args.output_fname == 0 )
            var = &args.output_fname;

        if ( var != 0 ) {
            if ( i == argc)
                error("Missing argument after %s", a);
            *var = &argv[i++];
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
        return usage();
    
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
                error("The output type %s not recognised.", args.output_type);
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
    kstring str = KSTRING_INIT;
    while ( hts_getline(fp, KS_SEP_LINE, &str) > 0 ) {
        if ( strncmp(str.s,"##INFO=", 7) ) {
            warning("Unknown format, check the line %s : %s.", args.header_fname, str.s);
            continue;
        }
        if ( bcf_hdr_append(hdr, str.s) )
            error("Could not parse %s : %s.", args.header_fname, str.s);
        
        str.l = 0;
    }
    hts_close(fp);
    
    htsFile *fp_input = hts_open(args.input_fname, "r");
    kstring_t str = KSTRING_INIT;
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
    if ( comment.l == 0 )
        error("No title found in file, %s.", args.input_fname);
    
    int n;
    int *splits = ksplit(&head, '\t', &n);
    if ( strcasencmp("#chr", head.s + splits[0], 3) ||
         ( strcasecmp("pos", head.s + splits[1]) && strcasecmp("start", head.s + splits[1])) ) {
        error("Failed to parse title of %s, %s.", args.input_fname, head.s);        
    }
    for ( i = 0; i < n; ++i ) {
        if ( args.n_cols == args.m_cols ) {
            args.m_cols += 8;
            args.cols = realloc(args.cols, args.m_cols *sizeof(struct tsv_col));
        }
        if ( tsv_register(head.s+splits[i], &args.cols[args.n_cols]) )
            continue;
        args.cols[args.n_cols].col = i;
        args.n_cols++;        
    }
    
    return 0;
}
int convert_tsv_vcf()
{
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    bcf_hdr_set_chrs(hdr, args.fai);
    if (init_columns(hdr) )
        error("Empty columns.");

    htsFile *fp_output = args.output_fname == 0 ? hts_open("-", hts_bcf_wmode(out_type)) :
        hts_open(args.output_fname, hts_bcf_wmode(out_type));

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
            n += args.cols[i].setter(args.cols[i], rec, &line);
        }
        
        if ( n )
            bcf_write(fp_output, hdr, rec);
        clear_line(&line);
    }
    if ( hts_close(fp_output))
        error("Failed to close %s", args.output_fname);

    free(line.string.s);
    bcf_hdr_destroy(hdr);
    hts_close(fp_input);
    hts_close(fp_output);
    return 0;    
}

int release_memory()
{
    int i;
    for (i = 0; i < args.n; ++i )
        free(args.cols[i].name);
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
