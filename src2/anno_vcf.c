#include "utils.h"
#include "anno_vcf.h"
#include "anno_col.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/bgzf.h"

static int match_allele(bcf1_t *line, bcf1_t *dat)
{
    int line_type = bcf_get_variant_types(line);
    int dat_type = bcf_get_variant_types(dat);

    if ( (line_type & dat_type) == 0 )
        return 1;
    int i, j;
    for ( i = 1; i < line->n_allele; i++) {
        for ( j = 1; j < dat->n_allele; j++) {
            if ( strcmp(line->d.allele[i], dat->d.allele[j]) == 0 )
                return 0;
        }
    }
    // if no matchs
    return 1;
}
// fill_buffer update returns
// return -1 on no change
//         0 on empty
//         else number of records
static int anno_vcf_update_buffer(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line)
{
    struct anno_vcf_buffer *b = f->buffer;
    if ( b->cached && b->buffer[0]->rid == line->rid && b->buffer[b->cached-1]->pos >= line->pos && b->buffer[0]->pos <= line->pos)
        return -1;
    else b->cached = 0;
    if ( b->last_rid != line->rid ) {
        b->no_such_chrom = 0;
        b->last_rid = line->rid;
    }
    else if ( b->no_such_chrom == 1 )
        return 0;

    if ( f->itr ) {
        hts_itr_destroy(f->itr);
        f->itr = NULL;
    }
    int l, i;
    for ( i = 1; i < line->n_allele; i++ )
        if ( l > line->d.var[i].n )
            l = line->d.var[i].n;
    int end_pos = l < 0 ? line->pos - l : line->pos;

    if ( f->tbx_idx ) {
        int tid = tbx_name2id(f->tbx_idx, bcf_seqname(hdr, line));
        if ( tid == -1 ) {
            if ( b->no_such_chrom == 0 ) {
                warnings("No chromosome %s found in %s.", bcf_seqname(hdr, line), f->fname);
                b->no_such_chrom = 1;
            }
            return 0;
        }
        f->itr = tbx_itr_queryi(f->tbx_idx, tid, line->pos, end_pos+1);
    }
    else if ( f->bcf_idx ) {
        // check id in header of database
        int tid = bcf_hdr_name2id(f->hdr, bcf_seqname(hdr, line));
        if ( tid == -1 ) {
            if ( b->no_such_chrom == 0 )  {
                warnings("No chromsome %s found in %s.", bcf_seqname(hdr, line), f->fname);
                b->no_such_chrom = 1;
            }
            return 0;
        }
        f->itr = bcf_itr_queryi(f->bcf_idx, tid, line->pos, end_pos+1);
    }
    else goto load_index_failed;

    // no record
    if ( f->itr == NULL )
        return 0;

    for ( ;; ) {
        if ( b->cached == b->max ) {
            b->max += 8;
            b->buffer = realloc(b->buffer, sizeof(void*)*b->max);
            int i;
            for ( i = 8; i > 0; --i) {
                b->buffer[b->max-i] = bcf_init();
                b->buffer[b->max-i]->pos = -1;
            }
        }
        if ( f->tbx_idx ) {
            kstring_t str = {0,0,0};
            if ( tbx_itr_next(f->fp, f->tbx_idx, f->itr, &str) < 0 ) {
                if ( str.m ) free(str.s);
                break;
            }
            vcf_parse1(&str, f->hdr, b->buffer[b->cached]);
            free(str.s);
        }
        else if ( f->bcf_idx ) {
            if ( bcf_itr_next(f->fp, f->itr, b->buffer[b->cached]) < 0 )
                break;
        }
        else goto load_index_failed;

        b->cached++;
    }

    if ( f->itr ) {
        hts_itr_destroy(f->itr);
        f->itr = NULL;
    }

    return b->cached;

  load_index_failed:
    error("Failed to reload index of %s. This error perhaps caused by BUGs. Please report this to shiquan@genomics.cn.", f->fname);
    
}
static struct anno_vcf_buffer *anno_vcf_buffer_init()
{
    struct anno_vcf_buffer *b = malloc(sizeof(*b));
    memset(b, 0, sizeof(*b));
    b->last_rid = -1;
    b->vcmp = vcmp_init();
    return b;
}
static void anno_vcf_buffer_destroy(struct anno_vcf_buffer *b)
{
    int i;
    for ( i = 0; i < b->max; ++i )
        bcf_destroy(b->buffer[i]);
    if ( b->vcmp )
        vcmp_destroy(b->vcmp);
    if ( b->tmpi )    free(b->tmpi);    
    if ( b->tmpi2 )   free(b->tmpi2);
    if ( b->tmpi3 )   free(b->tmpi3);
    if ( b->tmpf )    free(b->tmpf);
    if ( b->tmpf2 )   free(b->tmpf2);
    if ( b->tmpf3 )   free(b->tmpf3);
    if ( b->tmps )    free(b->tmps);
    if ( b->tmps2 )   free(b->tmps2);
    if ( b->tmpks.m ) free(b->tmpks.s);
    free(b);
}
struct anno_vcf_file *anno_vcf_file_init(bcf_hdr_t *hdr, const char *fname, char *column)
{
    assert(column); 
    kstring_t str = {0,0,0};
    kputs(column, &str);
    int n = 0;
    int *s = ksplit(&str, ',', &n);
    // if no tags specified
    if ( n == 0 ) {
        free(str.s);
        return NULL;
    }

    struct anno_vcf_file *f = malloc(sizeof(*f));
    memset(f, 0, sizeof(*f));
    f->fname = fname;
    f->fp = hts_open(fname, "r");
    if ( f->fp == NULL ) 
        error("%s : %s.", fname, strerror(errno));
    htsFormat type = *hts_get_format(f->fp);
    if ( type.format != vcf && type.format != bcf )
        error("Unsupport file type, only accept VCF/BCF. %s", fname);

    if ( f->fp->format.compression != bgzf )
        error("This file is NOT compressed by bgzip. %s", fname);

    BGZF *b = hts_get_bgzfp(f->fp);
    if ( b && bgzf_check_EOF(b) == 0 ) {
        warnings("No BGZF EOF marker, file may be truncated. %s", fname);
    }

    if ( type.format == bcf ) {
        f->bcf_idx = bcf_index_load(fname);
        if ( f->bcf_idx == NULL)
            error("Failed to load bcf index of %s.", fname);
    }
    else {
        f->tbx_idx = tbx_index_load(fname);
        if ( f->tbx_idx == NULL )
            error("Failed to load tabix index of %s.", fname);            
    }
    
    f->hdr = bcf_hdr_read(f->fp);
    
    f->cols = malloc(n*sizeof(struct anno_col));
    int i;
    kstring_t temp = {0,0,0};
    for ( i = 0; i < n; ++i ) {
        struct anno_col *col = &f->cols[i];
        char *ss = str.s + s[i];
        col->replace = REPLACE_MISSING;
        if ( *ss == '+' ) ss++;
        else if ( *ss == '-' ) {
            col->replace = REPLACE_EXISTING;
            ss++;
        }
        if ( ss == NULL || *ss == '\0' || *ss == ' ' || strlen(ss) == 0 ) {
            warnings("Bad columan format. %s", column);
            continue;
        }
        // Skip all build-in tags
        if ( strcasecmp("chrom", ss) == 0 || strcasecmp("pos", ss) == 0 || strcasecmp("ref", ss) == 0 || strcasecmp("alt", ss) == 0 || strcasecmp("qual", ss) == 0 ) {
            warnings("Skip build in tag : %s", ss);
            continue;
        }
        if ( strncasecmp("FORMAT/", ss, 7) == 0 || strncasecmp("FMT/", ss, 4) == 0 || strcasecmp("FOMAT",ss) == 0 ) {             
            warnings("DO NOT support annotate FORMAT tags.");
            continue;
        }
        if ( strcasecmp("INFO", ss) == 0 ) {
            warnings("DO NOT support all INFO tags. Please use INFO/TAG to specify target tags.");
            continue;
        }

        if ( strcasecmp("FILTER", ss) == 0 ) {
            col->func.vcf = vcf_setter_filter;
        }
        else if ( strcasecmp("ID", ss) == 0 ) {
            col->func.vcf = vcf_setter_id;
        }
        else {
            if ( strncasecmp("INFO/", ss, 5) == 0 ) 
                ss += 5;
            int id = bcf_hdr_id2int(hdr, BCF_DT_ID, ss);
            if ( !bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id) ) {
                bcf_hrec_t *rec = bcf_hdr_get_hrec(f->hdr, BCF_HL_INFO, "ID", ss, NULL);
                if ( rec == NULL ) {
                    warnings("Tag \"%s\" is not defined in header. %s", ss, fname);
                    continue;
                }
                temp.l = 0;
                bcf_hrec_format(rec, &temp);
                bcf_hdr_append(hdr, temp.s);
                bcf_hdr_sync(hdr);
                id = bcf_hdr_id2int(hdr, BCF_DT_ID, ss);
                assert (bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));
            }
            switch ( bcf_hdr_id2type(hdr, BCF_HL_INFO, id)) {
                case BCF_HT_FLAG: col->func.vcf = vcf_setter_info_flag; break;
                case BCF_HT_INT:  col->func.vcf = vcf_setter_info_int; break;
                case BCF_HT_REAL: col->func.vcf = vcf_setter_info_real; break;
                case BCF_HT_STR:  col->func.vcf = vcf_setter_info_str; break;
                default: error("Tag \"%s\" of type not recongized (%d). ", ss, bcf_hdr_id2type(hdr, BCF_HL_INFO, id)); 
            }
            col->number = bcf_hdr_id2length(hdr, BCF_HL_INFO, id);
        } // end else
        col->hdr_key = strdup(ss);
        f->n_col++;
    }
    free(s);
    if ( temp.m ) free(temp.s);
    free(str.s);
    f->buffer = anno_vcf_buffer_init();
    if ( f->n_col == 0 ) {
        anno_vcf_file_destroy(f);
        return NULL;
    }    
    return f;
}

struct anno_vcf_file *anno_vcf_file_duplicate(struct anno_vcf_file *f)
{
    struct anno_vcf_file *d = malloc(sizeof(*d));
    memset(d, 0, sizeof(*d));
    d->fname = f->fname;
    d->fp = hts_open(f->fname, "r");
    d->hdr = bcf_hdr_read(d->fp);
    if ( f->bcf_idx )
        d->bcf_idx = bcf_index_load(f->fname);
    else if ( f->tbx_idx )
        d->tbx_idx = tbx_index_load(f->fname);
    else
        error("Try to copy from a empty anno_vcf_file.");
    d->n_col = f->n_col;
    d->cols = malloc(d->n_col*sizeof(struct anno_col));
    int i;
    for ( i = 0; i < d->n_col; ++i )
        anno_col_copy(&f->cols[i], &d->cols[i]);
    d->buffer = anno_vcf_buffer_init();
    return d;
}

void anno_vcf_file_destroy(struct anno_vcf_file *f)
{
    hts_close(f->fp);
    bcf_hdr_destroy(f->hdr);
    if ( f->bcf_idx)
        hts_idx_destroy(f->bcf_idx);
    else if ( f->tbx_idx )
        tbx_destroy(f->tbx_idx);

    if ( f->itr )
        hts_itr_destroy(f->itr);

    int i;
    for ( i = 0; i < f->n_col; ++i ) free(f->cols[i].hdr_key);
    free(f->cols);
    anno_vcf_buffer_destroy(f->buffer);
    free(f);
}

int anno_vcf_core(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line)
{
    int i, j;
    
    if ( anno_vcf_update_buffer(f, hdr, line) == 0 )
        return 0;
    struct anno_vcf_buffer *b = f->buffer;
    for ( j = 0; j < b->cached; ++j ) {
        bcf1_t *d = b->buffer[j];
        if ( d->rlen != line->rlen )
            continue;
        if ( d->pos != line->pos )
            continue;
        if ( match_allele(line, d) )
            continue;
        
        for ( i = 0; i < f->n_col; ++i ) {
            struct anno_col *col = &f->cols[i];
            col->curr_name = bcf_seqname(hdr, line);
            col->curr_line = line->pos+1;
            if ( col->func.vcf(f, hdr, line, col, d) )
                warnings("Failed to annotate %s:%d with %s.", col->curr_name, col->curr_line, f->fname);
        }
    }
    
    return 0;
}


#ifdef ANNO_VCF_MAIN

#include "anno_thread_pool.h"
#include "anno_pool.h"
#include "number.h"
#include <unistd.h>

int usage()
{
    fprintf(stderr, "bcfanno_vcf [options] in.vcf\n");
    fprintf(stderr, " -data <data.bcf>     BED-like format database with header.\n");
    fprintf(stderr, " -tag  <tag,tag>      Specify tags.\n");    
    fprintf(stderr, " -t [1]               Threads.\n");
    fprintf(stderr, " -O <u|v|b|z>         Output format.\n");
    fprintf(stderr, " -o <output.vcf>      Output file.\n");
    fprintf(stderr, " -r [1000]            Record per thread per time.\n");
    return 1;
}
static const char *hts_bcf_wmode(int file_type)
{
    if ( file_type == FT_BCF ) return "wbu";    // uncompressed BCF
    if ( file_type & FT_BCF ) return "wb";      // compressed BCF
    if ( file_type & FT_GZ ) return "wz";       // compressed VCF
    return "w";                                 // uncompressed VCF
}
struct args {
    const char *input_fname;
    const char *output_fname;
    const char *data_fname;
    int n_thread;
    struct anno_vcf_file **files;
    htsFile *fp_input;
    htsFile *fp_out;
    bcf_hdr_t *hdr_out;
    int n_record;
} args = {
    .input_fname  = NULL,
    .output_fname = NULL,
    .data_fname   = NULL,
    .n_thread     = 1,
    .files        = NULL,
    .fp_input     = NULL,
    .fp_out       = NULL,
    .hdr_out      = NULL,
    .n_record     = 100,
};

void *anno_vcf(void *arg, int idx) {
    struct anno_pool *pool = (struct anno_pool*)arg;
    struct args *args = (struct args*)pool->arg;
    struct anno_vcf_file *f = args->files[idx];
    int i;
    for ( i = 0; i < pool->n_reader; ++i) {
        bcf1_t *line = pool->readers[i];
        if ( bcf_get_variant_types(line) == VCF_REF ) continue;
        bcf_unpack(line, BCF_UN_INFO);
        anno_vcf_core(f, args->hdr_out, line);
    }
    return pool;
}

int parse_args(int argc, char **argv)
{
    int i;
    if ( argc == 1 )
        return usage();

    const char *tags        = 0;
    const char *thread      = 0;
    const char *record      = 0;
    const char *output_type = 0;
    
    for ( i = 1; i < argc; ) {
        
        const char *a = argv[i++];
        const char **var = 0;
        
        if ( strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0 )
            return usage();
        else if ( strcmp(a, "-data") == 0 || strcmp(a, "-bed") == 0 )
            var = &args.data_fname;
        else if ( strcmp(a, "-tag") == 0 )
            var = &tags;
        else if ( strcmp(a, "-t") == 0 )
            var = &thread;
        else if ( strcmp(a, "-o") == 0 )
            var = &args.output_fname;
        else if ( strcmp(a, "-O") == 0 )
            var = &output_type;
        else if ( strcmp(a, "-r") == 0 )
            var = &record;

        if ( var != 0 ) {
            if ( i == argc)
                error("Missing an argument after %s", a);
            *var = argv[i++];
            continue;
        }

        if ( args.input_fname != 0 )
            error("Unknown argument, %s", a);

        args.input_fname = a;
    }

    if ( args.data_fname == 0 )
        error("Please specify bed database with -data.");

    if ( args.input_fname == 0 && (!isatty(fileno(stdin))) )
        args.input_fname = "-";

    if ( args.input_fname == 0 )
        error("No input file.");

    args.fp_input = hts_open(args.input_fname, "r");
    if ( args.fp_input == NULL )
        error("%s : %s.", args.input_fname, strerror(errno));

    htsFormat type = *hts_get_format(args.fp_input);
    if ( type.format != vcf && type.format != bcf )
        error("Unsupport input format, only accept VCF/BCF.");

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
            default:
                error("The output type \"%d\" unrecognised", out_type);
        }
    }
    args.fp_out = args.output_fname == 0 ? hts_open("-", hts_bcf_wmode(out_type)) : hts_open(args.output_fname, hts_bcf_wmode(out_type));
    if ( thread )
        args.n_thread = str2int((char*)thread);
    if ( record )
        args.n_record = str2int((char*)record);
    if ( args.n_thread < 1)
        args.n_thread = 1;
    if ( args.n_record < 1 )
        args.n_record = 100;

    args.hdr_out = bcf_hdr_read(args.fp_input);

    if ( args.hdr_out == NULL )
        error("Failed to parse header of input.");
    
    args.files = malloc(args.n_thread*sizeof(void*));
    args.files[0] = anno_vcf_file_init(args.hdr_out, args.data_fname, (char*)tags);
    
    for ( i = 1; i < args.n_thread; ++i ) 
        args.files[i] = anno_vcf_file_duplicate(args.files[0]);
    
    bcf_hdr_write(args.fp_out, args.hdr_out);
    
    return 0;
}

int bcfanno_vcf()
{
    struct thread_pool *p = thread_pool_init(args.n_thread);
    struct thread_pool_process *q = thread_pool_process_init(p, args.n_thread*2, 0);
    struct thread_pool_result  *r;
    
    for ( ;; ) {
        struct anno_pool *arg = anno_reader(args.fp_input, args.hdr_out, args.n_record);
        if ( arg->n_reader == 0 )
            break;
        arg->arg = &args;
        int block;
        do {
            block = thread_pool_dispatch2(p, q, anno_vcf, arg, 1);
            if ( ( r = thread_pool_next_result(q) ) ) {
                // generate output
                struct anno_pool *data = (struct anno_pool*)r->data;
                int i;
                for ( i = 0; i < data->n_reader; ++i ) {
                    bcf_write1(args.fp_out, args.hdr_out, data->readers[i]);
                    bcf_destroy(data->readers[i]);
                }
                free(data->readers);
                thread_pool_delete_result(r, 1);
            }
            // flush output
        } while ( block == -1);
    }
    
    thread_pool_process_flush(q);
    while ( (r = thread_pool_next_result(q)) ) {
        // generate output
        struct anno_pool *data = (struct anno_pool*)r->data;
        int i;
        for ( i = 0; i < data->n_reader; ++i ) {
            bcf_write1(args.fp_out, args.hdr_out, data->readers[i]);
            bcf_destroy(data->readers[i]);
        }
        free(data->readers);
        thread_pool_delete_result(r, 1);
    }
    thread_pool_process_destroy(q);
    thread_pool_destroy(p);
    return 0;
}

void release_memory()
{
    hts_close(args.fp_input);
    hts_close(args.fp_out);
    bcf_hdr_destroy(args.hdr_out);
    int i;
    for ( i = 0; i < args.n_thread; ++i )
        anno_vcf_file_destroy(args.files[i]);
    free(args.files);
}

int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    bcfanno_vcf();

    release_memory();

    return 0;
}

#endif
