#include "utils.h"
#include "anno_bed.h"
#include "anno_col.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "number.h"
#include "stack_lite.h"

// icol start from 0
static char *tsv_get_column(struct anno_bed_tsv *t, int icol)
{
    if ( icol >= t->n_field ) {
        warnings("Out of column.");
        return NULL;
    }
    return t->string.s + t->fields[icol];
}
static char *func_region_string_generate(struct anno_bed_file *file, bcf1_t *line, struct anno_col *col)
{
    struct anno_bed_buffer *buffer = file->buffer;
    if ( buffer->cached == 0 )
        return NULL;

    struct anno_stack *s = anno_stack_init();
    int i;
    int end = line->pos +line->rlen;
    for ( i = buffer->i; i < buffer->cached; ++i ) {
        struct anno_bed_tsv *t = buffer->buffer[i];
        if ( buffer->end_pos_for_skip == 0 || buffer->end_pos_for_skip < t->end ) buffer->end_pos_for_skip = t->end;
        if ( line->pos > buffer->end_pos_for_skip ) {
            buffer->i = i;
            continue;
        }

        if ( end < t->start ) break;

        if ( line->pos < t->start || line->pos >= t->end ) continue;
        
        char *name = tsv_get_column(t, col->icol);
        if ( name == NULL )
            warnings("Failed to retrieve record, %s, %s, %s, %d.", file->fname, col->hdr_key, col->curr_name, col->curr_line);
        else
            anno_stack_push(s, name);    
    }
    kstring_t string = {0,0,0};
    for ( i = 0; i < s->l; ++i ) {
        if ( i ) kputc(',', &string);
        kputs(s->a[i], &string);
    }
    anno_stack_destroy(s);
    return string.s;
}
static struct anno_bed_tsv *anno_bed_tsv_init()
{
    struct anno_bed_tsv *t = malloc(sizeof(*t));
    t->n_field = 0;
    t->fields = NULL;
    t->string.l = t->string.m =0;
    t->string.s = 0;
    return t;
}
static void anno_bed_tsv_clean(struct anno_bed_tsv *t)
{
    if ( t->n_field )
        free(t->fields);
    t->fields = NULL;
    t->n_field = 0;
    t->start = -1;
    t->end = -1;
    t->string.l = 0;
}
static void anno_bed_tsv_destroy(struct anno_bed_tsv *t)
{
    if ( t->n_field )
        free(t->fields);
    if ( t->string.m )
        free(t->string.s);
    free(t);
}
static int string2tsv(struct anno_bed_tsv *t)
{
    int *arr = ksplit(&t->string, '\t', &t->n_field);
    if ( arr )
        t->fields = arr;
    else {
        warnings("Failed to split %s, retry ..", t->string.s);
        arr = ksplit(&t->string, '\t', &t->n_field);
        if ( arr != NULL )
            t->fields = arr;
        else {
            warnings("Failed again, this may caused by insufficient memory.");
            goto failed_convert;
        }
    }
    assert(t->n_field > 3);
    t->start = str2int(t->string.s + t->fields[1]);
    t->end   = str2int(t->string.s + t->fields[2]);
    return 0;

  failed_convert:
    return 1;
}

static int anno_bed_update_buffer(struct anno_bed_file *file, bcf_hdr_t *hdr, bcf1_t *line)
{
    assert(file->idx);
    // fill buffer
    struct anno_bed_buffer *buffer = file->buffer;
    int tid;
    tid = tbx_name2id(file->idx, bcf_seqname(hdr, line));
    // if overlap flag is set and cached this region already, skip re-fill buffer
    if ( file->overlapped == 0 &&
         buffer->last_rid  == tid &&
         buffer->last_start <= line->pos +1 &&
         buffer->last_end > line->pos )
        return -1;

    if ( tid == -1 ) {
        if ( buffer->no_such_chrom == 0 ) {
            warnings("No chromosome %s found in database %s.", bcf_seqname(hdr, line), file->fname);
            buffer->no_such_chrom = 1;
        }
        return 1;
    }
    else buffer->no_such_chrom = 0; 

    // reset buffer
    buffer->cached = 0;
    buffer->i = 0;
    buffer->end_pos_for_skip = 0;
        
    hts_itr_t *itr = tbx_itr_queryi(file->idx, tid, line->pos, line->pos + line->rlen);
    if ( itr == NULL )
        return 1;

    buffer->last_rid = tid;
    buffer->last_start = -1;
    buffer->last_end = -1;    

    for ( ;; ) {

        if ( buffer->cached == buffer->max ) {
            buffer->max += 8;
            buffer->buffer = realloc(buffer->buffer, sizeof(void*)*buffer->max);
            int i;
            for ( i = 8; i > 0; --i)
                buffer->buffer[buffer->max-i] = anno_bed_tsv_init();
        }

        struct anno_bed_tsv *t = buffer->buffer[buffer->cached];
        anno_bed_tsv_clean(t);
        
        if ( tbx_itr_next(file->fp, file->idx, itr, &t->string) < 0 )
            break;
        
        if ( string2tsv(t) )
            continue;
        
        // Skip if variant located outside of target region.
        if ( line->pos < t->start || line->pos >= t->end )
            continue;
        if ( t->end - t->start == 1 && line->pos != t->start )
            continue;

        buffer->cached++;

        if ( buffer->last_end == -1 ) {
            buffer->last_end   = t->end;
            buffer->last_start = t->start;
            continue;
        }
        if ( buffer->last_end < t->end)
            buffer->last_end = t->end;
        if ( buffer->last_start > t->start )
            buffer->last_start = t->start;            
    }

    hts_itr_destroy(itr);
    
    return buffer->cached;
}

int anno_bed_update_buffer_chunk(struct anno_bed_file *f, bcf_hdr_t *hdr, struct anno_pool *pool)
{
    assert(pool->n_reader > 0);
    // first line
    bcf1_t *line = pool->curr_line;
    struct anno_bed_buffer *b = f->buffer;
    // reset all cached records
    b->cached = 0;
    b->i = 0;
    b->end_pos_for_skip = 0;
    
    int tid;
    tid = tbx_name2id(f->idx, bcf_seqname(hdr, line));    

    if ( tid == -1 ) {
        if ( b->no_such_chrom == 0 ) {
            warnings("No chromosome %s found in database %s.", bcf_seqname(hdr, line), f->fname);
            b->no_such_chrom = 1;
        }
        return 0;
    }
    else b->no_such_chrom = 0; 
    
    // re-fill buffer
    b->cached = 0;
    b->i = 0;
    
    hts_itr_t *itr = tbx_itr_queryi(f->idx, tid, pool->curr_start, pool->curr_end+1);
    if ( itr == NULL )
        return 0;

    b->last_rid = tid;
    b->last_start = -1;
    b->last_end = -1;    

    for ( ;; ) {

        if ( b->cached == b->max ) {
            b->max += 8;
            b->buffer = realloc(b->buffer, sizeof(void*)*b->max);
            int i;
            for ( i = 8; i > 0; --i)
                b->buffer[b->max-i] = anno_bed_tsv_init();
        }

        struct anno_bed_tsv *t = b->buffer[b->cached];
        anno_bed_tsv_clean(t);
        
        if ( tbx_itr_next(f->fp, f->idx, itr, &t->string) < 0 )
            break;
        
        if ( string2tsv(t) )
            continue;

        b->cached++;

        if ( b->last_end == -1 ) {
            b->last_end   = t->end;
            b->last_start = t->start;
            continue;
        }
        if ( b->last_end < t->end) b->last_end = t->end;
        if ( b->last_start > t->start) b->last_start = t->start;
    }

    hts_itr_destroy(itr);
    
    return b->cached;

}
int anno_bed_core(struct anno_bed_file *file, bcf_hdr_t *hdr, bcf1_t *line)
{
    int i;
    // no record found
    if ( anno_bed_update_buffer(file, hdr, line) == 0)
        return 0;
    
    for ( i = 0; i < file->n_col; ++i ) {
        struct anno_col *col = &file->cols[i];
        col->curr_name = bcf_seqname(hdr, line);
        col->curr_line = line->pos+1;
        if ( col->func.bed(file, hdr, line, col) ) {
            warnings("Failed to update %s for record %s:%d.", col->hdr_key, col->curr_name, col->curr_line);
        }
    }
    return 0;
}

int anno_bed_chunk(struct anno_bed_file *f, bcf_hdr_t *hdr, struct anno_pool *pool )
{
    int i = 0, j = 0;
    struct anno_bed_buffer *b = f->buffer;
    b->cached = 0;

    if ( anno_bed_update_buffer_chunk(f, hdr, pool) == 0 )
        return 0;
        
    for ( i = pool->i_chunk; i < pool->n_chunk; ++i) {
        bcf1_t *line = pool->readers[i];

        if ( bcf_get_variant_types(line) == VCF_REF ) continue;
        bcf_unpack(line, BCF_UN_INFO);
            
        for ( j = 0; j < f->n_col; ++j ) {
            struct anno_col *col = &f->cols[j];
            col->curr_name = bcf_seqname(hdr, line);
            col->curr_line = line->pos+1;
            if ( col->func.bed(f, hdr, line, col) ) 
                warnings("Failed to update %s for record %s:%d.", col->hdr_key, col->curr_name, col->curr_line);
        }
    }
    return 0;
}
// setter functions
static int anno_bed_setter_info_string(struct anno_bed_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col) {
    struct anno_bed_buffer *b = file->buffer;
    int ret;
    if ( col->replace == REPLACE_MISSING ) {
        ret = bcf_get_info_string(hdr, line, col->hdr_key, &b->tmps, &b->mtmps);
        if ( ret > 0 && (b->tmps[0]!='.'||b->tmps[1]!= 0)) return 0;
    }
    char *string = func_region_string_generate(file, line, col);
    if ( string == NULL ) return 0;
    ret = bcf_update_info_string_fixed(hdr, line, col->hdr_key, string);
    free(string);
    return ret;
}

static int anno_bed_setter_info_int32(struct anno_bed_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col) {
    struct anno_bed_buffer *b = file->buffer;
    int ret;
    if ( col->replace == REPLACE_MISSING ) {
        ret = bcf_get_info_string(hdr, line, col->hdr_key, &b->tmps, &b->mtmps);
        if ( ret > 0 && (b->tmps[0]!='.'||b->tmps[1]!= 0)) return 0;
    }
    char *string = func_region_string_generate(file, line, col);
    if ( string == NULL ) return 0;
    int v = str2int(string);
    ret = bcf_update_info_int32_fixed(hdr, line, col->hdr_key, &v, 1);
    free(string);
    return ret;
}

static int anno_bed_setter_info_flag(struct anno_bed_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col) {
    struct anno_bed_buffer *b = file->buffer;
    int ret;
    if ( col->replace == REPLACE_MISSING ) {
        ret = bcf_get_info_string(hdr, line, col->hdr_key, &b->tmps, &b->mtmps);
        if ( ret > 0 && (b->tmps[0]!='.'||b->tmps[1]!= 0)) return 0;
    }
    char *string = func_region_string_generate(file, line, col);
    if ( string == NULL ) return 0;
    if ( string[0] == '.' && string[1] == 0 ) { free(string); return 0; }
    
    ret = bcf_update_info_flag(hdr, line, col->hdr_key, NULL, 1);
    free(string);
    return ret;
}

static int anno_bed_setter_info_float(struct anno_bed_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col) {
    struct anno_bed_buffer *b = file->buffer;
    int ret;
    if ( col->replace == REPLACE_MISSING ) {
        ret = bcf_get_info_string(hdr, line, col->hdr_key, &b->tmps, &b->mtmps);
        if ( ret > 0 && (b->tmps[0]!='.'||b->tmps[1]!= 0)) return 0;
    }
    char *string = func_region_string_generate(file, line, col);
    if ( string == NULL ) return 0;
    float v = atof(string);
    ret = bcf_update_info_float_fixed(hdr, line, col->hdr_key, &v, 1);
    free(string);
    return ret;
}

struct anno_bed_file *anno_bed_file_init(bcf_hdr_t *hdr, const char *fname, char *column)
{
    struct anno_bed_file *f = malloc(sizeof(*f));
    memset(f, 0, sizeof(*f));
    f->fname = fname;
    f->fp = hts_open(fname, "r");
    if ( f->fp == NULL )
        error("%s : %s.", fname, strerror(errno));
    f->idx = tbx_index_load(fname);
    if ( f->idx == NULL )
        error("Failed to load index of %s.", fname);

    // do not check the last regions
    f->overlapped = 1;

    // init buffer
    struct anno_bed_buffer *b = malloc(sizeof(*b));
    b->no_such_chrom = 0;
    b->last_rid = -1;
    b->last_start = -1;
    b->last_end = -1;
    b->cached = 0;
    b->max = 0;
    b->buffer = NULL;

    f->buffer = b;

    int no_columns = 0;
    // if no column specified, annotate all tags
    if ( column == NULL ) {
        no_columns = 1;
    }
    else {
        kstring_t temp = {0,0,0};
        int i, n, *s;
        kputs(column, &temp);
        s = ksplit(&temp, ',', &n);
        f->cols = malloc(n *sizeof(struct anno_col));
        for ( i = 0; i < n; ++i ) {
            char *ss = temp.s + s[i];
            struct anno_col *col = &f->cols[f->n_col];
            col->icol = -1;
            col->replace = REPLACE_MISSING;
            if ( *ss == '+' ) ss++;
            else if ( *ss == '-' ) {
                col->replace = REPLACE_EXISTING;
                ss++;
            }
            // in case empty tag
            if (ss[0] == '\0') continue;            // emit capped INFO mark
            if ( strncmp(ss, "INFO/", 5) == 0 )
                ss += 5;
            col->hdr_key = strdup(ss);
            f->n_col++;
        }
        free(temp.s);
        free(s);
    }
    int i; // public iter    
    int m = 0;     // max cache size for cols, used if no_column == 1
    kstring_t string = {0,0,0};
    for ( ;; ) {
        string.l = 0;
        if ( hts_getline(f->fp, '\n', &string) < 0 )
            break;
        // access header line in the beginning of file
        if ( string.s[0] != '#' )
            break;
        if ( strncmp(string.s, "##INFO=", 7) == 0 ) {
            char *ss = string.s + 11;
            char *se = ss;
            while ( se && *se != ',' ) se++;
            struct anno_col *col = NULL;
            if ( no_columns == 1 ) {
                if ( f->n_col == m ) {
                    m = m == 0 ? 2 : m+2;
                    f->cols = realloc(f->cols, m*sizeof(struct anno_col));
                }
                col = &f->cols[f->n_col];
                col->hdr_key = strndup(ss, se-ss+1);
                col->hdr_key[se-ss] = '\0';
                f->n_col++;
            }
            else {
                for ( i = 0; i < f->n_col; ++i ) {
                    if ( strncmp(f->cols[i].hdr_key, ss, se-ss) == 0 && strlen(f->cols[i].hdr_key) == se -ss )
                        break;
                }
                // tag is not set in the column, skip this header line 
                if ( i == f->n_col )
                    continue;
                col = &f->cols[i];
            }
            col->icol = -1;
            bcf_hdr_append(hdr, string.s);
            bcf_hdr_sync(hdr);

            int hdr_id = bcf_hdr_id2int(hdr, BCF_DT_ID, col->hdr_key);
            int ret = bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, hdr_id);            
            if ( ret == 0 ) error("Assert failure %d. %s, %s", ret, col->hdr_key, string.s);

            switch ( bcf_hdr_id2type(hdr, BCF_HL_INFO, hdr_id) ) {
                case BCF_HT_FLAG:
                    col->func.bed = anno_bed_setter_info_flag;
                    break;

                case BCF_HT_INT:
                    col->func.bed = anno_bed_setter_info_int32;
                    break;
                    
                case BCF_HT_REAL:
                    col->func.bed = anno_bed_setter_info_float;
                    break;

                case BCF_HT_STR:
                    col->func.bed = anno_bed_setter_info_string;
                    break;

                default:
                    error("Tag \"%s\" of type not recongized (%d). ", col->hdr_key, bcf_hdr_id2type(hdr, BCF_HL_INFO, hdr_id)); 
            }
        }
        
        // check the column names
        if ( strncasecmp(string.s, "#chr", 4) == 0 ) {
            int n, k;
            int *s = ksplit(&string, '\t', &n);

            if ( n < 4 ) error("Bad header format of %s", fname);

            for ( k = 3; k < n; ++k ) {
                char *ss = string.s + s[k];
                for ( i = 0; i < f->n_col; ++i ) {
                    struct anno_col *col = &f->cols[i];
                    if ( strcmp(col->hdr_key, ss) == 0 )
                        break;
                }
                // if name line specify more names than column
                if ( i == f->n_col )
                    continue;

                struct anno_col *col = &f->cols[i];
                // specify column
                col->icol = k;
            }
            free(s);
        }        
    }

    // check inited columns
    for ( i = 0; i < f->n_col; ++i ) {
        struct anno_col *col = &f->cols[i];
        if ( col->hdr_key && col->icol == -1 )
            error("Column %s not found in %s", col->hdr_key, fname);

        int hdr_id = bcf_hdr_id2int(hdr, BCF_DT_ID, col->hdr_key);
        assert(hdr_id >-1);
        col->number = bcf_hdr_id2length(hdr, BCF_HL_INFO, hdr_id);
        if ( col->number == BCF_VL_A || col->number == BCF_VL_R || col->number == BCF_VL_G )
            error("Only support fixed INFO number for tag %s. Please reset type of it.", col->hdr_key);        
    }

    free(string.s);
    
    return f;
}

struct anno_bed_file *anno_bed_file_duplicate(struct anno_bed_file *f)
{
    struct anno_bed_file *d = malloc(sizeof(*d));
    memset(d, 0, sizeof(*d));
    
    d->fname = f->fname;

    // not directly copy the address, but reopen file and index, because tabix is NOT thread-safe !    
    d->fp = hts_open(f->fname, "r");
    assert(d->fp);
    d->idx = tbx_index_load(f->fname);
    assert(d->idx);
    d->overlapped = f->overlapped;

    // init buffer
    struct anno_bed_buffer *b = malloc(sizeof(*b));
    b->no_such_chrom = 0;
    b->last_rid = -1;
    b->last_start = -1;
    b->last_end = -1;
    b->cached = 0;
    b->max = 0;
    b->buffer = NULL;
    d->buffer = b;
    
    d->n_col = f->n_col;
    d->cols = malloc(d->n_col*sizeof(struct anno_col));
    int i;
    for ( i = 0; i < d->n_col; ++i )
        anno_col_copy(&f->cols[i], &d->cols[i]);

    return d;
}

void anno_bed_file_destroy(struct anno_bed_file *f)
{
    hts_close(f->fp);
    tbx_destroy(f->idx);
    int i;
    for ( i = 0; i < f->n_col; ++i ) free(f->cols[i].hdr_key);
    free(f->cols);
    struct anno_bed_buffer *b = f->buffer;
    for ( i = 0; i < b->max; ++i ) {
        struct anno_bed_tsv *t = b->buffer[i];
        anno_bed_tsv_destroy(t);
    }
    //if ( b->tmps ) free(b->tmps);
    if ( b->buffer ) free(b->buffer);
    free(b);
    free(f);
}

#ifdef ANNO_BED_MAIN

#include "anno_thread_pool.h"
#include "anno_pool.h"
#include <unistd.h>

int usage()
{
    fprintf(stderr, "bcfanno_bed [options] in.vcf\n");
    fprintf(stderr, " -bed <data.bed.gz>   BED-like format database with header.\n");
    fprintf(stderr, " -tag <tag,tag>       Specify tags.\n");    
    fprintf(stderr, " -t [1]               Threads.\n");
    fprintf(stderr, " -O <u|v|b|z>         Output format.\n");
    fprintf(stderr, " -o <output.vcf>      Output file.\n");
    fprintf(stderr, " -r [1000]            Record per thread per time.\n");
    fprintf(stderr, " -overlap             Record in database is overlapped.\n");
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
    struct anno_bed_file **files;
    htsFile *fp_input;
    htsFile *fp_out;
    bcf_hdr_t *hdr_out;
    int n_record;
    int overlapped;
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
    .overlapped   = 0,
};

void *anno_bed(void *arg, int idx) {
    struct anno_pool *pool = (struct anno_pool*)arg;
    struct args *args = (struct args*)pool->arg;
    struct anno_bed_file *f = args->files[idx];
    int i;
    for ( i = 0; i < pool->n_reader; ++i) 
        anno_bed_core(f, args->hdr_out, pool->readers[i]);
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
        else if ( strcmp(a, "-overlap") == 0 ) {
            args.overlapped = 1;
            continue;
        }

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
    args.fp_out = args.output_fname == 0 ?
        hts_open("-", hts_bcf_wmode(out_type)) :
        hts_open(args.output_fname, hts_bcf_wmode(out_type));

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

    bcf_hdr_write(args.fp_out, args.hdr_out);
    
    args.files = malloc(args.n_thread*sizeof(void*));
    args.files[0] = anno_bed_file_init(args.hdr_out, args.data_fname, (char*)tags);
    
    for ( i = 1; i < args.n_thread; ++i ) 
        args.files[i] = anno_bed_file_duplicate(args.files[0]);

    return 0;
}

int bcfanno_bed()
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
            block = thread_pool_dispatch2(p, q, anno_bed, arg, 1);
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
        anno_bed_file_destroy(args.files[i]);
    free(args.files);
}

int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    bcfanno_bed();

    release_memory();

    return 0;
}

#endif
