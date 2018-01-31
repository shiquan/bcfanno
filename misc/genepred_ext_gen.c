// Check the genepred data and transcripts sequences is consistant.
#include "utils.h"
#include "genepred.h"
#include "htslib/faidx.h"
#include "htslib/tbx.h"
#include "htslib/kstring.h"
#include "ksw.h"
#include "number.h"
#include "htslib/kseq.h"
#include "htslib/khash.h"
#include "anno_thread_pool.h"
#include <pthread.h>
#include <zlib.h>

KSTREAM_INIT(gzFile, gzread, 16384)

KHASH_MAP_INIT_STR(name, char*)

typedef kh_name_t hash_t;

// skip the version number in the transcript name
static int skip_version_flag = 0;

int usage()
{
    fprintf(stderr,
            "realign_trans - realgin transcripts to reference genome and generate gap information.\n"
            "options:\n"
            " -data     Gene prediction format database.\n"
            " -rna      Transcripts sequence database included in FASTA format, indexed by samtools faidx.\n"
            " -ref      Genome reference database in FASTA format, indexed by samtools faidx.\n"
            " -format   Format of gene prediction database, genepred is default. [genepred,refflat,refgene,genepredext]\n"
            " -skip-ver Skip version number of the transcript name.\n"
            " -p        Set threads [1].\n"
            " -chrs     Rename chromosome (contig) names in database. Renamed chroms consistent with reference genome.\n"
	    " -gapo     Penalty score for gap open. [5]\n"
	    " -gape     Penalty score for gap extension. [2]\n"
	    " -sa       Penalty score for match. [2]\n"
	    " -sb       Penalty score for mismatch. [2]\n"	   
            "\nwebsite: https://github.com/shiquan/bcfanno\n");
    return 1;
}
struct index {
    faidx_t *rna_fai;
    faidx_t *ref_fai;
    struct args *args;
};

struct args {
    const char *genepred_fname;
    const char *refrna_fname;
    const char *reference_fname;
    const char *chrs_fname;
    void *hash;
    const char *format;
    int gapo;
    int gape;
    int sa;
    int sb;
    int8_t mat[25];
    int threads;
    int stable_flag;
    int chunk_size;
    kstream_t *ks;

    // thread unsafe indexs, preloaded index for each thread
    struct index *indexs;
} args = {
    .genepred_fname = NULL,
    .refrna_fname = NULL,
    .reference_fname = NULL,
    .chrs_fname = NULL,
    .hash = NULL,
    .format = "genepred",
    .gapo = 5,
    .gape = 2,
    .sa = 2,
    .sb = 2,
    .threads = 1,
    .stable_flag = 0,
    .chunk_size = 10,
    .ks = NULL,
    .indexs = NULL,
};

void *load_chrs(const char *chrs_fname) {
    
    int n, i;
    char **names = hts_readlist(chrs_fname, 1, &n);
    if ( n == 0 )
	return NULL;

    hash_t *hash = kh_init(name);
    
    kstring_t string = {0, 0, 0};
    int nfields = 0;
    
    for ( i = 0; i < n; ++i ) {
	// skip comment lines
	if (names[i] == NULL || names[i][0] == '#')
	    continue;
	// init string
	string.l = 0;
	kputs(names[i], &string);
	// split string by tab
	int *splits = ksplit(&string, '\t', &nfields);
	if (nfields != 2) {
	    warnings("Error format; only accept two columns per line. %s", names[i]);
	    continue;
	}
	char *ss = string.s + splits[0];
	char *se = string.s + splits[1];
	if ( strlen(ss) > 4 && strncmp(ss, "CTG/", 4) == 0 )
	    ss += 4;
        char *key = strdup(ss);
        int ret;
        khiter_t k = kh_get(name, hash, key);
        if ( k == kh_end(hash) ) {
            k = kh_put(name, hash, key, &ret);
            char **val = &kh_value(hash, k);
            *val = strdup(se);
        } 
        else { // for duplicated record
            free(key);
        }
    }
    return (void*)hash;
}

int parse_args(int ac, char **av)
{
    if ( ac == 1 )
        return usage();

    int i;
    const char *gapo = NULL;
    const char *gape = NULL;
    const char *sa = NULL;
    const char *sb = NULL;
    const char *threads = NULL;
    for ( i = 1; i < ac; ) {
        const char *a = av[i++];
        const char **var = 0;
        if ( strcmp(a, "-h") == 0 )
            return usage();

        if ( strcmp(a, "-data") == 0 && args.genepred_fname == NULL )
            var = &args.genepred_fname;
        else if ( strcmp(a, "-rna") == 0 && args.refrna_fname == NULL )
            var = &args.refrna_fname;
        else if ( strcmp(a, "-ref") == 0 && args.reference_fname == NULL )
            var = &args.reference_fname;
        else if ( strcmp(a, "-format") == 0 )
            var = &args.format;
        else if ( strcmp(a, "-gape") == 0 )
            var = &gape;
        else if ( strcmp(a, "-gapo") == 0 )
            var = &gapo;
        else if ( strcmp(a, "-sa") == 0 )
            var = &sa;
        else if ( strcmp(a, "-sb") == 0 )
            var = &sb;
        else if ( strcmp(a, "-p") == 0 )
            var = &threads;
        else if ( strcmp(a, "-chrs") == 0 )
            var = &args.chrs_fname;
        
        if ( var != 0 ) {
            if ( i == ac )
                error("Missing an argument after %s.", a);
            *var = av[i++];
            continue;
        }

        if ( strcmp(a, "-skip-ver") == 0 ) {
            skip_version_flag = 1;
            continue;
        }
        else if ( strcmp(a, "-stable") == 0 ) {
            args.stable_flag = 1;
            continue;
        }
        
        error("Unknown argument: %s.", a);
        return 1;          
    }

    if ( args.genepred_fname == NULL )
        error("Specify genepred database with -data.");

    if ( args.refrna_fname == NULL )
        error("Specify RNA sequence database with -rna.");

    if ( args.reference_fname == NULL )
        error("Specify reference genome database with -ref.");

    if ( args.format ) {
        if ( strcmp(args.format, "genepred") == 0 )
            set_format_genepred();
        else if ( strcmp(args.format, "refgene") == 0 )
            set_format_refgene();
        else if ( strcmp(args.format, "refflat") == 0 )
            set_format_refflat();
        else if ( strcmp(args.format, "genepredext") == 0 )
            set_format_genepredext();
        else
            error("Unknown format %s.", args.format);
    }

    if ( sa != NULL ) {
        int temp = str2int((char*)sa);
        args.sa = temp > 0 ?  temp : args.sa;
    }

    if ( sb != NULL ) {
        int temp = str2int((char*)sb);
        args.sb = temp > 0 ?  temp : args.sb;
    }

    if ( gape != NULL ) {
        int temp = str2int((char*)gape);
        args.gape = temp > 0 ?  temp : args.gape;
    }

    if ( gapo != NULL ) {
        int temp = str2int((char*)gapo);
        args.gapo = temp > 0 ?  temp : args.gapo;
    }

    if ( args.chrs_fname ) {
        args.hash = load_chrs(args.chrs_fname);
    }
    
    if ( threads != NULL ) {
        args.threads = str2int((char*)threads);
        if ( args.threads < 1 )
            args.threads = 1;
    }
    
    // Since BGZF is NOT thread safe, here only accept thread == 1.

    int j, k;
    for ( i = k = 0; i < 4; ++i ) {
        for ( j = 0; j < 4; ++j ) args.mat[k++] = i == j ? args.sa : -args.sb;
        args.mat[k++] = 0; 
    }
    for ( j = 0; j < 5; ++j ) args.mat[k++] = 0;

    // init index structure for each thread
    args.indexs = malloc(args.threads * sizeof(struct index));
    for ( i = 0; i < args.threads; ++i ) {
        struct index *idx = &args.indexs[i];
        idx->rna_fai = fai_load(args.refrna_fname);
        idx->ref_fai = fai_load(args.reference_fname);
        idx->args = &args;
    }
    
    return 0;
}

struct data {
    struct args *args;
    int n_buffers;
    kstring_t *buffers;
};

kstring_t *read_buffer(kstream_t *ks, int chunk_size, int *_n)
{
    kstring_t *buffers = NULL;
    kstring_t str = {0, 0, 0};
    int ret;
    int n, m;
    m = n = 0;
    int size = 0;
    int dret;
    for ( ;; ) {
        ret = ks_getuntil2(ks, KS_SEP_LINE, &str, &dret, 0);
        if ( ret < 0 )
            break;
        if ( n == m ) {
            m = m ? m<<1 : 16;
            buffers = (kstring_t*)realloc(buffers, m*sizeof(kstring_t));
        }
        kstring_t *s = &buffers[n];
        s->l = str.l;
        s->m = str.m;
        s->s = strdup(str.s);
        s->s[s->l] = '\0';
        if ( s->s == NULL ) {
            error("failed to allocated memory.");
        }
        size += buffers[n++].l;
        if ( size > chunk_size )
            break;
    }
    *_n = n;
    if ( str.m )
        free(str.s);
    
    return buffers;
}

char *process(struct args *args, struct index *index, kstring_t *str)
{
    // extern function from faidx_def.c
    // retrieve transcript version number from transcript reference database
    // the format of transcript title in FASTA should be format as >TRANSCRIPT[space]VERSION 
    extern int trans_retrieve_version(void *_fai, const char *trans);
    
    static const int sntab[256] = {
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    };

    struct genepred_line *line = genepred_line_create();

    parse_line(str, line);

    if ( skip_version_flag ) {
        char *s = line->name1;
        for ( ; *s != '.' && *s; ++s); 
        if (*s == '.')
            *s = '\0';
    }

    int rlen = 0;
    int ver = trans_retrieve_version((void*)index->rna_fai, line->name1);
    char *trans = faidx_fetch_seq((const faidx_t*)index->rna_fai, line->name1, 0, INT_MAX, &rlen);

    if ( trans == 0 || rlen <= 0 ) {
        warnings("%s not found in the Transcripts FASTA.", line->name1);
        return NULL;
    }

    // rename chrs
    if ( args->hash != NULL ) {
        //int ret;
        khiter_t k = kh_get(name, (hash_t*)args->hash, line->chrom);
        if ( k != kh_end((hash_t*)args->hash) ) {
            char *t = kh_val((hash_t*)args->hash, k);
            free(line->chrom);
            line->chrom = strdup(t);
        }
    }

    int i;
    kstring_t tmp = { 0, 0, 0};
    kstring_t string = { 0, 0, 0};
    for ( i = 0; i < line->exon_count; ++i ) {
        int l = 0;
        char *seq = faidx_fetch_seq((const faidx_t*)index->ref_fai, line->chrom, line->exons[BLOCK_START][i]-1, line->exons[BLOCK_END][i]-1, &l);
        if ( seq == NULL || l <= 0 ) {
            tmp.l = 0;
            break;
        }
            
        if ( l != line->exons[BLOCK_END][i] - line->exons[BLOCK_START][i]+1)
            error("Inconsistant block length.");
        if ( seq == NULL ) {
            warnings("Failed to fetch sequence from reference. %s:%d-%d.", line->chrom, line->exons[BLOCK_START][i], line->exons[BLOCK_END][i]);
            break;
        }
        kputs(seq, &tmp);
        free(seq);  
    }

    if ( tmp.l == 0 || tmp.s == 0) {
        return NULL;
    }
    else if (tmp.l == rlen) {
        genepred2line(line, &string);
        ksprintf(&string, "\t%d\t%dM", ver, rlen);
        return string.s;
    }
    
    int qlen = tmp.l;
    //reference
    uint8_t *qseq = malloc(qlen *sizeof(*qseq));
    //transcript
    uint8_t *rseq = malloc(rlen *sizeof(*rseq));
        
    for ( i = 0; i < tmp.l; ++i )
        qseq[i] = sntab[(int)tmp.s[i]];
    for ( i = 0; i < rlen; ++i )
        rseq[i] = (uint8_t)sntab[(int)trans[i]];
    if ( line->strand == '-' ) {
        for ( i = 0; i < qlen/2; ++i ) {
            int t = qseq[i] == 4 ? 4 : 3 - qseq[i]; qseq[i] = 3 - qseq[qlen-i-1]; qseq[qlen-i-1] = t;
            if ( qlen & 1 ) qseq[qlen/2+1] = 3 - qseq[qlen/2+1];
        }
    }

    uint32_t *cigar = NULL;
    int n_cigar = 0;
    int score = 0;
    score = ksw_global(rlen, rseq, qlen, qseq, 5, args->mat, args->gapo, args->gape, 100 > abs(rlen - qlen) ? 100 : abs(rlen - qlen), &n_cigar, &cigar);
    int ret = 0;

    genepred2line(line, &string);
    kputc('\t', &string); kputw(ver, &string); kputc('\t', &string);
    
    if ( n_cigar ) {
        for ( i = 0; i < n_cigar; ++i ) {
            int c = cigar[i]&0xf;
            ksprintf(&string, "%d%c", cigar[i]>>4, "MIDSH"[c]);
        }
    }
    else {
        kputs("*", &string);
        warnings("%s not properly checked. score : %d", line->name1, score);
        ret = 1;
    }

    free(qseq);
    free(rseq);
    free(cigar);
    free(trans);
    free(tmp.s);
    return string.s;
}

static void *worker(void *_data, int idx)
{
    int i;
    struct data *data = (struct data*)_data;
    assert(idx >= 0);
    struct index *index = &data->args->indexs[idx];
    kstring_t string = {0,0,0};
    for ( i = 0; i < data->n_buffers; ++i ) {
        char *str =  process(data->args, index, &data->buffers[i]);
        if ( str ) {
            kputs(str, &string);
            kputc('\n', &string);
            free(str);
        }
        free(data->buffers[i].s);
    }

    return (void*)string.s;
}

int realign_trans1()
{
    gzFile fp = gzopen(args.genepred_fname, "r");
    if ( fp == NULL )
        error("%s : %s.", args.genepred_fname, strerror(errno));
    
    args.ks = ks_init(fp);
    if ( args.threads == 0 ) {
        int ret;
        int dret;
        kstring_t str = {0,0,0};
        for ( ;; ) {
            ret = ks_getuntil2(args.ks, KS_SEP_LINE, &str, &dret, 0);
            if ( ret < 0 )
                break;
            char *seq =  process(&args, &args.indexs[0], &str);
            puts(seq);
            free(seq);
        }
        if (str.m) free(str.s);
    }
    else {
        struct thread_pool *p = thread_pool_init(args.threads);
        struct thread_pool_process *q = thread_pool_process_init(p, args.threads*2, 0);
        struct thread_pool_result *r;

        for ( ;; ) {
            struct data *data = malloc(sizeof(*data));
            data->buffers = read_buffer(args.ks, args.chunk_size, &data->n_buffers);
            if ( data->n_buffers == 0 )
                break;
            
            data->args = &args;
            
            int block;
            do {
                block = thread_pool_dispatch2(p, q, worker, data, 1);
                if ( (r = thread_pool_next_result(q)) ) {
                    char *seq = (char*)r->data;
                    if ( seq ) {
                        printf("%s", seq);
                        fflush(stdout);
                    }
                    thread_pool_delete_result(r, 1);
                }
                if ( block == -1 )
                    fflush(stdout);            
            } while ( block == -1 );
        }
        thread_pool_process_flush(q);
        
        while (( r = thread_pool_next_result(q))) {
            char *seq = (char*)r->data;
            if ( seq ) {
                printf("%s", seq);
                fflush(stdout);
            }
            thread_pool_delete_result(r, 1);
        }
        thread_pool_process_destroy(q);
        thread_pool_destroy(p);
    }
    
    ks_destroy(args.ks);
    if ( args.hash ) {
        khiter_t k;
        for ( k = kh_begin((hash_t*)args.hash); k != kh_end((hash_t*)args.hash); ++k ) {
            if ( kh_exist((hash_t*)args.hash, k) ) {
                char *val = kh_value((hash_t*)args.hash,k);
                free(val);
                kh_del(name, (hash_t*)args.hash, k);
            }
        }
        kh_destroy(name, (hash_t*)args.hash);
    }
    gzclose(fp);
    int i;
    for ( i = 0; i < args.threads; ++i) {
        struct index *idx = &args.indexs[i];
        fai_destroy(idx->rna_fai);
        fai_destroy(idx->ref_fai);        
    }
    
    return 0;    
}
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    if ( realign_trans1() )
        return 1;

    return 0;
}
