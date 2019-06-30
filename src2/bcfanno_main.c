#include "utils.h"
#include "anno_pool.h"
#include "anno_bed.h"
#include "anno_vcf.h"
#include "anno_col.h"
#include "anno_thread_pool.h"
#include "config.h"
#include "htslib/hts.h"
#include "htslib/tbx.h"
#include "htslib/vcf.h"

#include "number.h"
#include "anno_flank.h"

// for genepredext format, this format has been instead by GenomeElementAnnotation file.
//#include "genepred.h"
//#include "anno_hgvs.h"

// for GenomeElementAnnotation file
#include "anno_seqon.h"

#include "bcfanno_version.h"
#include <unistd.h>

struct anno_index {
    // point to hdr_out, DO NOT free it
    bcf_hdr_t *hdr_out;
    // bed handlers
    int n_bed;
    struct anno_bed_file **bed_files;
    // vcf handlers
    int n_vcf;
    struct anno_vcf_file **vcf_files;
    // hgvs handler. for genepredext file, will be instead by genome element annotation file
    // struct anno_hgvs_file *hgvs;
    // struct to access GenomeElementAnnotation file.
    struct anno_mc_file *mc_file;
    // flank sequence
    struct seqidx *seqidx;
};

extern int bcf_add_flankseq(struct seqidx *idx, bcf_hdr_t *hdr, bcf1_t *line);

static const char *hts_bcf_wmode(int file_type)
{
    if ( file_type == FT_BCF ) return "wbu";    // uncompressed BCF
    if ( file_type & FT_BCF ) return "wb";      // compressed BCF
    if ( file_type & FT_GZ ) return "wz";       // compressed VCF
    return "w";                                 // uncompressed VCF
}

int usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About : Annotate VCF/BCF file.\n");
    fprintf(stderr, "Version : %s, build with htslib version : %s\n", BCFANNO_VERSION, hts_version());
    fprintf(stderr, "Usage : bcfanno -c config.json in.vcf.gz\n");
    fprintf(stderr, "   -c, --config <file>            configure file, include annotations and tags, see man page for details\n");
    fprintf(stderr, "   -o, --output <file>            write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "   -q                             quiet mode\n");
    fprintf(stderr, "   -r  [number]                   records per thread. Default is %d.\n", RECORDS_PER_CHUNK);    
    fprintf(stderr, "   -t, --thread                   thread\n");
    fprintf(stderr, "   --unsort                       set if input is not sorted by cooridinate, **bad performance**\n");
    fprintf(stderr, "   --flank                        if set this flag and reference genome specified in configure, FLKSEQ tag will be generated\n");
    fprintf(stderr, "   --mito                         set the mitochodrial sequence name, default is chrM. Human mito use a different genetic code map!\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Homepage: https://github.com/shiquan/bcfanno\n");
    fprintf(stderr, "\n");
    return 1;
}

struct args {

    int test_databases_only;
    
    // input vcf path
    const char *fname_input;
    
    // output vcf path, stdout in default
    const char *fname_output;
    
    // configure path in json format
    const char *fname_json;
    
    // vcf header of input
    bcf_hdr_t *hdr;
    
    //bcf_hdr_t *hdr_out;
    
    // file handler of input vcf
    htsFile *fp_input;
    
    // file handler of output vcf
    htsFile *fp_out;
    
    // output format, default is vcf
    int output_type;
    
    // cache all arguments
    kstring_t commands;

    // give warnings instead of abortion 
    int quiet;
    
    struct bcfanno_config *config;

    int input_unsorted;
    // if this flag and reference genome is set, FLKSEQ will be annotated
    int flank_seq_is_need;
    
    // records to cache per thread
    int n_record;
    
    int n_thread;
    struct anno_index **indexs;

    uint64_t total_record;
} args = {
    .test_databases_only = 0,
    .fname_input  = NULL,
    .fname_output = NULL,
    .fname_json   = NULL,
    .hdr          = NULL,
    //.hdr_out      = NULL,
    .fp_input     = NULL,
    .fp_out       = NULL,
    .output_type  = 0,
    .commands     = {0, 0, 0},
    .quiet        = 0,
    .n_thread     = 1,
    .input_unsorted = 0,
    .flank_seq_is_need = 0,
    .n_record     = RECORDS_PER_CHUNK,
    .indexs       = NULL,
    .total_record = 0,
};

static int annotation_file_is_gea_format = 0;

struct anno_index *anno_index_init(bcf_hdr_t *hdr, struct bcfanno_config *config)
{
    extern int bcf_header_add_flankseq(bcf_hdr_t *hdr);
    extern struct seqidx *sequence_index_duplicate(struct seqidx *idx);
    extern struct seqidx* load_sequence_index(const char *file);
    
    struct anno_index *idx = malloc(sizeof(*idx));
    memset(idx, 0, sizeof(*idx));
    struct bed_config *bed_config = &config->bed;
    struct vcf_config *vcf_config = &config->vcf;
    struct refgene_config *refgene_config = &config->refgene;
    int i;
    if ( bed_config->n_bed > 0 ) {
        idx->bed_files = malloc(bed_config->n_bed *sizeof(void*));
        for ( i = 0; i < bed_config->n_bed; ++i ) 
            idx->bed_files[i] = anno_bed_file_init(hdr, bed_config->files[i].fname, bed_config->files[i].columns);
        idx->n_bed = bed_config->n_bed;
    }
    else idx->n_bed = 0;
    
    if ( vcf_config->n_vcf > 0 ) {
        idx->vcf_files = malloc(vcf_config->n_vcf*sizeof(void*));
        for ( i = 0; i < vcf_config->n_vcf; ++i )
            idx->vcf_files[i] = anno_vcf_file_init(hdr, vcf_config->files[i].fname, vcf_config->files[i].columns);
        idx->n_vcf = vcf_config->n_vcf;
    }
    else idx->n_vcf = 0;

    // idx->hgvs = NULL;
    idx->mc_file = NULL;
    
    if ( refgene_config->genepred_fname && refgene_config->refseq_fname ) {
        if ( file_is_GEA(refgene_config->genepred_fname) == 0 ) {
            idx->mc_file = anno_mc_file_init(hdr, refgene_config->columns, refgene_config->genepred_fname, refgene_config->refseq_fname, config->reference_path, refgene_config->trans_list_fname);
            annotation_file_is_gea_format = 1;
        }
        else {
            // idx->hgvs = anno_hgvs_file_init(hdr, refgene_config->columns, refgene_config->genepred_fname, refgene_config->refseq_fname, config->reference_path, refgene_config->gene_list_fname, refgene_config->trans_list_fname);
            error("BCFANNO now use Genome Element Annotation database to predict gene and variant types. \nPlease download the updated databses from github.com/shiquan/bcfanno.");
        }
    }
    
    if ( config->reference_path ) {
        idx->seqidx = load_sequence_index(config->reference_path);
        if ( idx->seqidx ) bcf_header_add_flankseq(hdr);
    }
    else idx->seqidx = NULL;
    
    idx->hdr_out = hdr;
    
    return idx;
}
struct anno_index *anno_index_duplicate(struct anno_index *idx)
{
    extern struct seqidx *sequence_index_duplicate(struct seqidx *idx);
    
    struct anno_index *d = malloc(sizeof(*d));
    memset(d, 0, sizeof(*d));
    int i;
    d->hdr_out = idx->hdr_out;
    d->n_vcf = idx->n_vcf;
    d->n_bed = idx->n_bed;
    d->vcf_files = malloc(d->n_vcf*sizeof(void*));
    d->bed_files = malloc(d->n_bed*sizeof(void*));
    for ( i = 0; i < d->n_vcf; ++i ) d->vcf_files[i] = anno_vcf_file_duplicate(idx->vcf_files[i]);
    for ( i = 0; i < d->n_bed; ++i ) d->bed_files[i] = anno_bed_file_duplicate(idx->bed_files[i]);
    // if ( idx->hgvs ) d->hgvs = anno_hgvs_file_duplicate(idx->hgvs);
    if ( idx->mc_file ) d->mc_file = anno_mc_file_duplicate(idx->mc_file);
    if ( idx->seqidx ) d->seqidx = sequence_index_duplicate(idx->seqidx);
    return d;
}
void anno_index_destroy(struct anno_index *idx, int l)
{
    extern void sequence_index_destroy(struct seqidx *idx);
    int i;
    for ( i = 0; i < idx->n_vcf; ++i ) anno_vcf_file_destroy(idx->vcf_files[i]);
    for ( i = 0; i < idx->n_bed; ++i ) anno_bed_file_destroy(idx->bed_files[i]);
    if ( idx->vcf_files ) free(idx->vcf_files);
    if ( idx->bed_files ) free(idx->bed_files);
    // if ( idx->hgvs ) anno_hgvs_file_destroy(idx->hgvs);
    if ( idx->mc_file) anno_mc_file_destroy(idx->mc_file, l);
    if ( idx->seqidx ) sequence_index_destroy(idx->seqidx);
    free(idx);
}

// No detail message output.
static int quiet_mode = 0;

int parse_args(int argc, char **argv)
{
    int i;
    for (i = 0; i < argc; ++i ) {
	if ( i ) kputc(' ', &args.commands);
	kputs(argv[i], &args.commands);
    }    
    const char *output_fname_type = 0;
    const char *thread = 0;
    const char *record = 0;
    const char *mito = 0;
    for (i = 1; i < argc; ) {
	const char *a = argv[i++];
	if ( strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0)
	    return usage();
        
	// quiet mode
	if ( strcmp(a, "-q") == 0 || strcmp(a, "--quiet") == 0 ) {
	    quiet_mode = 1;
	    continue;
	}
        
        if ( strcmp(a, "--test_only") == 0 ) {
	    args.test_databases_only = 1;
	    continue;
	}
        if ( strcmp(a, "--unsort") == 0 ) {
            args.input_unsorted = 1;
            continue;
        }
        if ( strcmp(a, "--flank") == 0 ) {
            args.flank_seq_is_need = 1;
            continue;
        }
            
        const char **var = 0;
	if ( strcmp(a, "-c") == 0 || strcmp(a, "--config") == 0 ) 
	    var = &args.fname_json;
	else if ( strcmp(a, "-o") == 0 || strcmp(a, "--output") == 0)
	    var = &args.fname_output;
	else if ( strcmp(a, "-O") == 0 || strcmp(a, "--output-type") == 0 )
	    var = &output_fname_type;
        else if ( strcmp(a, "-t") == 0 || strcmp(a, "-thread") == 0 )
            var = &thread;
        else if ( strcmp(a, "-r") == 0 || strcmp(a, "-record") == 0 )
            var = &record;
        else if ( strcmp(a, "--mito") == 0 )
            var = &mito;
        
	if ( var != 0 ) {
	    if (i == argc) error("Missing an argument after %s", a);
	    *var = argv[i++];
	    continue;
	}
        
        if ( a[0] == '-' && a[1] ) error("Unknown parameter. %s", a);

        if ( args.fname_input == 0 ) {
	    args.fname_input = a;
	    continue;
	}
	error("Unknown argument : %s, use -h see help information.", a);
    }
    
    if ( quiet_mode == 0 ) {
        LOG_print("Version: %s + htslib-%s", BCFANNO_VERSION, hts_version());
        LOG_print("Homepage: https://github.com/shiquan/bcfanno");
	LOG_print("Args: %s", args.commands.s);
    }
    
    if ( args.fname_json == 0 ) {
	fprintf(stderr, "[error] No configure file is specified. Use -h for help message.\n");
	//fprintf(stderr, "[notice] %s.\n", DONOT_POST_ERR_STRING);
	return 1;
    }

    args.config = bcfanno_config_init();
    if ( bcfanno_load_config(args.config, args.fname_json) != 0 ) 
	error("Failed to load configure file. %s : %s", args.fname_json, strerror(errno));    
    
    if ( quiet_mode == 0 ) {
	//LOG_print("Load configure file success.");
	bcfanno_config_debug(args.config);
    }

    // if input file is not set, use stdin
    if ( args.fname_input == 0 && (!isatty(fileno(stdin))) )
        args.fname_input = "-";

    // if no detect stdin, go error
    if ( args.fname_input == 0)
        error("No input file! bcfanno only accept one BCF/VCF input file. Use -h for more informations.");

    // read input file
    args.fp_input = hts_open(args.fname_input, "r");
    if ( args.fp_input == NULL )
        error("Failed to open %s.", args.fname_input);

    // check input type is VCF/BCF or not
    htsFormat type = *hts_get_format(args.fp_input);
    if ( type.format  != vcf && type.format != bcf )
        error("Unsupported input format, only accept BCF/VCF format. %s", args.fname_input);

    if ( thread ) {
        args.n_thread = str2int((char*)thread);
        if ( args.n_thread < 1 ) args.n_thread = 1;
    }
    if ( record ) {
        args.n_record = str2int((char*)record);
        if ( args.n_record < 0 ) args.n_record = 1000;
    }
        
    // init output type
    int out_type = FT_VCF;
    if ( output_fname_type != 0 ) {
	switch (output_fname_type[0]) {
	    case 'b':
		out_type = FT_BCF_GZ; break;
	    case 'u':
		out_type = FT_BCF; break;
	    case 'z':
		out_type = FT_VCF_GZ; break;
	    case 'v':
		out_type = FT_VCF; break;
	    default :
		error("The output type \"%d\" not recognised\n", out_type);
	};
    }
    // init output file handler
    args.fp_out = args.fname_output == 0 ? hts_open("-", hts_bcf_wmode(out_type)) : hts_open(args.fname_output, hts_bcf_wmode(out_type));

//    if ( annotation_file_is_gea_format == 0 ) { // assume it is genepredext format
        // set genepredExt format
        // set_format_genepredext();
    //  }
    
    // read bcf header from input bcf/vcf
    args.hdr = bcf_hdr_read(args.fp_input);
    if ( args.hdr == NULL)
	error("Failed to parse header of input.");


    // set Mito environment
    if ( mito == NULL ) 
        setenv("BCFANNO_MITOCHR", "chrM", 1);
    else 
        setenv("BCFANNO_MITOCHR", mito, 1);

    if ( quiet_mode == 0 ) {
        const char *mito_par =  getenv("BCFANNO_MITOCHR");
        LOG_print("Set environment parameter BCFANNO_MITOCHR to %s", mito_par);
    }
        
    /******************
         INIT indexs   
     ******************/

    // allocate memeory for threads
    args.indexs = malloc(args.n_thread*sizeof(struct anno_index));    
    args.indexs[0] = anno_index_init(args.hdr, args.config);
    for ( i = 1; i < args.n_thread; ++i )
        args.indexs[i] = anno_index_duplicate(args.indexs[0]);
    
    kstring_t str = {0,0,0};
    ksprintf(&str, "##bcfannoVersion=%s+htslib-%s\n", BCFANNO_VERSION, hts_version());
    bcf_hdr_append(args.hdr, str.s);
    str.l = 0;
    ksprintf(&str, "##bcfannoCommand=%s\n", args.commands.s);
    bcf_hdr_append(args.hdr, str.s);
    
    // write header to output    
    bcf_hdr_write(args.fp_out, args.hdr);
    free(args.commands.s);
    free(str.s);
    return 0;
}

void memory_release()
{
    hts_close(args.fp_input);
    hts_close(args.fp_out);
    bcfanno_config_destroy(args.config);
    bcf_hdr_destroy(args.hdr);
    int i;
    for ( i = 0; i < args.n_thread; ++i ) anno_index_destroy(args.indexs[i], i);
    free(args.indexs);
}

void *anno_core(void *arg, int idx)
{

    assert(idx >= 0);

    struct anno_index *index = args.indexs[idx];
    struct anno_pool  *pool  = (struct anno_pool*) arg;
    
    int i, j ;
    // IMPROVE HERE: read line by line may not require sorted input but highly CPU consume, read a chunk of records
    // based on the start and end of pool will highly improve the performance
    if ( args.input_unsorted == 1 ) {
        for ( i = 0; i < pool->n_reader; ++i ) {
            bcf1_t *line = pool->readers[i];
            debug_print("%d",line->pos);
            if ( bcf_get_variant_types(line) == VCF_REF )
                continue;

            // if ( index->hgvs ) 
            // anno_hgvs_core(index->hgvs, index->hdr_out, line);

            // if ( index->mc_file ) do not support unsorted input

            for ( j = 0; j < index->n_vcf; ++j )
                anno_vcf_core(index->vcf_files[j], index->hdr_out, line);
        
            for ( j = 0; j < index->n_bed; ++j )
                anno_bed_core(index->bed_files[j], index->hdr_out, line);

            if ( args.flank_seq_is_need == 1 && index->seqidx )
                bcf_add_flankseq(index->seqidx, index->hdr_out, line);
        }
    }
    // retrieve attributes in chunk
    else {
        for ( ;; ) {
            if ( pool->n_chunk == pool->n_reader ) break;
            update_chunk_region(pool);
            
            //if ( index->hgvs )
            // anno_hgvs_chunk(index->hgvs, index->hdr_out, pool);
            if ( index->mc_file )
                anno_mc_chunk(index->mc_file, index->hdr_out, pool);
            
            for ( i = 0; i < index->n_vcf; ++i )
                anno_vcf_chunk(index->vcf_files[i], index->hdr_out, pool);
            for ( i = 0; i < index->n_bed; ++i )                
                anno_bed_chunk(index->bed_files[i], index->hdr_out, pool);
        }
        if ( args.flank_seq_is_need == 1 && index->seqidx ) {
            for ( i = 0; i < pool->n_reader; ++i) 
                bcf_add_flankseq(index->seqidx, index->hdr_out, pool->readers[i]);
        }
    }
    
    return pool;
}

int annotate_light()
{
    struct anno_index *idx = args.indexs[0];

    if ( args.input_unsorted == 1 ) {
        bcf1_t *line = bcf_init();
        int j;
        while ( bcf_read(args.fp_input, args.hdr, line) == 0 ) {
            args.total_record ++;
            if ( line->rid == -1 ) goto output_line;
            if ( bcf_get_variant_types(line) == VCF_REF) goto output_line;
            debug_print("%d",line->pos);
            //if ( idx->hgvs )
            //  anno_hgvs_core(idx->hgvs, idx->hdr_out, line);
            
            for ( j = 0; j < idx->n_vcf; ++j )
                anno_vcf_core(idx->vcf_files[j], idx->hdr_out, line);
            
            for ( j = 0; j < idx->n_bed; ++j )
                anno_bed_core(idx->bed_files[j], idx->hdr_out, line);

            if ( args.flank_seq_is_need == 1 && idx->seqidx ) bcf_add_flankseq(idx->seqidx, idx->hdr_out, line);
          output_line:
            bcf_write1(args.fp_out, args.hdr, line);
        }
    }
    else {
        for ( ;; ) {
            struct anno_pool *pool = anno_reader(args.fp_input, args.hdr, args.n_record);
            args.total_record += (uint64_t)pool->n_reader;
            if ( pool == 0 || pool->n_reader == 0) break;
            int i;
            for ( ;; ) {
                if ( pool->n_chunk == pool->n_reader ) break;
                update_chunk_region(pool);

                //  if ( idx->hgvs )
                //  anno_hgvs_chunk(idx->hgvs, idx->hdr_out, pool);
                if ( idx->mc_file )
                    anno_mc_chunk(idx->mc_file, idx->hdr_out, pool);
                
                for ( i = 0; i < idx->n_vcf; ++i )
                    anno_vcf_chunk(idx->vcf_files[i], idx->hdr_out, pool);
                for ( i = 0; i < idx->n_bed; ++i )                
                    anno_bed_chunk(idx->bed_files[i], idx->hdr_out, pool);
            }
            if ( args.flank_seq_is_need == 1 && idx->seqidx ) {
                for ( i = 0; i < pool->n_reader; ++i) 
                    bcf_add_flankseq(idx->seqidx, idx->hdr_out, pool->readers[i]);                
            }
            for ( i = 0; i < pool->n_reader; ++i) {
                bcf_write1(args.fp_out, args.hdr, pool->readers[i]);
                bcf_destroy(pool->readers[i]);
            }
            free(pool->readers);                
        }
    }
    
    return 0;
}

int annotate()
{
    if ( args.test_databases_only == 1) return 0;

    // lightweight mode
    if ( args.n_thread == 1 ) return annotate_light();

    // keep 1 thread to maintain main stream    
    args.n_thread = args.n_thread-1;
    
    // multi thread mode
    struct thread_pool *p = thread_pool_init(args.n_thread);
    struct thread_pool_process *q = thread_pool_process_init(p, args.n_thread*2, 0);
    struct thread_pool_result  *r;

    for ( ;; ) {
        struct anno_pool *arg = anno_reader(args.fp_input, args.hdr, args.n_record);
        args.total_record += (uint64_t)arg->n_reader;
        if ( arg->n_reader == 0 )
            break;
        
        int block;
        do {
            block = thread_pool_dispatch2(p, q, anno_core, arg, 1);
            if ( (r = thread_pool_next_result(q))) {
                struct anno_pool *d = (struct anno_pool*)r->data;
                int i;
                for ( i = 0; i < d->n_reader; ++i) {
                    bcf_write1(args.fp_out, args.hdr, d->readers[i]);
                    bcf_destroy(d->readers[i]);
                }
                free(d->readers);                
                thread_pool_delete_result(r, 1);
            }
        } while (block == -1);
    }

    thread_pool_process_flush(q);
    while (( r = thread_pool_next_result(q) )) {
        struct anno_pool *d = (struct anno_pool*)r->data;
        int i;
        for ( i = 0; i < d->n_reader; ++i) {
            bcf_write1(args.fp_out, args.hdr, d->readers[i]);
            bcf_destroy(d->readers[i]);
        }
        free(d->readers);
        thread_pool_delete_result(r, 1);
    }
    thread_pool_process_destroy(q);
    thread_pool_destroy(p);

    return 0;
}

#include <time.h>

int main(int argc, char **argv)
{
    clock_t t = clock();
    
    if ( parse_args(argc, argv) )
        return 1;

    if ( annotate() )
        return 1;

    memory_release();

    if ( quiet_mode == 0 ) {
        t = clock() -t;
        double time_taken = ((double)t)/CLOCKS_PER_SEC;
        LOG_print("Annotate %lld records in %.2f seconds.", args.total_record, time_taken);
    }
    return 0;
}
