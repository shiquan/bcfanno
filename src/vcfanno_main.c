/*  vcfanno - a core program to annotate vcf/bcf by using current stat databases
 *
 */

#include "utils.h"
#include "anno.h"
#include "config.h"
#include "htslib/hts.h"
#include "htslib/tbx.h"
#include "htslib/vcf.h"
#include "anno_bed.h"
#include "anno_vcf.h"
//#include "anno_hgvs.h"
#include "version.h"
#include "hgvs.h"
#include "hgvs_vcf.h"
#include "genepred.h"

// const char *vcfanno_version = "version 0.0.5";
// cache ANNOCORE_BUFFER_LINES lines into buffers; for each buffer pool all lines come from one chromosome,
// if no enough lines, just put as much as possible
#define ANNOCORE_BUFFER_LINES 1000

// time stat
#include <sys/time.h>

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
    fprintf(stderr, "Version : %s, build with htslib version : %s\n", VCFANNO_VERSION, hts_version());
    fprintf(stderr, "Usage : vcfanno -c config.json in.vcf.gz\n");
    fprintf(stderr, "   -c, --config <file>            configure file, include annotations and tags, see man page for details\n");
    fprintf(stderr, "   -o, --output <file>            write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "   -q                             quiet mode\n");
    fprintf(stderr, "   -t                             warning instead of abortion for inconsistant position\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Homepage: https://github.com/shiquan/vcfanno\n");
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
    // vcf header of output, hdr_out should also kept by beds_options, vcfs_options, and refgene_options
    bcf_hdr_t *hdr_out;
    // file handler of input vcf
    htsFile *fp_input; 
    // file handler of output vcf
    htsFile *fp_out;
    // output directory, default is "-" for stdout
    const char *out_fname;
    // output format, default is vcf
    int output_type;
    // cache all arguments
    kstring_t commands;
    // options for annotate bed format databases, usually with four columns, and suggest to put the header
    // of tags in the comment regions
    struct beds_options bed_opts;
    struct vcfs_options vcf_opts;
    // struct refgene_options hgvs_opts;

    // warning instead of abortion flag
    int tol;
};

struct args args = {
    .test_databases_only = 0,
    .fname_input = 0,
    .fname_output = 0,
    .hdr = NULL,
    .hdr_out = NULL,
    .fp_input = NULL,
    .fp_out = NULL,
    .output_type = 0,
    .commands = KSTRING_INIT,
    .tol = 0,
};

void args_destroy(){
    // debug_print("destroy header");
    bcf_hdr_destroy(args.hdr);
    // debug_print("destroy header out");
    bcf_hdr_destroy(args.hdr_out);
    // debug_print("close input");
    hts_close(args.fp_input);
    // debug_print("close output");
    bcf_close(args.fp_out);
    // debug_print("close commands");
    if (args.commands.m)
	free(args.commands.s);
    // debug_print("close vcfs");
    vcfs_options_destroy(&args.vcf_opts);
    // debug_print("close beds");
    beds_options_destroy(&args.bed_opts);
    // debug_print("close hgvs");
    // refgene_options_destroy(&args.hgvs_opts);
    close_hgvs_anno();
}

static int quiet_mode = 0;

int test_databases_framework()
{
    
    // test success
    return 0;
}
int parse_args(int argc, char **argv)
{
    int i;
    for (i = 0; i < argc; ++i ) {
	if ( i ) kputc(' ', &args.commands);
	kputs(argv[i], &args.commands);
    }
    const char *output_fname_type = 0;
    for (i = 0; i < argc; ) {
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
        if ( strcmp(a, "-t") == 0 ) {
            args.tol = 1;
            continue;
        }
	const char **var = 0;
	if ( strcmp(a, "-c") == 0 || strcmp(a, "--config") == 0 ) 
	    var = &args.fname_json;
	else if ( strcmp(a, "-o") == 0 || strcmp(a, "--output") == 0)
	    var = &args.fname_output;
	else if ( strcmp(a, "-O") == 0 || strcmp(a, "--output-type") == 0 )
	    var = &output_fname_type;

	if ( var != 0 ) {
	    if (i == argc)
		error("Missing an argument after %s", a);
	    *var = argv[i++];
	    continue;
	}
	if ( args.fname_input == 0 ) {
	    args.fname_input = a;
	    continue;
	}
	error("Unknown argument : %s, use -h see help information.", a);
    }
    
    if ( quiet_mode == 0 ) {
        LOG_print("Version: %s + htslib-%s", VCFANNO_VERSION, hts_version());
        LOG_print("Homepage: https://github.com/shiquan/vcfanno");
	LOG_print("Args: %s", args.commands.s);
    }
    
    if ( args.fname_json == 0 ) {
	fprintf(stderr, "[error] No configure file is specified. Use -h for help message.\n");
	fprintf(stderr, "[notice] %s.\n", BE_SMART_STRING);
	return 1;
    }

    struct vcfanno_config *con = vcfanno_config_init();
    if ( vcfanno_load_config(con, args.fname_json) != 0 ) {
	error("Failed to load configure file. %s : %s", args.fname_json, strerror(errno));
    }
    
    if ( quiet_mode == 0 ) {
	LOG_print("Load configure file success.");
	vcfanno_config_debug(con);
    }
    
    if ( args.test_databases_only == 1)
	return test_databases_framework();

    // if input file is not set, use stdin
    if ( args.fname_input == 0 && (!isatty(fileno(stdin))) )
	args.fname_input = "-";
    // if no detect stdin, go error
    if ( args.fname_input == 0)
	error("No input file! vcfanno only accept one BCF/VCF input file. Use -h for more informations.");
    // read input file
    args.fp_input = hts_open(args.fname_input, "r");
    if ( args.fp_input == NULL )
	error("Failed to open %s.", args.fname_input);    
    // check input type is VCF/BCF or not
    htsFormat type = *hts_get_format(args.fp_input);
    if ( type.format  != vcf && type.format != bcf )
	error("Unsupported input format, only accept BCF/VCF format. %s", args.fname_input);

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

    args.bed_opts.beds_is_inited = 0;
    args.vcf_opts.vcfs_is_inited = 0;
    // args.hgvs_opts.refgene_is_inited = 0;
    // set genepred format
    set_format_genepred();
    beds_options_init( &args.bed_opts );
    vcfs_options_init( &args.vcf_opts );
    // refgene_options_init( &args.hgvs_opts );

    // read bcf header from input bcf/vcf
    args.hdr = bcf_hdr_read(args.fp_input);
    if ( args.hdr == NULL)
	error("Failed to parse header of input.");
    // duplicate input header to generate output header
    args.hdr_out = bcf_hdr_dup(args.hdr);
    // alias hdr_out
    args.vcf_opts.hdr_out = args.hdr_out;
    args.bed_opts.hdr_out = args.hdr_out;
    // args.hgvs_opts.hdr_out = args.hdr_out;
    
    if ( con->refgene.refgene_is_set == 1) {

        if ( init_hgvs_anno(con->refgene.genepred_fname, con->refgene.refseq_fname, args.hdr_out) )
            return 1;

        if ( con->refgene.trans_list_fname != NULL ) {
            if ( set_transcripts_list( con->refgene.trans_list_fname ) ) {
                warnings("Failed to load transcripts list : %s", con->refgene.trans_list_fname);
            }
        }

        if ( con->refgene.gene_list_fname != NULL ) {
            if ( set_transcripts_list( con->refgene.gene_list_fname ) ) {
                warnings("Failed to load genes list : %s", con->refgene.gene_list_fname);
            }
        }
    }    

    for ( i = 0; i < con->vcfs.n_vcfs; ++i ) {
	vcfs_database_add(&args.vcf_opts, con->vcfs.files[i].fname, con->vcfs.files[i].columns);
    }

    for ( i = 0; i < con->beds.n_beds; ++i ) {
	beds_database_add(&args.bed_opts, con->beds.files[i].fname, con->beds.files[i].columns);
    }

    kstring_t str = { 0, 0, 0};
    ksprintf(&str, "##vcfannoVersion=%s+htslib-%s\n", VCFANNO_VERSION, hts_version());
    bcf_hdr_append(args.hdr_out, str.s);
    str.l = 0;
    ksprintf(&str, "##vcfannoCommand=%s\n", args.commands.s);
    bcf_hdr_append(args.hdr_out, str.s);
    free(str.s);
    // write header to output    
    bcf_hdr_write(args.fp_out, args.hdr_out);
    
    vcfanno_config_destroy(con);
    return 0;
}

struct varstat {
    uint64_t all_vars;
    uint64_t snps;
    uint64_t indels;
    uint64_t ti;
    uint64_t tv;
    uint64_t het_x;
    uint64_t var_x;
} varstat = {
    .all_vars = 0,
    .snps = 0,
    .indels = 0,
    .ti = 0,
    .tv = 0,
    .het_x = 0,
    .var_x = 0,
};

int anno_core(bcf1_t *line)
{
    
#ifdef DEBUG_MODE
    debug_print("%s:%d",bcf_seqname(args.hdr_out, line), line->pos+1);
#endif
    
    // do nothing for reference positions
    if ( bcf_get_variant_types(line) == VCF_REF )
	return 0;

    // Annotate hgvs name
    // anno_refgene_core(&args.hgvs_opts, line);
    if (args.tol == 0 && setter_hgvs_vcf(args.hdr_out, line) == 1)
        return 1;
    
    // todo: stat type module, ti,tv etc
    
    // annotate vcf files
    if (args.tol == 0 && anno_vcfs_core(&args.vcf_opts, line) == 1)
        return 1;
 
    // annotate bed format datasets
    if (args.tol == 0 && anno_beds_core(&args.bed_opts, line) == 1)
        return 1;

    // annotate transcript related bed format datasets
    // anno_trans_core();

    // filter set module
    // summary
    
    return 0;
}
void export_reports()
{
}

int main(int argc, char **argv)
{
    // parse arguments first, if failure or just do test will return 1, else return 0
    if ( parse_args(--argc, ++argv) == 1 )
	return 1;

    // todo: multi threads support
    // read input vcf file by line
    bcf1_t *line = bcf_init();
    while ( bcf_read(args.fp_input, args.hdr, line) == 0) {
	// skip uninited line
	if (line->rid == -1)
	    continue;
	// annotate vcf line function
        int ret = anno_core(line);
        if ( ret ) {
            fprintf(stderr, "Failed to update bcf line, %s : %d\n", bcf_seqname(args.hdr, line), line->pos+1);
            // if return 1 go abort, else give a warning
            if ( ret == 1 )
                return 1;
        }
	bcf_write1(args.fp_out, args.hdr_out, line);
    }
    if ( quiet_mode == 0 )
        LOG_print("Annotate finished. Close.");
    bcf_destroy(line);
    export_reports();
    args_destroy();
    return 0;
}
