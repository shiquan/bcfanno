#include "anno.h"
#include "plugin.h"
#include "version.h"

struct args {
    const char *fname_input;
    const char *fname_output;
    bcf_hdr_t *hdr, hdr_out;
    htsFile *fp;
    htsFile *fp_out;
    int output_type;
    struct anno_aux *anno;    
};

struct args args = {
    .fname_input = 0,
    .fname_output = 0,
    .hdr = NULL,
    .hdr_out = NULL,
    .fp = NULL,
    .fp_out = NULL,
    .output_type = 0,
    .anno = NULL,
};


int destroy_cols_vector(struct annot_cols_vector *v)
{
    if (v == NULL)
	return 0;
    int i, j;
    for ( i=0; i<v->n; i++ ) {
	struct annot_cols *cols = &v->vcols[i];
	for (j=0; j<cols->ncols; ++j) {
	    free(cols->cols[j].hdr_key);
	}
	free(cols->cols);
    }
    free(v->vcols);
    free(v);
    return 0;    
}

const char *hts_bcf_wmode(int file_type)
{
    if ( file_type == FT_BCF ) return "wbu";    // uncompressed BCF
    if ( file_type & FT_BCF ) return "wb";      // compressed BCF
    if ( file_type & FT_GZ ) return "wz";       // compressed VCF
    return "w";                                 // uncompressed VCF
}
int release_handler()
{
    bcf_sr_destroy(hand.files);
    bcf_hdr_destroy(hand.hdr_out);
    destroy_cols_vector(hand.vcf_cols);
    destroy_cols_vector(hand.sql_cols);
    hts_close(hand.out_fh);
    ignore_free(hand.sample_map);
    ignore_free(hand.tmpi);
    ignore_free(hand.tmpi2);
    ignore_free(hand.tmpi3);
    ignore_free(hand.tmpf);
    ignore_free(hand.tmpf2);
    ignore_free(hand.tmpf3);
    ignore_free(hand.tmps);
    ignore_free(hand.tmps2);
    ignore_free(hand.tmpp);
    ignore_free(hand.tmpp2);
    ignore_free(hand.tmpks.s);
    return 0;    
}
struct annot_cols_vector * acols_vector_init()
{
    struct annot_cols_vector *v = (struct annot_cols_vector*)malloc(sizeof(struct annot_cols_vector));
    v->m = v->n = 0;
    v->vcols = 0;
    return v;
}
int parse_columns(struct vcf_sql_api *api)
{
    if (api == NULL) 
	error("empty api");

    assert(api->columns);
    assert(hand.files->readers[hand.ti].header);
    assert(hand.hdr_out);

    if (api->type == anno_is_vcf) {

	int ncols = 0;
	struct annot_cols_vector *v = hand.vcf_cols;

	if (v->m == v->n) {
	    v->m = v->n == 0 ? 2: v->n << 1;
	    v->vcols = (struct annot_cols*)realloc(v->vcols, v->m *sizeof(struct annot_cols));
	}
	struct annot_cols *ac = &v->vcols[v->n++];

	ac->cols = init_columns((const char*)api->columns, hand.files->readers[hand.ti].header, hand.hdr_out, &ncols, anno_is_vcf);
	
	if (ac->cols == NULL || ncols == 0) {
	    warnings("failed to prase %s",api->columns);
	    ac->cols = NULL;
	    ac->ncols = 0;
	}
	ac->ncols = ncols;
    } else {
	/* TODO: api.type == api_is_sql */
	assert(1);
    }
    return 0;
}
int init_data(const char *json, const char *fname)
{
    if ( load_config(json) != 0 ) {
	error("Load JSON file failed");
    }
    // assume input is VCF/BCF file
    hand.files = bcf_sr_init();
    hand.files->require_index = 1;
    hand.vcf_cols = acols_vector_init();
    hand.sql_cols = acols_vector_init();
    
    if ( !bcf_sr_add_reader(hand.files, fname) ) {
	error("Failed to open %s: %s\n", fname, bcf_sr_strerror(hand.files->errnum));
    }

    hand.hdr =  hand.files->readers[0].header;
    hand.hdr_out = bcf_hdr_dup(hand.hdr);
    hand.vcmp = vcmp_init();
    hand.tmpks.l = hand.tmpks.m = 0;

    int i;
    for (i=0; i<anno_config_file.n_apis; ++i) {
	/* only accept vcf/bcf for now*/
	if (anno_config_file.apis[i].type == anno_is_vcf) {
	    
	    if ( !bcf_sr_add_reader(hand.files, anno_config_file.apis[i].vfile) ) 
		error("Failed to open %s: %s\n", fname, bcf_sr_strerror(hand.files->errnum));
	    
	    hand.ti = i;
	    parse_columns(&anno_config_file.apis[i]);
	    
	}
	/* else if (anno_config_file.apis[i].type == api_is_sql) { */
	/*     assert(1); */
	/*     if (!add_sql_reader(hand.connects, anno_config_file.apis[i].dynlib))  */
	/* 	error("Failed to open %s: %s\n", fname, str_errno()); */
	    
	/*     parse_columns(anno_config_file.apis[i]); */
	/* } */
    }
    return 0;
}

/* void init_buffers(int start, int end_pos) */
/* { */
/* } */
/* // check the environment parameters */
/* void check_environ_paras() */
/* { */
/* } */

#define VARSTAT(x) do \
    {				\
	sites_stat.all_sites++; \
	if (x == VCF_REF) { \
	    sites_stat.refs++;\
	} else if (x == VCF_SNP) {\
	    sites_stat.snps++;\
	} else if (x == VCF_MNP) {\
	    sites_stat.mnps++;\
	} else if (x == VCF_INDEL) {\
	    sites_stat.indels++;\
	} else if (x == VCF_OTHER) {\
	    sites_stat.other++;\
	} else {\
	    error("VARSTAT: Unknown flag %d", x);\
	}\
    } while(0)
	    
bcf1_t *anno_core(bcf1_t *line)
{
    int i, j;


    for (i=0; i < hand.vcf_cols->n; ++i) {
	hand.ti = i;
	if ( bcf_sr_has_line(hand.files, hand.ti) ) {
	    bcf1_t *aline = bcf_sr_get_line(hand.files, hand.ti);
	    struct annot_cols *cols = &hand.vcf_cols->vcols[i];

	    for (j=0; j<cols->ncols; ++j) {
		annot_col_t *col = &cols->cols[j];
		if ( col->setter(&hand, line, col, aline) ) 
		    warnings("failed to annotate %s at %s:%d",col->hdr_key, bcf_seqname(hand.hdr, line), line->pos+1);
	    }
	}
    }
    return line;
}
#undef VARSTAT

#include <getopt.h>

static struct option opts[] = {
    {"output", required_argument, NULL, 'o'},
    {"output-type", required_argument, NULL, 'O'},
    {"config", required_argument, NULL, 'c'},
    {"help", required_argument, NULL, 'h'},
    {NULL, 0, NULL, 0},    
};

int usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About : Annotate VCF/BCF file.\n");
    fprintf(stderr, "Version : %s, build with htslib version : %s\n", VCFANNO_VERSION, hts_version());
    fprintf(stderr, "Usage : vcfanno -c config.json in.vcf.gz\n");
    fprintf(stderr, "   -c, --config <file>            configure file, include annotations and tags, see man page for details\n");
    fprintf(stderr, "   -o, --output <file>            write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Homepage: https://github.com/shiquan/vcfanno\n");
    fprintf(stderr, "\n");
    return 1;
}
int paras_praser(int argc, char **argv)
{
    int c;
    char *config_file = NULL;
    char *input_fname = NULL;
    hand.out = "-";
    if (argc == 1)
	return usage();
    while ( (c=getopt_long(argc, argv, "O:o:c:h", opts, NULL))>=0 ) {
	switch (c) {
	    case 'o':
		hand.out = strdup(optarg);
		break;
	    case 'O':
		switch (optarg[0]) {
		    case 'b': hand.output_type = FT_BCF_GZ; break;
		    case 'u': hand.output_type = FT_BCF; break;
		    case 'z': hand.output_type = FT_VCF_GZ; break;
		    case 'v': hand.output_type = FT_VCF; break;
		    default : error("The output type \"%s\" not recognised\n", optarg);
		};
		break;
	    case 'c':
		config_file = strdup(optarg);
		break;
		
	    case 'h':
		return usage();
		
	    default:
		error("Unknown argument : %s", optarg);
	}
    }

    if (argc == optind) 
	input_fname = "-";
    else
	input_fname = argv[optind];

    if (argc - optind > 1)
	error("Only accept one VCF/BCF file. Please use `bcftools merge` to merge all the VCF/BCF files first.");

    htsFile *fp = hts_open(input_fname, "r");

    if (fp == NULL)
	error("Failed to open : %s", input_fname);
    
    htsFormat type = *hts_get_format(fp);
    
    if (type.format != vcf && type.format != bcf) {
	error("Unsupported input format! %s", input_fname);
    }

    hts_close(fp);
    hand.out_fh = hts_open(hand.out, hts_bcf_wmode(hand.output_type));
    if (hand.out_fh == NULL)
	error("cannot write to %s", hand.out);
    
    return 0;
}

void bcf_srs_regions_update(bcf_sr_regions_t *reg, const char *chr, int start, int end)
{
    if (start == -1 || end == -1) return;

    start--; end--; // 1based to 0based
    if ( !reg->seq_hash )
	reg->seq_hash = khash_str2int_init();

    int iseq;
    if ( khash_str2int_get(reg->seq_hash, chr, &iseq) < 0 ) {
	iseq = reg->nseqs++;
	
    }
}

int main(int argc, char **argv)
{
    // init argvs
    if ( paras_praser(argc, argv) != 0)
	return EXIT_FAILURE;

    
    /* hand.files = bcf_sr_init(); */
    /* hand.files->require_index = 1; */
    /* hand.vcf_cols = acols_vector_init(); */
    /* hand.sql_cols = acols_vector_init(); */
    /* hand.vcmp = vcmp_init(); */
    /* hand.tmpks.l = hand.tmpks.m = 0; */


    args.fp = hts_open(args.input_fname, "r");
    args.hdr = bcf_hdr_read(args.fp);
    hand.hdr_out = bcf_hdr_dup(args.hdr);

    int i;
    for (i = 0; i > anno_config_file.n_apis; ++i) {
	if (anno_config_file.apis[i].type == anno_is_vcf) {
	    if ( !bcf_sr_add_reader(hand.files, anno_config_file.apis[i].vfile) )
		error("Failed to open %s : %s !", anno_config_file.apis[i].vfile, bcf_sr_strerror(hand.files->errnum));

	    hand.ti = i;
	    parse_columns(&anno_config_file.apis[i]);
	    
	}	
    }

    
    if ( bcf_hdr_write(hand.out_fh, hand.hdr_out) != 0 )
	error("failed to write header.");

    bcf1_t *line = bcf_init();

    while ( bcf_read(fp, hdr, line) == 0 ) {

	if (line->errcode)
	    error("Encountered error, cannot proceed. Please check the error output above.");

	
		  
    }
    
    while ( bcf_sr_next_line(hand.files)) {
	
	if ( !bcf_sr_has_line(hand.files, 0) ) continue;
	bcf1_t *line = bcf_sr_get_line(hand.files, 0);

	if ( line->errcode ) {
	    error("Encountered error, cannot proceed. Please check the error output above.\n");
	}

	anno_core(line);
		
	bcf_write1(hand.out_fh, hand.hdr_out, line);
    }
    
    release_handler();
    return EXIT_SUCCESS;
}

