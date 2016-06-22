#include "anno.h"
#include "plugin.h"

/* stat variants positions */
struct {
    uint64_t all_sites;
    uint64_t refs;
    uint64_t snps;
    uint64_t mnps;
    uint64_t ti;
    uint64_t tv;
    uint64_t indels;
    uint64_t other;
} sites_stat = { 0, 0, 0, 0, 0, 0, 0, 0 };

//extern struct configs anno_config_file;

struct anno_handler hand = {
    .files = NULL,
    .hdr = NULL,
    .hdr_out = NULL,
    .ti = 0,
    //.sql_anno_count=0,
    .vcf_cols = NULL,
    .sql_cols = NULL,
    //.connects=NULL
    /* .sample_map = NULL, */
    /* .nsample_map = 0, */
    /* .sample_is_file = 0, */
    /* .mtmpi = 0, */
    /* .mtmpi2 = 0, */
    /* .mtmpi3 = 0, */
    /* .mtmps = 0, */
    /* .mtmps2 = 0, */
    /* .mtmps3 = 0, */
    /* .mtmpf = 0, */
    /* .mtmpf2 = 0, */
    /* .mtmpf3 = 0, */
    /* .tmpi = NULL, */
    /* .tmpi2 = NULL, */
    /* .tmpi3 = NULL, */
    /* .tmpf = NULL, */
    /* .tmpf2 = NULL, */
    /* .tmpf3 = NULL, */
    /* .tmps = NULL, */
    /* .tmps2 = NULL, */
    /* .tmpp = NULL, */
    /* .tmpp2 = NULL, */
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
	free(v->vcols);
    }
    free(v);
    return 0;    
}

int release_handler()
{
    bcf_sr_destroy(hand.files);
    bcf_hdr_destroy(hand.hdr);
    bcf_hdr_destroy(hand.hdr_out);
    destroy_cols_vector(hand.vcf_cols);
    destroy_cols_vector(hand.sql_cols);
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
	int length = 0;

	//debug_print("%s", bcf_hdr_fmt_text(hand.files->readers[hand.ti].header, 1, &length));
	
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
    if ( load_config(json) ) {
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
	    
	    hand.ti = i+1;
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

    //int x = bcf_get_variant_types(line); // get the variant type
    //VARSTAT(x);
    // get the start position of variant
    /* int len = 0; */
    /* for (i=1; i<line->n_allele; i++) { */
    /* 	if ( len > line->d.var[i].n ) { */
    /* 	    len = line->d.var[i].n; */
    /* 	} */
    /* } */
    
    /* int end_pos = len<0 ? line->pos - len: line->pos; */
    /* init_buffers(line->pos, end_pos); */
    /* for (i=0; i<hand.sql_anno_count; ++i) { */
    /* 	anno_sql_line(hand.connects[i], line); */
    /* } */

    for (i=0; i < hand.vcf_cols->n; ++i) {
	hand.ti = i+1;
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

static char *input_fname = NULL;
static char *config_fname = NULL;

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
    fprintf(stderr, "Usage: vcfanno -c config.json in.vcf.gz\n");
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
    
    return init_data(config_file, input_fname);
}


int main(int argc, char **argv)
{
    if ( paras_praser(argc, argv) ) {
	return EXIT_FAILURE;
    }
    //check_environ_paras();

    while ( bcf_sr_next_line(hand.files)) {
	
	if ( !bcf_sr_has_line(hand.files, 0) ) continue;
	bcf1_t *line = bcf_sr_get_line(hand.files, 0);

	if ( line->errcode ) {
	    error("Encountered error, cannot proceed. Please check the error output above.\n");
	}

	anno_core(line);
	bcf_write1(hand.out_fn, hand.hdr_out, line);
    }
    
    release_handler();
    return EXIT_SUCCESS;
}

