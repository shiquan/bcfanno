#include "anno.h"
#include "plugin.h"

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

struct anno_handler hand = { .files=NULL, .hdr=NULL, .hdr_out=NULL, .sql_anno_count=0, .vcf_anno_count=0, .connects=NULL };

void init_data(const char *json, const char *fname)
{
    if ( load_config(json) ) {
	error("Load JSON file failed");
    }
    // assume input is VCF/BCF file
    if ( !bcf_sr_add_reader(hand->files, fname) ) {
	error("Failed to open %s: %s\n", fname, bcf_sr_strerror(hand->files->errnum));
    }
    int i;
    for (i=0; i<anno_config_file.n_apis; ++i) {
	if (anno_config_file.apis[i].type == api_is_vcf) {
	    if (!bcf_sr_add_reader(hand->files, anno_config_file.apis[i].vfile)) {
		error("Failed to open %s: %s\n", fname, bcf_sr_strerror(hand->files->errnum));
	    }
	    parse_columns(anno_config_file.apis[i]);
	} else if (anno_config_file.apis[i].type == api_is_sql) {
	    if (!add_sql_reader(hand->connects, anno_config_file.apis[i].dynlib)) {
		error("Failed to open %s: %s\n", fname, str_errno());
	    }
	    parse_columns(anno_config_file.apis[i]);
	}
    }
}

void init_buffers(int start, int end_pos)
{
}
// check the environment parameters
void check_environ_paras()
{
}

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
	    
bcf1_t * anno_core(bcf1_t *line)
{
    int i, j;
    int x = bcf_get_variant_types(line); // get the variant type
    VARSTAT(x);
    // get the start position of variant
    int len = 0;
    for (i=1; i<line->n_allele; i++) {
	if ( len > line->d.var[i].n ) {
	    len = line->d.var[i].n;
	}
    }
    
    int end_pos = len<0 ? line->pos - len: line->pos;
    init_buffers(line->pos, end_pos);
    for (i=0; i<hand->sql_anno_count; ++i) {
	anno_sql_line(hand->connects[i], line);
    }
    for (i=0; i<hand->vcf_anno_count; ++i) {
	anno_vcf_line(i, line);
    }
    return line;
}
#undef VARSTAT

static char *input_fname = NULL;
static char *config_fname = NULL;

int paras_praser(int argc, char **argv)
{

    assert(input_fname != NULL);
    htsFile *fp = hts_open(input_fname, "r");
    htsFormat type = *hts_get_format(fp);
    hts_close(fp);
    if (type.format!=vcf && type.format!=bcf) {
	error("Unsupport input format! %s", input_fname);
    }
    return 0;
}

int main(int argc, char **argv)
{
    if ( para_praser(argc, argv) ) {
	return EXIT_FAILURE;
    }
    check_environ_paras();
    init_data(config_fname, input_fname);
    while( bcf_sr_next_line(hand->files))
    {
	if ( !bcf_sr_has_line(hand->files, 0) ) continue;
	bcf1_t *line = bcf_sr_get_line(hand->files, 0);
	if ( line->errcode ) {
	    error("Encountered error, cannot proceed. Please check the error output above.\n");
	}
	anno_core(line);
	bcf_write1(hand->out_fh, hand->hdr_out, line);
    }
    release_handler();
    return EXIT_SUCCESS;
}
