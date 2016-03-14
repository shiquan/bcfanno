#include "anno.h"
#include "plugin.h"

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
bcf1_t * anno_core(bcf1_t *line)
{
    int i, j;
    bcf_get_variant_types(line); // get the variant type
    // get the start position of variant
    int len = 0;
    for (i=1; i<line->n_allele; i++) {
	if ( len > line->d.var[i].n ) {
	    len = line->d.var[i].n;
	}
    }
    int end_pos = len<0 ? line->pos - len: line->pos;
    init_buffers(line->pos, end_pos);
    
}

int main(int argc, char **argv)
{
    
    return EXIT_SUCCESS;
}
