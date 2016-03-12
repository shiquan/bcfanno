#include "anno.h"
#include "plugin.h"

struct anno_handler hand = { .files=NULL, .hdr=NULL, .hdr_out=NULL, .sql_anno_count=0, .connects=NULL };

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
    for (i=0; i<anno_config_file->n_apis; ++i) {
	
    }
}
void anno_core(bcf1_t *line)
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
    
	
}

int main(int argc, char **argv)
{
    
    return EXIT_SUCCESS;
}
