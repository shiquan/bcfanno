#include "anno.h"
#include "plugin.h"


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
    int end_pos = len<0 ? line->pos - len; line->pos;
    
	
}

int main(int argc, char **argv)
{
    
    return EXIT_SUCCESS;
}
