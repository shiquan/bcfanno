#include "anno.h"
#include "anno_bed.h"

struct name_list {
    int l, m;
    char **a;
};

static void name_list_push(struct name_list *names, char *name)
{
    int i;
    for ( i = 0; i < names->l; ++i ) {
	if ( strcmp(names->a[i], name) == 0 )
	    return;
    }
    if ( names->m == names->l ) {
	names->m = names->m == 0 ? 2 : names->m + 2;
	names->a = (char**)realloc(names->a, sizeof(char*)*names->m);
    }
    names->a[names->l++] = (char*)strdup(name);
}

static char *generate_funcreg_string (struct anno_bed_file *file, struct anno_col *col)
{
    
    
}   

int beds_parse_header(bcf_hdr_t *hdr, struct anno_bed_file *file)
{
}

int beds_fill_buffer(struct anno_bed_file *file, bcf_hdr_t *hdr_out, bcf1_t *line)
{
}

int beds_options_init(struct beds_options *opts)
{
}
int beds_options_destroy(struct beds_options *opts)
{
}
int beds_databases_add(struct beds_options *opts, const char *fname, char *columns)
{
}
bcf1_t *anno_beds_core(struct beds_options *opts, bcf1_t *line)
{
}



#ifdef _BED_ANNO_MAIN
int main(int argc, char **argv)
{
    prase_args(argc, argv);

    // read line
    // annotate function region
}
#endif
