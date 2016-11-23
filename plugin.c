#include "plugin.h"
#include "anno.h"
#include <errno.h>
#include <htslib/kstring.h>

struct plugin_specs specs = {
    .n = 0,
    .m = 0,
    .specs = NULL,
};

int parse_plugin_specs(const char *fname, char *columns, int *n)
{
    if ( columns == NULL )
        return 1;
    
    if (specs.n == specs.m) {
        specs.m += 2;
        specs.specs = (struct plugin_spec*)realloc(specs.specs, specs.m * sizeof(struct plugin_spec));
    }
    struct plugin_spec *spec = &specs.specs[specs.n];
    spec->handle = dlopen(fname, RTLD_NOW);
    if ( spec->handle == NULL ) {
        fprintf(stderr, "Failed to load %s : %s.", fname, dlerror());
        return 1;
    }    
    return 0;    
}

bcf1_t *plugins_process(bcf_hdr_t *hdr, bcf1_t *line)
{
    int i, j;
    // for each plugin
    char *name = bcf_seqname(hdr, line);
    int start = line->pos+1;
    int end = line->pos + line->rlen;
    for (i = 0; i < specs.n; ++i ) {
        // for each function
        struct plugin_spec *spec = &specs.specs[i];

        for (j = 0; j < spec->ncols; ++j) 
            spec->cols[j].setter(line);
        
    }
    return line;
}
