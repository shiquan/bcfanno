/*  
    Copyright (C) 2016,2017  BGI Research

    Author: Shi Quan (shiquan@genomics.cn)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE. 
*/

#include "plugin.h"
#include "anno.h"
#include <errno.h>
#include <htslib/kstring.h>

struct plugin_specs specs = {
    .n = 0,
    .m = 0,
    .specs = NULL,
};

static int pl_info_update_int(struct anno_col *col, bcf1_t *line, void *data)
{
}
static int pl_info_update_string(struct anno_col *col, bcf1_t *line, void *data)
{
}
static int pl_info_update(struct anno_col *col, bcf1_t *line, void *_data)
{
    struct process_results *data = (struct process_results*)_data;
    int i;
    for (i =0; i < data->n_row; ++i) {
        // if return empty, skip
        if ( data->data[i][col->id] == NULL)
            continue;
        
        
    }
    return 0;
}

// update BCF header structure
int init_header(bcf_hdr_t *hdr, struct header_cols *cols, struct plugin_spec *spec)
{
    return 0;
}

int parse_plugin_specs(bcf_hdr_t *hdr, const char *fname, char *columns, int *n)
{
    if ( columns == NULL )
        return 1;
    
    if (specs.n == specs.m) {
        specs.m += 2;
        specs.specs = (struct plugin_spec*)realloc(specs.specs, specs.m * sizeof(struct plugin_spec));
    }
    struct plugin_spec *spec = &specs.specs[specs.n];
    void *handle = dlopen(fname, RTLD_NOW);
    if ( handle == NULL ) {
        fprintf(stderr, "Failed to load %s : %s.", fname, dlerror());
        return 1;
    }
    dlerror();
    char *ret;
    spec->funcs.init = (dl_init_func)dlsym(handle, "init");
    ret = dlerror();
    if ( ret ) {
        fprintf(stderr, "Failed to find init() : %s", ret);
        return 1;
    }
    spec->funcs.about = (dl_init_func)dlsym(handle, "about");
    ret = dlerror();
    if ( ret ) {
        fprintf(stderr, "Failed to find about() : %s", ret);
        return 1;
    }
    spec->funcs.process = (dl_init_func)dlsym(handle, "process");
    ret = dlerror();
    if ( ret ) {
        fprintf(stderr, "Failed to find process() : %s", ret);
        return 1;
    }
    spec->funcs.final = (dl_init_func)dlsym(handle, "final");
    ret = dlerror();
    if ( ret ) {
        fprintf(stderr, "Failed to find final() : %s", ret);
        return 1;
    }

    // initize
    struct header_cols *cols = (struct header_cols*)spec->funcs.init();
    if ( init_header(hdr, cols, spec) ) {
        fprintf(stderr, "Failed to update BCF header.\n");
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
        void *data = spec->funcs.process(name, start, end, NULL);
        for (j = 0; j < spec->ncols; ++j) 
            spec->cols[j].setter(&spec->cols[i], line, data);
    }
    return line;
}

int close_plugins(void)
{
    int i;
    for (i = 0; i < specs.n; ++i ) {
        struct plugin_spec *spec = &specs.specs[i];
        spec->funcs.final();
    }
    if ( specs.m )
        free(specs.specs);
    return 0;    
}

