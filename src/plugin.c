/*  
    Copyright (C) 2016,2017  BGI Research

    Author: Shi Quan <shiquan@genomics.cn>

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

#include "utils.h"
#include "anno.h"
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include <ctype.h>
#include <dlfcn.h>
#include "plugin.h"

// BCFANNO self-defined APIs
// about() return help message
// init() must be called before process each record, this function will return header
// information, bcf header structure will be updated with this information.
// process() function return key value pairs for each record
// final() must be called at the end to release memory and safe return.

typedef void* (*mod_init_function)(const char *fname);
typedef int   (*mod_about_function)(void);
typedef void* (*mod_process_function)(char *name, int start, int end, int n_allele, char *alleles);
typedef int   (*mod_final_function)(void);

struct module_functions {
    mod_init_function init;
    mod_about_function about;
    mod_process_function process;
    mod_final_function final;
};

struct module_spec {
    const char *fname;
    struct module_functions mod_funcs;
};

struct args {    
    bcf_hdr_t *hdr_out;
    int m, n;
    struct module_spec *specs;    
} args = {
    .hdr_out = NULL,
    .n = 0,
    .m = 0,
    .specs = NULL,
};

static int module_init_hdr(struct module_spec *spec, bcf_hdr_t *hdr)
{
    struct pl_header *hdr_comp = (struct pl_header*)spec->mod_funcs.init();
    if ( hdr_comp == NULL )
        return 1;
    
    int i;
    for ( i = 0; i < hdr_comp->n_records; ++i) {
        char *str = hdr_comp->records[i];
        bcf_hdr_append(hdr, str);
        bcf_hdr_sync(hdr);
    }
    return i > 0 ? 0 : 1;
}

static int module_init_library(struct module_spec *spec, const char *dynapi_fname)
{
    
    spec->fname = dynapi_fname;
    
    void *handle = dlopen(dynapi_fname, RTLD_NOW);
    if ( handle == NULL )
        error("%s : %s.", dynapi_fname, strerror(errno));

    dlerror();
    char *ret;
    spec->mod_funcs.init = (mod_init_function)dlsym(handle, "init");
    ret = dlerror();
    if ( ret ) {
        error_print("Failed to find init() function. %s", ret);
        return 1;
    }

    spec->mod_funcs.about = (mod_about_function)dlsym(handle, "about");
    ret = dlerror();
    if ( ret )
        spec->mod_funcs.about = NULL;

    spec->mod_funcs.process = (mod_process_function)dlsym(handle, "process");
    ret = dlerror();
    if ( ret ) {
        error_print("Failed to find process() function. %s", ret);
        return 1;
    }

    spec->mod_funcs.final = (mod_final_function)dlsym(handle, "final");
    ret = dlerror();
    if ( ret ) {
        error_print("Failed to find final() function. %s", ret);
        return 1;
    }
    
    return 0;
}

int module_init(struct bcfanno_config *config, bcf_hdr_t *hdr)
{
    struct modules_config *mod_conf = &config->modules;
    int i;
    int nfiles = 0;
    args.m = mod_conf->n_modules;
    args.specs = (struct module_spec*)calloc(mod_conf->n_modules, sizeof(struct module_spec));
    
    for ( i = 0; i < mod_conf->n_modules; ++i ) {

        // For self-defined API, may disrupt for connection or premisson problem .. 
        // Only skip to access failed APIs instead of abort directly.
        if ( module_init_library(&args.specs[nfiles], mod_conf->files[i].fname) ) {
            warnings("Failed to init API : %s", mod_conf->files[i].fname);
            continue;
        }

        if ( module_init_hdr(&args.spec[nfiles], hdr) ) {
            warnings("Failed to init header : %s", mod_conf->files[i].fname);
            continue;
        }
        
        nfiles++;
    }

    args.n = nfiles;
    args.hdr_out = hdr;
    
    return nfields == 0 ? 1 : 0;
}

static int mod_setter_info_real()
{
}

static int mod_setter_info_int()
{
}

static int mod_setter_info_str()
{
}

static int mod_setter_info_flag()
{
}

static int mod_setter_id()
{
}

static void generate_vcf_line(struct module_data *data, bcf1_t *line)
{
    
}

int module_process(bcf_hdr_t *hdr, bcf1_t *line)
{
    int i, j;
    // for each plugin
    char *name = bcf_seqname(hdr, line);
    int start = line->pos+1;
    int end = line->pos + line->rlen;
    
    for ( i = 0; i < args.n; ++i ) {
        struct module_spec *spec = &args.specs[i];
        struct module_data *data = (struct module_data*)spec->mod_funcs.process(name, start, end, line->n_allele, line->d.allele);
        
        generate_vcf_line(data, line);
    }

    return 0;
}

int module_close(void)
{
    int i;
    for ( i = 0; i < args.n; ++i ) {
        struct module_spec *spec = &args.specs[i];
        spec->mod_funcs.final();
    }

    if ( args.m )
        free(args.specs);
    
    return 0;    
}

