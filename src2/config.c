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

// prase configure file in JSON format
// todo: stable improvement

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>

#include "kson.h"
#include "json_config.h"
#include "config.h"
#include "utils.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"

#ifndef KSTRING_INIT
#define KSTRING_INIT { 0, 0, 0}
#endif

// return 0 for locate file
//        1 for absolute directory
int check_if_fullpath(char *file)
{
    char *path = strdup(file);
    char *dir = dirname(path);
    int ret = 1;
    if ( dir == NULL || (dir[0] == '.' && dir[1] == '\0')) ret = 0;
    free(path);
    return ret;
}
char *get_dirpath(const char *file)
{
    char *path = strdup(file);
    char *dir = dirname(path);
    char *p_dir = dir == NULL ? NULL : strdup(dir);
    free(path);
    return p_dir;
}
char *generate_fullpath_file(char *file, char *dir)
{
    if ( dir == NULL ) return strdup(file);
    kstring_t str = {0,0,0};
    kputs(dir, &str);
    if ( str.l > 0 && str.s[str.l-1] != '/') kputc('/', &str);
    kputs(file, &str);
    return str.s;
}
struct bcfanno_config temp_config_init = {
    .author = 0,
    .config_id = 0,
    .reference_version = 0,
    .vcf = { 0, 0 },
    .bed = { 0, 0 },
    .refgene = { 0, 0, 0, 0, 0},
    .module = { 0, 0 },
};

struct bcfanno_config *bcfanno_config_init()
{
    struct bcfanno_config *config = (struct bcfanno_config*)malloc(sizeof(struct bcfanno_config));
    memset(config, 0, sizeof(struct bcfanno_config));
    return config;
}

void bcfanno_config_destroy(struct bcfanno_config *config)
{
    if ( config->author )
	free (config->author );
    if ( config->config_id )
	free (config->config_id );
    if ( config->reference_version )
	free (config->reference_version);
    if ( config->reference_path )
        free (config->reference_path);

    int i;
    for ( i = 0; i < config->vcf.n_vcf; ++i ) {
	free(config->vcf.files[i].fname);
	free(config->vcf.files[i].columns);
    }
    if ( i )
	free(config->vcf.files);

    for ( i = 0; i < config->bed.n_bed; ++i ) {
	free(config->bed.files[i].fname);
	if (config->bed.files[i].columns)
	    free(config->bed.files[i].columns);
    }
    if ( i )
	free(config->bed.files);

    for ( i = 0; i < config->module.n_module; ++i ) 
        free(config->module.files[i].fname);
    
    if ( config->refgene.genepred_fname )
        free(config->refgene.genepred_fname);
    if ( config->refgene.refseq_fname )
        free(config->refgene.refseq_fname);
    if ( config->refgene.trans_list_fname )
        free(config->refgene.trans_list_fname);
    if ( config->refgene.gene_list_fname )
        free(config->refgene.gene_list_fname);

    free(config);    
}

static int load_config_core(struct bcfanno_config *config, kson_t *json, char *dir)
{
    const kson_node_t *root = json->root;
    if (root == NULL)
	error("Format error. Root node is empty, check configure file.");


#define BRANCH_INIT(_node) ( (_node)->v.str == NULL ? NULL : strdup((_node)->v.str) )
#define BRANCH(_f) do {\
        if ( dir && check_if_fullpath(_f) == 0 ) {     \
            char *new = generate_fullpath_file(_f, dir);\
            free(_f);\
            _f = new;\
        }\
    } while(0)
    
    int i;
    for ( i = 0; i < root->n; ++i ) {
	const kson_node_t *node = kson_by_index(root, i);
	if ( node == NULL )
	    continue;
	if ( node->key == NULL) {
	    warnings("Format error. Node key is empty. skip..");
	    continue;
	}
	// the summary informations should in the top level
	if ( strcmp(node->key, "Author") == 0 || strcmp(node->key, "author") == 0) {
	    config->author = BRANCH_INIT(node);
	}
        else if ( strcmp(node->key, "ID") == 0 || strcmp(node->key, "id") == 0) {
	    config->config_id = BRANCH_INIT(node);
	}
        else if ( strcmp(node->key, "ref") == 0 || strcmp(node->key, "reference") == 0) {
	    config->reference_path = BRANCH_INIT(node);
            BRANCH(config->reference_path);
	}
        else if ( strcmp(node->key, "HGVS") == 0 || strcmp(node->key, "hgvs") == 0) {
	    if ( node->type != KSON_TYPE_BRACE)
		error("Format error. Configure for HGVS format should looks like :\n"
		      "\"hgvs\":{\n \"gene_data\":\"genepred.txt.gz\",\n}"
		    );
	    struct refgene_config *refgene_config = &config->refgene;
	    int j;
	    for (j = 0; j < node->n; ++j) {
		const kson_node_t *node1 = kson_by_index(node, j);
		if (node1 == NULL)
		    error("Format error. HGVS node is empty, check configure file.");

		if ( strcmp(node1->key, "gene_data") == 0 || strcmp(node1->key, "refgene") == 0) {
		    refgene_config->genepred_fname = BRANCH_INIT(node1);
                    BRANCH(refgene_config->genepred_fname);
                    /* if ( check_if_fullpath(refgene_config->genepred_fname) ){ */
                    /*     char *new = generate_fullpath_file(refgene_config->genepred_fname, dir); */
                    /*     free(refgene_config->genepred_fname); */
                    /*     refgene_config->genepred_fname = new; */
                    /* } */
                }
		else if ( strcmp(node1->key, "refseq") == 0 || strcmp(node1->key, "transcript") == 0 || strcmp(node1->key, "trans") == 0) {
		    refgene_config->refseq_fname = BRANCH_INIT(node1);
                    BRANCH(refgene_config->refseq_fname);
                }
		else if ( strcmp(node1->key, "transcripts_list") == 0 || strcmp(node1->key, "trans_list") == 0 ) {
		    refgene_config->trans_list_fname = BRANCH_INIT(node1);
                    BRANCH(refgene_config->trans_list_fname);
                }
		else if ( strcmp(node1->key, "genes_list") == 0 ) {
		    refgene_config->gene_list_fname = BRANCH_INIT(node1);
                    BRANCH(refgene_config->gene_list_fname);
                }
                else if ( strcmp(node1->key, "column") == 0 || strcmp(node1->key, "columns") == 0 ) 
                    refgene_config->columns = BRANCH_INIT(node1);
		else
		    warnings("Unknown key : %s. skip ..", node1->key);		
	    }
	    /* if ( refgene_config->columns == NULL || refgene_config->columns[0] == '\0' ) */
	    /*     error("No columns specified in HGVS configure."); */
	    if ( refgene_config->genepred_fname == NULL || refgene_config->genepred_fname[0] == '\0' )
		error("No genepred databases specified in HGVS configure.");
            if ( refgene_config->refseq_fname == NULL || refgene_config->refseq_fname[0] == '\0' )
		error("No refseq.fa specified in HGVS configure.");
            
	    refgene_config->refgene_is_set = 1;
	}
        else if ( strcasecmp(node->key, "vcf") == 0 || strcasecmp(node->key, "vcfs") == 0) {
	    if ( node->type != KSON_TYPE_BRACKET )
		error("Format error, configure for vcf databases should looks like :\n"
		      "\"vcf\":[\n{\n\"file\":\"file.vcf.gz\",\"columns\":\"TAGS,TAGS\",\n},\n{},\n]"
		    );
	    if ( node->n == 0) {
		warnings("Empty vcf configure. Skip");
		continue;
	    }
	    struct vcf_config *vcf_config = &config->vcf;
	    vcf_config->n_vcf = (int)node->n;
	    vcf_config->files = (struct file_config*)malloc(vcf_config->n_vcf*sizeof(struct file_config));
	    int n_files = 0;
	    int j;

	    for ( j = 0; j < (long)node->n; ++j ) {
		const kson_node_t *node1 = node->v.child[j];
		if ( node1 == NULL )
		    continue;
		int k;
		struct file_config *file_config = &vcf_config->files[n_files];
		file_config->fname = NULL;
		file_config->columns = NULL;
		for ( k = 0; k < node1->n; ++k ) {
		    const kson_node_t *node2 = kson_by_index(node1, k);
		    if ( node2 == NULL || node2->key == NULL)
			continue;

		    if ( strcmp(node2->key, "file") == 0 ) {
			file_config->fname = BRANCH_INIT(node2);
                        BRANCH(file_config->fname);
                    }
		    else if ( strcmp(node2->key, "columns") == 0 )
			file_config->columns = BRANCH_INIT(node2);
		    else
			warnings("Unknown key : %s. skip ..", node2->key);
		}
		// if only set vcf file, go abort
		if ( file_config->columns == NULL || file_config->columns[0] == '\0')
		    error("No columns specified for vcf. %s", file_config->fname);
		// if only set columns, skip it
		if ( file_config->fname == NULL || file_config->fname[0] == '\0' ) {
		    free(file_config->columns);
		    file_config->columns = NULL;
		    continue;
		}		
		n_files++;		
	    }
	    if (n_files == 0 && vcf_config->n_vcf)
		free(vcf_config->files);
	    vcf_config->n_vcf = n_files;
	}
        else if ( strcasecmp(node->key, "bed") == 0 || strcasecmp(node->key, "beds") == 0 ) {
	    if ( node->type != KSON_TYPE_BRACKET )
		error("Format error, configure for vcf databases should looks like :\n"
		      "\"bed\":[\n{\n\"file\":\"file.bed.gz\",\"header\":\"TAGS,TAGS\",\n},\n{},\n]"
		    );
	    if ( node->n == 0) {
		warnings("Empty bed configure. Skip ..");
		continue;
	    }
	    struct bed_config *bed_config = &config->bed;
	    bed_config->n_bed = (int)node->n;
	    bed_config->files = (struct file_config*)malloc(bed_config->n_bed*sizeof(struct file_config));
	    int n_files = 0;
	    int j;

	    for ( j = 0; j < (long)node->n; ++j ) {
		const kson_node_t *node1 = node->v.child[j];
		if ( node1 == NULL )
		    continue;
		int k;
		struct file_config *file_config = &bed_config->files[n_files];
		file_config->fname = NULL;
		file_config->columns = NULL;
		for ( k = 0; k < node1->n; ++k ) {
		    const kson_node_t *node2 = kson_by_index(node1, k);
		    if (node2 == NULL || node2->key == NULL)
			continue;
		    if ( strcmp(node2->key, "file") == 0 ) {
			file_config->fname = BRANCH_INIT(node2);
                        BRANCH(file_config->fname);
                    }
		    else if ( strcmp(node2->key, "columns") == 0 )
			file_config->columns = BRANCH_INIT(node2);
		    else
			warnings("Unknown key : %s. skip ..", node1->key);
		}
		// if only set vcf file, check the header of bed file, description information keep in the header in default
			    
		// if only set columns, skip it
		if ( file_config->columns && file_config->fname == NULL ) {
		    free(file_config->columns);
		    file_config->columns = NULL;
		    continue;
		}
		n_files++;		
	    }
	    if ( n_files == 0 && bed_config->n_bed )
		free(bed_config->files);
	    bed_config->n_bed = n_files;	    
	}
        else if ( strcmp(node->key, "module") == 0 || strcmp(node->key, "module") == 0 || strcmp(node->key, "plugins") == 0 || strcmp(node->key, "plugin") == 0 ) {
	    if ( node->type != KSON_TYPE_BRACKET )
		error("Format error, configure for vcf databases should looks like :\n"
		      "\"module\":[\n\"api1\",\"api2\",\n]"
		    );
	    if ( node->n == 0) {
		warnings("Empty module configure. Skip ..");
		continue;
	    }
            struct module_config *module_config = &config->module;
            module_config->n_module = (int)node->n;
            module_config->files = (struct file_config*)malloc(module_config->n_module*sizeof(struct file_config));
            int n_files = 0;
            int j;
            for ( j = 0; j < node->n; ++j) {
                const kson_node_t *node1 = node->v.child[j];
                if ( node1 == NULL )
                    continue;
                
                struct file_config *file_config = &module_config->files[n_files];
                file_config->fname = BRANCH_INIT(node1);
                BRANCH(file_config->fname);
                file_config->columns = NULL;
                if ( file_config->fname == NULL )
                    continue;
                n_files++;
            }
        }
        else {
	    warnings("Unknown key : %s. skip ..", node->key);
	}
    }
#undef BRANCH    
#undef BRANCH_INIT
    return 0;
}

int bcfanno_load_config(struct bcfanno_config *config, const char *config_fname)
{
    // char *string = skip_comments(config_fname);
    char *dir = get_dirpath(config_fname);
    char *string = json_config_open(config_fname);
    if (string == NULL)
	error("Failed to parse configure file %s", config_fname);

    kson_t *json = NULL;
    json = kson_parse(string);
    free(string);
    int ret = load_config_core(config, json, dir);
    if (dir) free(dir);
    kson_destroy(json);
    return ret;
}

int bcfanno_config_debug(struct bcfanno_config *config)
{
    int i;
    LOG_print("configure file writer : %s", config->author == NULL ? "Unknown" : config->author);
    LOG_print("configure file ID : %s", config->config_id == NULL ? "Unknown" : config->config_id);
    LOG_print("reference sequence : %s", config->reference_path == NULL ? "Unknown" : config->reference_path);    

    if ( config->refgene.refgene_is_set == 1) {
	struct refgene_config *refgene = &config->refgene;	
	LOG_print("[refgene] GenePredPlus database: %s", refgene->genepred_fname);	
        LOG_print("[refgene] columns : %s", refgene->columns);
	if ( refgene->refseq_fname )
	    LOG_print("[refgene] transcript fasta : %s", refgene->refseq_fname);
	if ( refgene->trans_list_fname )
	    LOG_print("[refgene] trans_list : %s", refgene->trans_list_fname);
	if ( refgene->gene_list_fname )
	    LOG_print("[refgene] gene_list : %s", refgene->gene_list_fname);	
    }
    
    for ( i = 0; i < config->vcf.n_vcf; ++i ) {
	LOG_print("[vcf] %d", i);
	LOG_print("[vcf] file : %s", config->vcf.files[i].fname);
	LOG_print("[vcf] columns : %s", config->vcf.files[i].columns);	    
    }
    
    for ( i = 0; i < config->bed.n_bed; ++i ) {	
	LOG_print("[bed] %d", i);
	LOG_print("[bed] file : %s", config->bed.files[i].fname);
	if (config->bed.files[i].columns != NULL)
	    LOG_print("[bed] columns : %s", config->bed.files[i].columns);	    
    }
    return 0;
}


#ifdef _MAIN_CONFIG
int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "%s input.json\n", argv[0]);
	return 1;
    }
    struct bcfanno_config *con = bcfanno_config_init();
    bcfanno_load_config(con, argv[1]);
    bcfanno_config_debug(con);
    bcfanno_config_destroy(con);
    return 0;
}
#endif
