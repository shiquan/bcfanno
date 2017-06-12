// prase configure file in JSON format
// todo: stable improvement

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

struct vcfanno_config temp_config_init = {
    .author = 0,
    .config_id = 0,
    .reference_version = 0,
    .vcfs = { 0, 0 },
    .beds = { 0, 0 },
    .refgene = { 0, 0, 0, 0, 0},
};

struct vcfanno_config *vcfanno_config_init()
{
    struct vcfanno_config *config = (struct vcfanno_config*)malloc(sizeof(struct vcfanno_config));
    memset(config, 0, sizeof(struct vcfanno_config));
    return config;
}

void vcfanno_config_destroy(struct vcfanno_config *config)
{
    if ( config->author )
	free (config->author );
    if ( config->config_id )
	free (config->config_id );
    if ( config->reference_version )
	free (config->reference_version);
    int i;
    for (i = 0; i < config->vcfs.n_vcfs; ++i ) {
	free(config->vcfs.files[i].fname);
	free(config->vcfs.files[i].columns);
    }
    if ( i )
	free(config->vcfs.files);
    for (i = 0; i < config->beds.n_beds; ++i ) {
	free(config->beds.files[i].fname);
	if (config->beds.files[i].columns)
	    free(config->beds.files[i].columns);
    }
    if ( i )
	free(config->beds.files);

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


static int load_config_core(struct vcfanno_config *config, kson_t *json)
{
    const kson_node_t *root = json->root;
    if (root == NULL)
	error("Format error. Root node is empty, check configure file.");
    //kson_format(root);

#define BRANCH_INIT(_node) ( (_node)->v.str == NULL ? NULL : strdup((_node)->v.str) )
    int i;
    for (i = 0; i < root->n; ++i ) {
	const kson_node_t *node = kson_by_index(root, i);
	if ( node == NULL )
	    continue;
	if ( node->key == NULL) {
	    warnings("Foramt error. Node key is empty. skip..");
	    continue;
	}
	// the summary informations should in the top level
	if ( strcmp(node->key, "Author") == 0 || strcmp(node->key, "author") == 0) {
	    config->author = BRANCH_INIT(node);
	} else if ( strcmp(node->key, "ID") == 0 || strcmp(node->key, "id") == 0) {
	    config->config_id = BRANCH_INIT(node);
	} else if ( strcmp(node->key, "ref") == 0 || strcmp(node->key, "reference") == 0) {
	    config->reference_version = BRANCH_INIT(node);
	} else if ( strcmp(node->key, "HGVS") == 0 || strcmp(node->key, "hgvs") == 0) {
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

		if ( strcmp(node1->key, "gene_data") == 0 || strcmp(node1->key, "refgene") == 0)
		    refgene_config->genepred_fname = BRANCH_INIT(node1);
		else if ( strcmp(node1->key, "refseq") == 0 )
		    refgene_config->refseq_fname = BRANCH_INIT(node1);
		else if ( strcmp(node1->key, "transcripts_list") == 0 || strcmp(node1->key, "trans_list") == 0 )
		    refgene_config->trans_list_fname = BRANCH_INIT(node1);
		else if ( strcmp(node1->key, "genes_list") == 0 )
		    refgene_config->gene_list_fname = BRANCH_INIT(node1);
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
	} else if ( strcmp(node->key, "vcfs") == 0) {
	    if ( node->type != KSON_TYPE_BRACKET )
		error("Format error, configure for vcf databases should looks like :\n"
		      "\"vcfs\":[\n{\n\"file\":\"file.vcf.gz\",\"columns\":\"TAGS,TAGS\",\n},\n{},\n]"
		    );
	    if ( node->n == 0) {
		warnings("Empty vcfs configure. Skip");
		continue;
	    }
	    struct vcfs_config *vcfs_config = &config->vcfs;
	    vcfs_config->n_vcfs = (int)node->n;
	    vcfs_config->files = (struct file_config*)malloc(vcfs_config->n_vcfs*sizeof(struct file_config));
	    int n_files = 0;
	    int j;

	    for ( j = 0; j < (long)node->n; ++j ) {
		const kson_node_t *node1 = node->v.child[j];
		if ( node1 == NULL )
		    continue;
		int k;
		struct file_config *file_config = &vcfs_config->files[n_files];
		file_config->fname = NULL;
		file_config->columns = NULL;
		for ( k = 0; k < node1->n; ++k ) {
		    const kson_node_t *node2 = kson_by_index(node1, k);
		    if ( node2 == NULL || node2->key == NULL)
			continue;

		    if ( strcmp(node2->key, "file") == 0 )
			file_config->fname = BRANCH_INIT(node2);
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
	    if (n_files == 0 && vcfs_config->n_vcfs)
		free(vcfs_config->files);
	    vcfs_config->n_vcfs = n_files;
	} else if ( strcmp(node->key, "beds") == 0 ) {
	    if ( node->type != KSON_TYPE_BRACKET )
		error("Format error, configure for vcf databases should looks like :\n"
		      "\"beds\":[\n{\n\"file\":\"file.bed.gz\",\"header\":\"TAGS,TAGS\",\n},\n{},\n]"
		    );
	    if ( node->n == 0) {
		warnings("Empty beds configure. Skip ..");
		continue;
	    }
	    struct beds_config *beds_config = &config->beds;
	    beds_config->n_beds = (int)node->n;
	    beds_config->files = (struct file_config*)malloc(beds_config->n_beds*sizeof(struct file_config));
	    int n_files = 0;
	    int j;

	    for ( j = 0; j < (long)node->n; ++j ) {
		const kson_node_t *node1 = node->v.child[j];
		if ( node1 == NULL )
		    continue;
		int k;
		struct file_config *file_config = &beds_config->files[n_files];
		file_config->fname = NULL;
		file_config->columns = NULL;
		for ( k = 0; k < node1->n; ++k ) {
		    const kson_node_t *node2 = kson_by_index(node1, k);
		    if (node2 == NULL || node2->key == NULL)
			continue;
		    if ( strcmp(node2->key, "file") == 0 )
			file_config->fname = BRANCH_INIT(node2);
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
	    if ( n_files == 0 && beds_config->n_beds )
		free(beds_config->files);
	    beds_config->n_beds = n_files;	    
	} else {
	    warnings("Unknown key : %s. skip ..", node->key);
	}
    }
#undef BRANCH_INIT
    return 0;
}

int vcfanno_load_config(struct vcfanno_config *config, const char * config_fname)
{
    // char *string = skip_comments(config_fname);
    char *string = json_config_open(config_fname);
    if (string == NULL)
	error("Failed to parse configure file %s", config_fname);

    kson_t *json = NULL;
    json = kson_parse(string);
    free(string);
    int ret = load_config_core(config, json);
    kson_destroy(json);
    return ret;
}

int vcfanno_config_debug(struct vcfanno_config *config)
{
    int i;
    LOG_print("configure file writer : %s", config->author == NULL ? "Unknown" : config->author);
    LOG_print("configure file ID : %s", config->config_id == NULL ? "Unknown" : config->config_id);
    LOG_print("reference sequence : %s", config->reference_version == NULL ? "Unknown" : config->reference_version);    

    if ( config->refgene.refgene_is_set == 1) {
	struct refgene_config *refgene = &config->refgene;	
	LOG_print("[refgene] gene_data : %s", refgene->genepred_fname);	
	// LOG_print("[refgene] columns : %s", refgene->columns);
	if ( refgene->refseq_fname )
	    LOG_print("[refgene] refseq : %s", refgene->refseq_fname);
	if ( refgene->trans_list_fname )
	    LOG_print("[refgene] trans_list : %s", refgene->trans_list_fname);
	if ( refgene->gene_list_fname )
	    LOG_print("[refgene] gene_list : %s", refgene->gene_list_fname);	
    }
    
    for ( i = 0; i < config->vcfs.n_vcfs; ++i ) {
	LOG_print("[vcfs] %d", i);
	LOG_print("[vcfs] file : %s", config->vcfs.files[i].fname);
	LOG_print("[vcfs] columns : %s", config->vcfs.files[i].columns);	    
    }
    
    for ( i = 0; i < config->beds.n_beds; ++i ) {	
	LOG_print("[beds] %d", i);
	LOG_print("[beds] file : %s", config->beds.files[i].fname);
	if (config->beds.files[i].columns != NULL)
	    LOG_print("[beds] columns : %s", config->beds.files[i].columns);	    
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
    struct vcfanno_config *con = vcfanno_config_init();
    vcfanno_load_config(con, argv[1]);
    vcfanno_config_debug(con);
    vcfanno_config_destroy(con);
    return 0;
}
#endif
