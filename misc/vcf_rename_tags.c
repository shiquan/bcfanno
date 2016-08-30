// rename tags for vcf/bcf file, only support INFO fields
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "utils.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#define KSTRING_INIT {0, 0, 0}

const char *list_fname = 0;

struct args  {
    htsFile *fp_input;
    htsFile *fp_output;
    bcf_hdr_t *hdr_in;
    bcf_hdr_t *hdr_out;
};

struct args args = {
    .fp_input = NULL,
    .fp_output = NULL,
    .hdr_in = NULL,
    .hdr_out = NULL,
};

void args_destroy() {
    hts_close(args.fp_input);
    hts_close(args.fp_output);
    bcf_hdr_destroy(args.hdr_in);
    bcf_hdr_destroy(args.hdr_out);
}
// tags list with two columns, old name -> new name

int usage() {
    fprintf(stderr, "vcf_renames_tags -list tags.txt -O <type> -o output_file input_file\n");
    return 1;
}
static const char *hts_bcf_wmode(int file_type)
{
    if ( file_type == FT_BCF )
	return "wbu";    // uncompressed BCF
    if ( file_type & FT_BCF )
	return "wb";      // compressed BCF
    if ( file_type & FT_GZ )
	return "wz";       // compressed VCF
    return "w";                                 // uncompressed VCF
}

int parse_args(int argc, char **argv)
{
    if (argc == 0)
	error("Failed to parse arguments. Use -h for more information.");
    const char *input_fname = 0;
    const char *output_fname = 0;
    const char *out_type_string = 0;
    int i;
    for ( i = 0; i < argc; ) {
	const char *a = argv[i++];
	if ( strcmp(a, "-h") == 0 )
	    return usage();
	const char **arg_var = 0;
	if ( strcmp(a, "-list") == 0 )
	    arg_var = &list_fname;
	else if ( strcmp(a, "-O") == 0 && out_type_string == 0 )
	    arg_var = &out_type_string;
	else if ( strcmp(a, "-o") == 0 && output_fname == 0 )
	    arg_var = &output_fname;

	if ( arg_var != 0 ) {
	    if ( i == argc )
		error("Missing arg after %s.", a);
	    *arg_var = argv[i++];
	    continue;
	}
	if ( input_fname == 0 ) {
	    input_fname = a;
	    continue;
	}
	error("Unknown arg %s.", a);
    }

    if ( input_fname == 0 && (!isatty(fileno(stdin))) )
	input_fname = "-";
    if ( input_fname == 0 )
	error("No input file!");
    if ( list_fname == 0 )
	error("No list file specified.");

    args.fp_input = hts_open(input_fname, "r");
    if (args.fp_input == NULL)
	error("Failed to open %s : %s.", input_fname, strerror(errno));

    htsFormat type = *hts_get_format(args.fp_input);
    if (type.format  != vcf && type.format != bcf)
	error("Unsupported input format! %s", input_fname);

    int out_type = FT_VCF;
    if (out_type_string != 0) {
	switch (out_type_string[0]) {
	    case 'b':
		out_type = FT_BCF_GZ; break;
	    case 'u':
		out_type = FT_BCF; break;
	    case 'z':
		out_type = FT_VCF_GZ; break;
	    case 'v':
		out_type = FT_VCF; break;
	    default :
		error("The output type \"%d\" not recognised\n", out_type);
	};
    }
    args.hdr_in = bcf_hdr_read(args.fp_input);
    htsFile *fout = output_fname == 0 ? hts_open("-", hts_bcf_wmode(out_type)) : hts_open(output_fname, hts_bcf_wmode(out_type));
    args.hdr_out = bcf_hdr_dup(args.hdr_in);
    
    return 0;
}

int tags_convert()
{
    if (list_fname == NULL)
	error("No list file specified. Cannot init.");
    int n, i;
    char **names = hts_readlist(list_fname, 1, &n);
    if (n == 0)
	return 1;

    kstring_t string = KSTRING_INIT;
    int nfields = 0;
    for ( i = 0; i < n; ++i ) {
	// skip comment lines
	if (names[i] == NULL || names[i][0] == '#')
	    continue;
	// init string
	string.l = 0;
	kputs(names[i], &string);
	// split string by tab
	int *splits = ksplit(&string, '\t', &nfields);
	if (nfields != 2) {
	    warnings("Error format; only accept two columns per line. %s", names[i]);
	    continue;
	}
	char *ss = string.s + splits[0];
	char *se = string.s + splits[1];
	if ( strlen(ss) > 4 && strncmp(ss, "CTG/", 4) == 0 ) {
	    ss += 4;
	    if ( strlen(se) > 4 && strncmp(se, "CTG/", 4) == 0 )
		ss += 4;
	    int rid = bcf_hdr_name2id(args.hdr_out, ss);
	    bcf_hrec_t *rec = bcf_hdr_get_hrec(args.hdr_out, BCF_HL_CTG, "ID", ss, NULL);
	    if ( !rec )
		continue;
	    int j = bcf_hrec_find_key(rec, "ID");
	    assert(j >= 0);
	    free(rec->vals[j]);
	    rec->vals[j] = strdup(se);
	    args.hdr_out->id[BCF_DT_CTG][rid].key = rec->vals[j];
	    free(splits);
	    continue;
	}
	    
	if ( strlen(ss) > 5 && strncmp(ss, "INFO/", 5) == 0 )
	    ss += 5;
	if ( strlen(se) > 5 && strncmp(se, "INFO/", 5) == 0 )
	    se += 5;
	int rid = bcf_hdr_id2int(args.hdr_out, BCF_DT_ID, ss);
	bcf_hrec_t *rec = bcf_hdr_get_hrec(args.hdr_out, BCF_HL_INFO, "ID", ss, NULL);
	if ( !rec )
	    continue;
	int j = bcf_hrec_find_key(rec, "ID");
	assert(j >= 0);
	free(rec->vals[j]);
	rec->vals[j] = strdup(se);
	args.hdr_out->id[BCF_DT_ID][rid].key = rec->vals[j];	
	free(splits);
    }
    for ( i = 0; i < n; ++i )
	free(names[i]);
    free(names);
    if ( string.m )
	free(string.s);
    return 0;
}
int main(int argc, char **argv)
{
    if ( parse_args(--argc, ++argv) == 1)
	return 1;
    tags_convert();
    bcf1_t *line = bcf_init();
    while ( bcf_read(args.fp_input, args.hdr_in, line) == 0 ) {
	bcf_write1(args.fp_output, args.hdr_out, line);
    }
    bcf_destroy(line);
    args_destroy();
    return 0;
}
