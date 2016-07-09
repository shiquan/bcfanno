#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "utils.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/khash.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "htslib/faidx.h"


KHASH_MAP_INIT_STR(list, char*)

typedef khash_t(list) glist_t;
		      
struct genepred {
    int rid;
    uint32_t txstart;
    uint32_t txend;
    char strand;
    char *name;
    char *name1; // transcript name or null for some genepred file
    uint32_t cdsstart;
    uint32_t cdsend;
    uint32_t exoncount;
    uint32_t *exonstarts;
    uint32_t *exonends;
};

struct genepred_mempool {
    int rid;
    uint32_t start;
    uint32_t end;
    int n, m;
    struct genepred *pools;    
};

struct args {
    char *genepred_file;
    tbx_t *gpidx;
    faidx_t *faidx;
    glist_t *glist;
};


static glist_t *glist = NULL;

int init_gene_list(const char *list)
{
    int i, n = 0;
    char **names = hts_readlist(list, 1, &n);
    if (n== 0) return 0;
    if (glist == NULL) glist = kh_init(list);
    
}
int usage(int argc, char **argv)
{
    fprintf(stderr, "Usage: %s [-h] -data genepred.txt.gz|refgene.txt.gz -O [u|v|z|b] -o output_file in_file.vcf.gz \n", argv[0]);
    return 0;
}

const char *hts_bcf_wmode(int file_type)
{
    if ( file_type == FT_BCF ) return "wbu";    // uncompressed BCF
    if ( file_type & FT_BCF ) return "wb";      // compressed BCF
    if ( file_type & FT_GZ ) return "wz";       // compressed VCF
    return "w";                                 // uncompressed VCF
}
void anno_hgvs_core(bcf_hdr_t *hdr, bcf1_t *line)
{
    
}
int main(int argc, char **argv)
{
    const char *in_file = 0; // stdin in default
    const char *out_file = 0; // stdout in default
    const char *out_type = 0; // output format, u in default
    const char *gp_file = 0; // error abort if not set genepred file

    for (int i = 1; i< argc;) {
	const char *a = argv[i++];
	if (strcmp(a, "-h") == 0) 
	    return usage(argc, argv);

	const char **arg_var = 0;
	if (strcmp(a, "-o") == 0 && out_file == 0)
	    arg_var = &out_file;
	else if (strcmp(a, "-O") == 0 && out_type == 0)
	    arg_var = &out_type;
	else if (strcmp(a, "-data") == 0 && gp_file == 0)
	    arg_var = &gp_file;
	
	if (arg_var != 0) {
	    if (i == argc)
		error("Missing arg after %s", a);
	    *arg_var = argv[i++];
	    continue;
	}

	if (in_file == 0) {
	    in_file = a;
	    continue;
	}
	error("Unknown arg : %s", a);
    }
    /* assuming input file is stdin, only vcf, gzipped vcf, bcf file is accept,
       err msg will output if file type unrecongnized*/
    if (in_file == 0 && (!isatty(fileno(stdin)))) in_file = "-";

    htsFile *fp = hts_open(in_file, "r");
    if (fp == NULL)
	error("Failed to open %s.", in_file);
    htsFormat type = *hts_get_format(fp);
    if (type->format  != vcf && type->format != bcf)
	error("Unsupported input format! %s", input_fname);
    
    const bcf_hdr_t *hdr = bcf_hdr_read(fp);
    
    /* duplicate header file and sync new tag in the output header .
     * assuming output is stdout in Vcf format in default. */
    bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);	
    int type = FT_VCF;
    if (out_type != 0) {
	switch (out_type[0]) {
	    case 'b':
		type = FT_BCF_GZ; break;
	    case 'u':
		type = FT_BCF; break;
	    case 'z':
		type = FT_VCF_GZ; break;
	    case 'v':
		type = FT_VCF; break;
	    default :
		error("The output type \"%s\" not recognised\n", out_type);
	};
    }
    htsFile *fout = out_file == 0 ? hts_open("-", hts_bcf_wmode(type)) : hts_open(out_file, hts_bcf_wmode(type));
    
    /* init genepred or refgene database, hold tabix index cache and faidx cache in memory */
    bcf1_t *line = bcf_init();
    while ( bcf_read(fp, hdr, line) ) {
	anno_hgvs_core(line);
	bcf_write1(fout, hdr_out, line);
    }
    
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    bcf_hdr_destroy(hdr_out);
    hts_close(fout);
    hts_close(fp);
    return 0;
}
