#include "anno.h"
#include "vcmp.h"
#include "plugin.h"
//#include "hgvs.h"

// only if annotation database is VCF/BCF file, header_in has values or else header_in == NULL
annot_col_t *init_columns(const char *rules, bcf_hdr_t *header_in, bcf_hdr_t *header_out, int *ncols, enum anno_type type)
{
    assert(rules != NULL);
    if (type == anno_is_vcf && header_in == NULL) {
	error("Inconsistent file type!");
    }
    char *ss = (char*)rules, *se = ss;
    int nc = 0;
    annot_col_t *cols = NULL;
    kstring_t tmp = KSTRING_INIT;
    kstring_t str = KSTRING_INIT;
    int i = -1;

    while (*ss) {
	if ( *se && *se!=',' ) {
	    se++;
	    continue;
	}
	int replace = REPLACE_ALL;
	if ( *ss=='+') {
	    replace = REPLACE_MISSING;
	    ss++;
	} else if (*ss=='-') {
	    replace = REPLACE_EXISTING;
	    ss++;
	}
	i++;
	str.l = 0;
	kputsn(ss, se-ss, &str); 
	if ( !str.s[0] ) {
	    warnings("Empty tag in %s", rules);
	} else if ( !strcasecmp("CHROM", str.s) || !strcasecmp("POS", str.s) || !strcasecmp("FROM", str.s) || !strcasecmp("TO", str.s) || !strcasecmp("REF", str.s) || !strcasecmp("ALT", str.s) || !strcasecmp("FILTER", str.s) || !strcasecmp("QUAL", str.s)) {
	    warnings("Skip tag %s", str.s);
	} else if ( !strcasecmp("ID", str.s) ) {
	    nc++;
            cols = (struct annot_col*) realloc(cols, sizeof(struct annot_col)* (nc));
            struct annot_col *col = &cols[nc-1];
            col->icol = i;
            col->replace = replace;
            col->setter = type == anno_is_vcf ? vcf_setter_id : setter_id;
            col->hdr_key = strdup(str.s);
        } else if (!strcasecmp("INFO", str.s) || !strcasecmp("FORMAT", str.s) ) {
	    error("do not support annotate all INFO,FORMAT fields. todo INFO/TAG instead\n");
	} else if (!strncasecmp("FORMAT/", str.s, 7) || !strncasecmp("FMT/", str.s, 4)) {
            char *key = str.s + (!strncasecmp("FMT", str.s, 4) ? 4 : 7);
            if (!strcasecmp("GT", key)) 
		error("It is not allowed to change GT tag.");

	    int hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, str.s);
	    
	    if ( !bcf_hdr_idinfo_exists(header_out, BCF_HL_FMT, hdr_id) ) {
		
		if ( type == anno_is_vcf ) {
		    bcf_hrec_t *hrec = bcf_hdr_get_hrec(header_in, BCF_HL_FMT, "ID", str.s, NULL);
		    if ( !hrec )
			error("The tag \"%s\" is not defined in header: %s\n", str.s, rules);
		    tmp.l = 0;
		    bcf_hrec_format(hrec, &tmp);
		    bcf_hdr_append(header_out, tmp.s);
		    bcf_hdr_sync(header_out);
		    hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, str.s);
		    assert( bcf_hdr_idinfo_exists(header_out, BCF_HL_FMT, hdr_id) );
		} else {
		    error("The tag \"%s\" is not defined in header: %s\n", str.s, rules);
		}
	    }

            //int hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, key);
            nc++;
	    cols = (struct annot_col*) realloc(cols, sizeof(struct annot_col)*(nc));
            struct annot_col *col = &cols[nc-1];
            col->icol = -1;
            col->replace = replace;
            col->hdr_key = strdup(key);

            switch ( bcf_hdr_id2type(header_out, BCF_HL_FMT, hdr_id) ) {

		case BCF_HT_INT:
		    col->setter = type == anno_is_vcf ? vcf_setter_format_int : setter_format_int;
		    break;

		case BCF_HT_REAL:
		    col->setter = type == anno_is_vcf ? vcf_setter_format_real : setter_format_real;
		    break;

		case BCF_HT_STR:
		    col->setter = type == anno_is_vcf ? vcf_setter_format_str : setter_format_str;
		    break;

		default :
		    error("The type of %s not recognised (%d)\n", str.s, bcf_hdr_id2type(header_out, BCF_HL_FMT, hdr_id));
            }

	} else if ( !strncasecmp("INFO/", str.s, 5) ) {
	    memmove(str.s, str.s+5, str.l-4);
	    str.l -= 4;

	    int hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, str.s);

	    if ( !bcf_hdr_idinfo_exists(header_out, BCF_HL_INFO, hdr_id) ) {
		if ( type == anno_is_vcf ) {
		    bcf_hrec_t *hrec = bcf_hdr_get_hrec(header_in, BCF_HL_INFO, "ID", str.s, NULL);

		    if ( !hrec )
			error("The tag \"%s\" is not defined in header: %s\n", str.s, rules);
		    tmp.l = 0;

		    bcf_hrec_format(hrec, &tmp);
		    bcf_hdr_append(header_out, tmp.s);
		    bcf_hdr_sync(header_out);
		    hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, str.s);
		    assert( bcf_hdr_idinfo_exists(header_out, BCF_HL_INFO, hdr_id) );
		} else {
		    error("The tag \"%s\" is not defined in header: %s\n", str.s, rules);
		}
	    }
	    nc++;
	    cols = (struct annot_col*) realloc(cols, sizeof(struct annot_col)*(nc));
	    struct annot_col *col = &cols[nc-1];
	    col->icol = i;
	    col->replace = replace;
	    col->hdr_key = strdup(str.s);

	    col->number = bcf_hdr_id2length(header_out, BCF_HL_INFO, hdr_id);

	    switch ( bcf_hdr_id2type(header_out, BCF_HL_INFO, hdr_id) ) {

		case BCF_HT_FLAG:
		    col->setter = type == anno_is_vcf ? vcf_setter_info_flag : setter_info_flag;
		    break;

		case BCF_HT_INT:
		    col->setter = type == anno_is_vcf ? vcf_setter_info_int : setter_info_int;
		    break;

		case BCF_HT_REAL:
		    col->setter = type == anno_is_vcf ? vcf_setter_info_real : setter_info_real;
		    break;

		case BCF_HT_STR:
		    col->setter = type == anno_is_vcf ? vcf_setter_info_str : setter_info_str;
		    break;

		default:
		    error("The type of %s not recognised (%d)\n", str.s, bcf_hdr_id2type(header_out, BCF_HL_INFO, hdr_id));
	    }
	}
	if ( !*se ) break;
        ss = ++se;
    }
    *ncols = nc;
    if (str.m) free(str.s);
    if (tmp.m) free(tmp.s);
    return cols;
}

void print_annot_cols(annot_col_t *cols, int n)    
{
    int i;
    for (i=0; i<n; ++i) {
	annot_col_t *col = &cols[i];
	fprintf(stderr, "[%s] %d : %s, %p\n",__func__, col->icol, col->hdr_key, col->setter);
    }
    return;
}

#ifdef _SETTER_MAIN

#include <htslib/hts.h>
#include <htslib/vcf.h>

int main(int argc, char **argv)
{
    if (argc != 3) {
	fprintf(stderr,"anno_setter <in.vcf.gz> <columns_string>\n");
	return 1;
    }
    bcf_hdr_t *h = NULL; //bcf_hdr_init();
    htsFile *fp = hts_open(argv[1], "r");
    if (fp == NULL)
	error("%s : %s", argv[1], strerror(errno));
    
    h = bcf_hdr_read(fp);
    if (h == NULL)
	error("failed to prase header");
    bcf_hdr_t *out = bcf_hdr_dup(h);
    char *string = strdup(argv[2]);
    int ncols = 0;
    annot_col_t *cols = init_columns(string, h, out, &ncols, anno_is_vcf);
    print_annot_cols(cols, ncols);
    hts_close(fp);
    return 0;
}

#endif
