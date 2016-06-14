#include "anno.h"
#include "vcmp.h"
#include "plugin.h"

// only if annotation database is VCF/BCF file, header_in has values or else header_in == NULL
annot_col_t * init_columns(char *rules, bcf_hdr_t *header_in, bcf_hdr_t *header_out, int *ncols, enum anno_type type)
{
    assert(rules != NULL);
    if (type == anno_file_is_vcf && header_in == NULL) {
	error("Inconsistent file type!");
    }
    char *ss = rules, *se = ss;
    *ncols = 0;
    annot_col_t *cols = NULL;
    kstring_t tmp = KSTRING_INIT;
    kstring_t str = KSTRING_INIT; 
    int i = -1;
    while (*ss)
    {
	if ( *se && *se!=',' ) {
	    se++;
	    continue;
	}
	int replace = REPLACE_ALL;
	if ( *ss=='+') {
	    replace == REPLACE_MISSING;
	    ss++;
	} else if (*ss=='-') {
	    replace == REPLACE_EXISTING;
	    ss++;
	}
	i++;
	str.l = 0;
	kputsn(ss, se-ss, &str); 
	if ( !str.s[0] ); // empty skip
	else if ( !strcasecmp("CHROM", str.s) || !strcasecmp("POS", str.s) ||
                  !strcasecmp("FROM", str.s) || !strcasecmp("TO", str.s) ||
                  !strcasecmp("REF", str.s) || !strcasecmp("ALT", str.s) ||
                  !strcasecmp("FILTER", str.s) || !strcasecmp("QUAL", str.s)
            ); // skip, for consistent with vcfannotate.c only
        else if ( !strcasecmp("ID", str.s) ) {
	    *ncols++;
            cols = (struct annot_col*) realloc(cols, sizeof(struct annot_col)* (*ncols));
            struct annot_col *col = &cols[*ncols-1];
            col->icol = i;
            col->replace = replace;
            col->setter = type == anno_file_is_vcf ? vcf_setter_id : anno_setter_id;
            col->hdr_key = strdup(str.s);
        } else if (!strcasecmp("INFO", str.s) || !strcasecmp("FORMAT", str.s) ) {
	    error("do not support annotate all INFO/FORMAT fields. todo INFO/TAG instead\n");
	} else if (!strncasecmp("FORMAT/", str.s, 7) || !strncasecmp("FMT/", str.s, 4)) {
            char *key = str.s + (!strncasecmp("FMT", str.s, 4) ? 4 : 7);
            if (!strcasecmp("GT", key) ) {
		error("It is not allowed to change GT tag.");
	    }
            //bcf_hrec_t *hrec = bcf_hdr_get_hrec(acols.header, BCF_HL_FMT, "ID", key, NULL);
            //tmp.l = 0;
            //bcf_hrec_format(hrec, &tmp);
            //bcf_hdr_append(header_out, tmp.s);
            //bcf_hdr_sync(header_out);
	    int hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, str.s);
	    if ( !bcf_hdr_idinfo_exists(header_out, BCF_HL_FMT, hdr_id) ) {
		if ( type == anno_file_is_vcf ) {
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

            int hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, key);
            *ncols++;
	    cols = (struct annot_col*) realloc(cols, sizeof(struct annot_col)*(*ncols));
            struct annot_col *col = &cols[*ncols-1];
            col->icol = -1;
            col->replace = replace;
            col->hdr_key = strdup(key);
            switch ( bcf_hdr_id2type(header_out, BCF_HL_FMT, hdr_id) )
	    {
		case BCF_HT_INT:
		    col->setter = acols->type == anno_file_is_vcf ? vcf_setter_format_int ? anno_setter_format_int;
		    break;
                case BCF_HT_REAL:
		    col->setter = acols->type == anno_file_is_vcf ? vcf_setter_format_real ? anno_setter_format_real;
		    break;
                case BCF_HT_STR:
		    col->setter = acols->type == anno_file_is_vcf ? vcf_setter_format_str ? anno_setter_format_str;
		    break;
                default :
		    error("The type of %s not recognised (%d)\n", str.s, bcf_hdr_id2type(args->hdr_out, BCF_HL_FMT, hdr_id));
		    
            }
	} else if (!strncasecmp("INFO/RefGene", str.s, 12) || !strncasecmp("RefGene", str.s, 7)) {
	    // RefGene.FuncRegion      RefGene.ExIn   RefGene.Function RefGene.Gene RefGene.HGVS
	    // All tags with RefGene started name will construct from one or more transcripts and be splited by '|'
	    if ( !strncasecmp("INFO/", str.s, 5) ) {
		memmove(str.s, str.s+5, str.l-4);
		str.l -= 4;
	    }
	    int hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, str.s);

	    if ( !strcasecmp("RefGene.HGVS", str.s, 12)) {
		if ( !bcf_hdr_idinfo_exists(header_out, BCF_HL_INFO, hdr_id) ) {
		    bcf_hdr_append("##INFO=<ID=RefGene.HGVS,Number=A,Type=String,Description=\"HGVS nomenclature annotated by VCFANNO\">");
		    bcf_hdr_sync(header_out);
		    hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, str.s);
		    assert( bcf_hdr_idinfo_exists(header_out, BCF_HL_INFO, hdr_id) );
		}
		*ncols++;
		cols = (struct annot_col*) realloc(cols, sizeof(struct annot_col)*(*ncols));
		struct annot_col *col = &cols[*ncols-1];
		col->icol = i;
		col->replace = replace;
		col->hdr_key = strdup(str.s);
		col->number = bcf_hdr_id2length(header_out, BCF_HL_INFO, hdr_id);
		col->setter = local_setter_refgene_hgvs;
	    } else if ( !strcasecmp("RefGene.FuncRegion", str.s, 17)) {
		if ( !bcf_hdr_idinfo_exists(header_out, BCF_HL_INFO, hdr_id) ) {
		    bcf_hdr_append("##INFO=<ID=RefGene.FuncRegion,Number=A,Type=String,Description=\"Function region annotated by VCFANNO\">");
		    bcf_hdr_sync(header_out);
		    hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, str.s);
		    assert( bcf_hdr_idinfo_exists(header_out, BCF_HL_INFO, hdr_id) );
		}
		*ncols++;
		cols = (struct annot_col*) realloc(cols, sizeof(struct annot_col)*(*ncols));
		struct annot_col *col = &cols[*ncols-1];
		col->icol = i;
		col->replace = replace;
		col->hdr_key = strdup(str.s);
		col->number = bcf_hdr_id2length(header_out, BCF_HL_INFO, hdr_id);
		col->setter = local_setter_refgene_funcreg;
	    } else if ( !strcasecmp("RefGene.ExIn", str.s, 12)) {
		if ( !bcf_hdr_idinfo_exists(header_out, BCF_HL_INFO, hdr_id) ) {
		    bcf_hdr_append("##INFO=<ID=RefGene.FuncRegion,Number=A,Type=String,Description=\"Exon or intron id annotated by VCFANNO\">");
		    bcf_hdr_sync(header_out);
		    hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, str.s);
		    assert( bcf_hdr_idinfo_exists(header_out, BCF_HL_INFO, hdr_id) );
		}
		*ncols++;
		cols = (struct annot_col*) realloc(cols, sizeof(struct annot_col)*(*ncols));
		struct annot_col *col = &cols[*ncols-1];
		col->icol = i;
		col->replace = replace;
		col->hdr_key = strdup(str.s);
		col->number = bcf_hdr_id2length(header_out, BCF_HL_INFO, hdr_id);
		col->setter = local_setter_refgene_exin;
	    } else if ( !strcasecmp("RefGene.Fuction", str.s, 15)) {

		if ( !bcf_hdr_idinfo_exists(header_out, BCF_HL_INFO, hdr_id) ) {
		    bcf_hdr_append("##INFO=<ID=RefGene.FuncRegion,Number=A,Type=String,Description=\"Function annotated by VCFANNO\">");
		    bcf_hdr_sync(header_out);
		    hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, str.s);
		    assert( bcf_hdr_idinfo_exists(header_out, BCF_HL_INFO, hdr_id) );
		}
		*ncols++;
		cols = (struct annot_col*) realloc(cols, sizeof(struct annot_col)*(*ncols));
		struct annot_col *col = &cols[*ncols-1];
		col->icol = i;
		col->replace = replace;
		col->hdr_key = strdup(str.s);
		col->number = bcf_hdr_id2length(header_out, BCF_HL_INFO, hdr_id);
		col->setter = local_setter_refgene_function;	      
	    } else {
		error("The tag %s not support yet (send an email to the mail list for function request please.)", str.s);
	    }
	} else if ( !strcasecmp("Gene", str.s, 4)) {
	    // Usually a mutation localed in one gene
	    if ( !bcf_hdr_idinfo_exists(header_out, BCF_HL_INFO, hdr_id) ) {
		bcf_hdr_append("##INFO=<ID=RefGene.Gene,Number=1,Type=String,Description=\"Gene name annotated by VCFANNO\">");
		bcf_hdr_sync(header_out);
		hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, str.s);
		assert( bcf_hdr_idinfo_exists(header_out, BCF_HL_INFO, hdr_id) );
	    }
		
	    *ncols++;
            cols = (struct annot_col*) realloc(cols, sizeof(struct annot_col)*(*ncols));
            struct annot_col *col = &cols[*ncols-1];
            col->icol = i;
            col->replace = replace;
            col->hdr_key = strdup(str.s);
            col->number = bcf_hdr_id2length(header_out, BCF_HL_INFO, hdr_id);
	    col->setter = local_setter_refgene_gene;
	} else {
	    if ( !strncasecmp("INFO/", str.s, 5) ) { memmove(str.s, str.s+5, str.l-4); }
	    int hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, str.s);
	    if ( !bcf_hdr_idinfo_exists(header_out, BCF_HL_INFO, hdr_id) ) {
		if ( type == anno_file_is_vcf ) {
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
            *ncols++;
            cols = (struct annot_col*) realloc(cols, sizeof(struct annot_col)*(*ncols));
            struct annot_col *col = &cols[*ncols-1];
            col->icol = i;
            col->replace = replace;
            col->hdr_key = strdup(str.s);
            col->number = bcf_hdr_id2length(header_out, BCF_HL_INFO, hdr_id);
            switch ( bcf_hdr_id2type(header_out, BCF_HL_INFO, hdr_id) )
            {
                case BCF_HT_FLAG:  col->setter = acols->type == anno_file_is_vcf ? vcf_setter_info_flag : anno_setter_info_flag; break;
                case BCF_HT_INT:   col->setter = acols->type == anno_file_is_vcf ? vcf_setter_info_int : anno_setter_info_int; break;
                case BCF_HT_REAL:  col->setter = acols->type == anno_file_is_vcf ? vcf_setter_info_real : anno_setter_info_real; break;
                case BCF_HT_STR:   col->setter = acols->type == anno_file_is_vcf ? vcf_setter_info_str : anno_setter_info_str; break;
                default: error("The type of %s not recognised (%d)\n", str.s, bcf_hdr_id2type(args->hdr_out, BCF_HL_INFO, hdr_id));
            }
        }
        if ( !*se ) break;
        ss = ++se;
    }
    if (str.m) free(str.s);
    if (tmp.m) free(tmp.s);
    return cols;
}


/* int anno_setter_id(anno_setters_t *handler, bcf1_t *line, struct annot_col *col, void *data) */
/* { */
/*     annot_line_t *tab = (annot_line_t*)data; */
/*     if ( tab->cols[col->icol] && tab->cols[col->icol][0]=='.' && !tab->cols[col->icol][1]) */
/* 	return 0; // donot replace with '.' */
/*     if ( col->replace != REPLACE_MISSING) */
/* 	return bcf_update_id(handler->hdr_out, line, tab->cols[col->icol]); */

/*     if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) ) */
/* 	return bcf_update_id() */
/* } */
