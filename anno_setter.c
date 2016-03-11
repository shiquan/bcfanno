#include "anno.h"
#include "vcmp.h"



struct annot_cols_pack {
    annot_col *cols;
    int ncols;
    int type;
    char *columns;
};

struct annot_cols_pack * annopacks;

void init_columns(struct annot_cols_pack *acols, char *rules, bcf_hdr_t *header_out)
{
    if ( rules==NULL) return;
    acols->columns = strdup(rules);
    char *ss = rules, *se = ss;
    cols->n = 0;
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
	kputsn(ss, se-ss, &str); // empty skip
	if ( !str.s[0] );
	else if ( !strcasecmp("CHROM", str.s) || !strcasecmp("POS", str.s) ||
                  !strcasecmp("FROM", str.s) || !strcasecmp("TO", str.s) ||
                  !strcasecmp("REF", str.s) || !strcasecmp("ALT", str.s) ||
                  !strcasecmp("FILTER", str.s) || !strcasecmp("QUAL", str.s)
            ); // skip, for consistent with vcfannotate.c only
        else if ( !strcasecmp("ID", str.s) ) {
	    acols->ncols++;
            acols->cols = (struct annot_col*) realloc(acols->cols, sizeof(struct annot_col*)*acols->ncols);
            struct annot_col *col = &acols->cols[acols->ncols-1];
            col->icol = i;
            col->replace = replace;
            col->setter = acols->type == anno_file_is_vcf ? vcf_setter_id : anno_setter_id;
            col->hdr_key = strdup(str.s);
        } else if (!strcasecmp("INFO", str.s) || !strcasecmp("FORMAT", str.s) ) {
	    error("do not support annotate all INFO/FORMAT fields. todo INFO/TAG instead\n");
	} else if (!strncasecmp("FORMAT/", str.s, 7) || !strncasecmp("FMT/", str.s, 4) {
            char *key = str.s + (!strncasecmp("FMT", str.s, 4) ? 4 : 7);
            if (!strcasecmp("GT", key) ) {
		error("It is not allowed to change GT tag.");
	    }
            bcf_hrec_t *hrec = bcf_hdr_get_hrec(acols.header, BCF_HL_FMT, "ID", key, NULL);
            tmp.l = 0;
            bcf_hrec_format(hrec, &tmp);
            bcf_hdr_append(header_out, tmp.s);
            bcf_hdr_sync(header_out);
            int hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, key);
            acols->ncols++;
	    acols->cols = (struct annot_col*) realloc(acols->cols, sizeof(struct annot_col)*acols->ncols);
            struct annot_col *col = &acols->cols[acols->ncols-1];
            col->icol = -1;
            col->replace = replace;
            col->hdr_key = strdup(key);
            switch ( bcf_hdr_id2type(header_out, BCF_HL_FMT, hdr_id) ) {
		case BCF_HT_INT:  col->setter = acols->type == anno_file_is_vcf ? vcf_setter_format_int ? anno_setter_format_int; break;
                case BCF_HT_REAL: col->setter = acols->type == anno_file_is_vcf ? vcf_setter_format_real ? anno_setter_format_real; break;
                case BCF_HT_STR:  col->setter = acols->type == anno_file_is_vcf ? vcf_setter_format_str ? anno_setter_format_str; break;
                default : error("The type of %s not recognised (%d)\n", str.s, bcf_hdr_id2type(args->hdr_out, BCF_HL_FMT, hdr_id));
            }
        }
        else
        {
            if ( !strncasecmp("INFO/", str.s, 5) ) { memmove(str.s, str.s+5, str.l-4); }
            int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, str.s);
            if ( !bcf_hdr_idinfo_exists(args->hdr_out, BCF_HL_INFO, hdr_id) )
            {
                bcf_hrec_t *hrec = bcf_hdr_get_hrec(acols.header, BCF_HL_INFO, "ID", str.s, NULL);
                if ( !hrec ) error("The tag \"%s\" is not defined in %s\n", str.s, args->files->readers[1].fname);
                tmp.l = 0;
                bcf_hrec_format(hrec, &tmp);
                bcf_hdr_append(header_out, tmp.s);
                bcf_hdr_sync(header_out);
                hdr_id = bcf_hdr_id2int(header_out, BCF_DT_ID, str.s);
                assert( bcf_hdr_idinfo_exists(header_out, BCF_HL_INFO, hdr_id) );
            }

            acols->ncols++;
            acols->cols = (struct annot_col*) realloc(acols->cols, sizeof(struct annot_col)*acols->ncols);
            struct annot_col *col = &acols->cols[acols->ncols-1];
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
    free(str.s);
    free(tmp.s);
}


int anno_setter_id(anno_setters_t *handler, bcf1_t *line, struct annot_col *col, void *data)
{
    anno_line_t *tab = (anno_line_t*)data;
    if ( tab->cols[col->icol] && tab->cols[col->icol][0]=='.' && !tab->cols[col->icol][1])
	return 0; // donot replace with '.'
    if ( col->replace != REPLACE_MISSING)
	return bcf_update_id(handler->hdr_out, line, tab->cols[col->icol]);

    if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) )
	return bcf_update_id()
}
