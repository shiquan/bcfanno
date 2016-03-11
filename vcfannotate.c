#include "anno.h"
#include "plugin.h"
#include "vcmp.h"
//#include "filter.h"
//#include "convert.h"

/* struct _args_t; */

/* Typedef struct _rm_tag_t */
/* { */
/*     char *key; */
/*     int hdr_id; */
/*     void (*handler)(struct _args_t *, bcf1_t *, struct _rm_tag_t *); */
/* } */
/* rm_tag_t; */

/* typedef struct */
/* { */
/*     char **cols; */
/*     int ncols, mcols; */
/*     char **als; */
/*     int nals, mals; */
/*     kstring_t line; */
/*     int rid, start, end; */
/* } */
/* annot_line_t; */

#define REPLACE_MISSING  0  // replace only missing values
#define REPLACE_ALL      1  // replace both missing and existing values
#define REPLACE_EXISTING 2  // replace only if tgt is not missing
#define SET_OR_APPEND    3  // set new value if missing or non-existent, append otherwise

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define MARK_LISTED   1
#define MARK_UNLISTED 2

/* typedef struct _hand_t */
/* { */
/*     bcf_srs_t *files; */
/*     bcf_hdr_t *hdr, *hdr_out; */
/*     htsFile *out_fh; */
/*     int output_type, n_threads; */
/*     //bcf_sr_regions_t *tgts; */

/*     filter_t *filter; */
/*     char *filter_str; */
/*     int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE */

/*     /\* rm_tag_t *rm;           // tags scheduled for removal *\/ */
/*     /\* int nrm; *\/ */
/*     /\* int flt_keep_pass;      // when all filters removed, reset to PASS *\/ */

/*     vcmp_t *vcmp;           // for matching annotation and VCF lines by allele */
/*     annot_line_t *alines;   // buffered annotation lines */
/*     int nalines, malines; */
/*     int ref_idx, alt_idx, chr_idx, from_idx, to_idx;   // -1 if not present */
/*     annot_col_t *cols;      // column indexes and setters */
/*     int ncols; */

/*     /\* char *set_ids_fmt; *\/ */
/*     /\* convert_t *set_ids; *\/ */
/*     /\* int set_ids_replace; *\/ */

/*     /\* int *sample_map, nsample_map, sample_is_file;   // map[idst] -> isrc *\/ */
/*     int mtmpi, mtmpf, mtmps; */
/*     int mtmpi2, mtmpf2, mtmps2; */
/*     int mtmpi3, mtmpf3, mtmps3; */
/*     int32_t *tmpi, *tmpi2, *tmpi3; */
/*     float *tmpf, *tmpf2, *tmpf3; */
/*     char *tmps, *tmps2, **tmpp, **tmpp2; */
/*     kstring_t tmpks; */

/*     /\* char **argv, *output_fname, *targets_fname, *regions_list, *header_fname; *\/ */
/*     /\* char *remove_annots, *columns, *rename_chrs, *sample_names, *mark_sites; *\/ */
/*     /\* int argc, drop_header, record_cmd_line, tgts_is_vcf, mark_sites_logic; *\/ */
/* } */
/* hand_t; */

//char *msprintf(const char *fmt, ...);

/*
void remove_id(hand_t *hand, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_id(hand->hdr,line,NULL);
}
void remove_filter(hand_t *hand, bcf1_t *line, rm_tag_t *tag)
{
    if ( !tag->key ) bcf_update_filter(hand->hdr, line, NULL, hand->flt_keep_pass);
    else bcf_remove_filter(hand->hdr, line, tag->hdr_id, hand->flt_keep_pass);
}
void remove_qual(hand_t *hand, bcf1_t *line, rm_tag_t *tag)
{
    bcf_float_set_missing(line->qual);
}
void remove_info(hand_t *hand, bcf1_t *line, rm_tag_t *tag)
{
    // remove all INFO fields
    if ( !(line->unpacked & BCF_UN_INFO) ) bcf_unpack(line, BCF_UN_INFO);

    int i;
    for (i=0; i<line->n_info; i++)
    {
        bcf_info_t *inf = &line->d.info[i];
        if ( inf->vptr_free )
        {
            free(inf->vptr - inf->vptr_off);
            inf->vptr_free = 0;
        }
        line->d.shared_dirty |= BCF1_DIRTY_INF;
        inf->vptr = NULL;
    }
}
void remove_info_tag(hand_t *hand, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_info(hand->hdr, line, tag->key, NULL, 0, BCF_HT_INT);  // the type does not matter with n=0
}
void remove_format_tag(hand_t *hand, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_format(hand->hdr, line, tag->key, NULL, 0, BCF_HT_INT);  // the type does not matter with n=0
}
void remove_format(hand_t *hand, bcf1_t *line, rm_tag_t *tag)
{
    // remove all FORMAT fields except GT
    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);

    int i;
    for (i=0; i<line->n_fmt; i++)
    {
        bcf_fmt_t *fmt = &line->d.fmt[i];
        const char *key = bcf_hdr_int2id(hand->hdr,BCF_DT_ID,fmt->id);
        if ( key[0]=='G' && key[1]=='T' && !key[2] ) continue;

        if ( fmt->p_free )
        {
            free(fmt->p - fmt->p_off);
            fmt->p_free = 0;
        }
        line->d.indiv_dirty = 1;
        fmt->p = NULL;
    }
}

static void remove_hdr_lines(bcf_hdr_t *hdr, int type)
{
    int i = 0, nrm = 0;
    while ( i<hdr->nhrec )
    {
        if ( hdr->hrec[i]->type!=type ) { i++; continue; }
        bcf_hrec_t *hrec = hdr->hrec[i];
        if ( type==BCF_HL_FMT )
        {
            // everything except FORMAT/GT
            int id = bcf_hrec_find_key(hrec, "ID");
            if ( id>=0 && !strcmp(hrec->vals[id],"GT") ) { i++; continue; }
        }
        nrm++;
        hdr->nhrec--;
        if ( i < hdr->nhrec )
            memmove(&hdr->hrec[i],&hdr->hrec[i+1],(hdr->nhrec-i)*sizeof(bcf_hrec_t*));
        bcf_hrec_destroy(hrec);
    }
    if ( nrm ) bcf_hdr_sync(hdr);
}

static void init_remove_annots(hand_t *hand)
{
    int keep_info = 0, keep_fmt = 0, keep_flt = 0;
    void *keep = khash_str2int_init();
    kstring_t str = {0,0,0};
    char *ss = hand->remove_annots;
    while ( *ss )
    {
        hand->nrm++;
        hand->rm = (rm_tag_t*) realloc(hand->rm,sizeof(rm_tag_t)*hand->nrm);
        rm_tag_t *tag = &hand->rm[hand->nrm-1];
        tag->key = NULL;

        int type = BCF_HL_GEN;
        if ( !strncasecmp("INFO/",ss,5) ) { type = BCF_HL_INFO; ss += 5; }
        else if ( !strncasecmp("INF/",ss,4) ) { type = BCF_HL_INFO; ss += 4; }
        else if ( !strncasecmp("FORMAT/",ss,7) ) { type = BCF_HL_FMT; ss += 7; }
        else if ( !strncasecmp("FMT/",ss,4) ) { type = BCF_HL_FMT; ss += 4; }
        else if ( !strncasecmp("FILTER/",ss,7) ) { type = BCF_HL_FLT; ss += 7; }
        else if ( !strncasecmp("^INFO/",ss,6) ) { type = BCF_HL_INFO; ss += 6; keep_info = 1; }
        else if ( !strncasecmp("^INF/",ss,5) ) { type = BCF_HL_INFO; ss += 5; keep_info = 1; }
        else if ( !strncasecmp("^FORMAT/",ss,8) ) { type = BCF_HL_FMT; ss += 8; keep_fmt = 1; }
        else if ( !strncasecmp("^FMT/",ss,5) ) { type = BCF_HL_FMT; ss += 5; keep_fmt = 1; }
        else if ( !strncasecmp("^FILTER/",ss,8) ) { type = BCF_HL_FLT; ss += 8; keep_flt = 1; }

        char *se = ss;
        while ( *se && *se!=',' ) se++;
        str.l = 0;
        kputsn(ss, se-ss, &str);

        if ( type==BCF_HL_FLT )
        {
            if ( !keep_flt )
            {
                hand->flt_keep_pass = 1;
                tag->handler = remove_filter;
                tag->key = strdup(str.s);
                tag->hdr_id = bcf_hdr_id2int(hand->hdr, BCF_DT_ID, tag->key);
                if ( !bcf_hdr_idinfo_exists(hand->hdr,BCF_HL_FLT,tag->hdr_id) ) error("Cannot remove %s, not defined in the header.\n", str.s);
                bcf_hdr_remove(hand->hdr_out,BCF_HL_FLT,tag->key);
            }
            else
            {
                int value, ret = khash_str2int_get(keep, str.s, &value);
                if ( ret==-1 ) khash_str2int_set(keep, strdup(str.s),1<<BCF_HL_FLT);
                else khash_str2int_set(keep, str.s, value | 1<<BCF_HL_FLT);
                hand->nrm--;
            }
        }
        else if ( type!=BCF_HL_GEN )
        {
            int id = bcf_hdr_id2int(hand->hdr,BCF_DT_ID,str.s);
            if ( !bcf_hdr_idinfo_exists(hand->hdr,type,id) )
            {
                fprintf(stderr,"Warning: The tag \"%s\" not defined in the header\n", str.s);
                hand->nrm--;
            }
            else if ( (type==BCF_HL_FMT && keep_fmt) || (type==BCF_HL_INFO && keep_info) )
            {
                int value, ret = khash_str2int_get(keep, str.s, &value);
                if ( ret==-1 ) khash_str2int_set(keep, strdup(str.s),1<<type);
                else khash_str2int_set(keep, str.s, value | 1<<type);
                hand->nrm--;
            }
            else
            {
                tag->key = strdup(str.s);
                if ( type==BCF_HL_INFO ) tag->handler = remove_info_tag;
                else if ( type==BCF_HL_FMT ) tag->handler = remove_format_tag;
                bcf_hdr_remove(hand->hdr_out,type,tag->key);
            }
        }
        else if ( !strcasecmp("ID",str.s) ) tag->handler = remove_id;
        else if ( !strcasecmp("FILTER",str.s) )
        {
            tag->handler = remove_filter;
            remove_hdr_lines(hand->hdr_out,BCF_HL_FLT);
        }
        else if ( !strcasecmp("QUAL",str.s) ) tag->handler = remove_qual;
        else if ( !strcasecmp("INFO",str.s) ) 
        {
            tag->handler = remove_info;
            remove_hdr_lines(hand->hdr_out,BCF_HL_INFO);
        }
        else if ( !strcasecmp("FMT",str.s) || !strcasecmp("FORMAT",str.s) )
        {
            tag->handler = remove_format;
            remove_hdr_lines(hand->hdr_out,BCF_HL_FMT);
        }
        else if ( str.l )
        {
            if ( str.s[0]=='#' && str.s[1]=='#' )
                bcf_hdr_remove(hand->hdr_out,BCF_HL_GEN,str.s+2);
            else
                bcf_hdr_remove(hand->hdr_out,BCF_HL_STR,str.s);
            hand->nrm--;
        }

        ss = *se ? se+1 : se;
    }
    free(str.s);
    if ( keep_flt || keep_info || keep_fmt )
    {
        int j;
        for (j=0; j<hand->hdr->nhrec; j++)
        {
            bcf_hrec_t *hrec = hand->hdr->hrec[j];
            if ( hrec->type!=BCF_HL_FLT && hrec->type!=BCF_HL_INFO && hrec->type!=BCF_HL_FMT ) continue;
            if ( !keep_flt && hrec->type==BCF_HL_FLT ) continue;
            if ( !keep_info && hrec->type==BCF_HL_INFO ) continue;
            if ( !keep_fmt && hrec->type==BCF_HL_FMT ) continue;
            int k = bcf_hrec_find_key(hrec,"ID");
            assert( k>=0 ); // this should always be true for valid VCFs
            int value, ret = khash_str2int_get(keep,hrec->vals[k],&value);
            if ( ret==0 && value>>hrec->type ) // keep
            {
                if ( hrec->type==BCF_HL_FLT && !strcmp("PASS",hrec->vals[k]) ) hand->flt_keep_pass = 1;
                continue;
            }
            hand->nrm++;
            hand->rm = (rm_tag_t*) realloc(hand->rm,sizeof(rm_tag_t)*hand->nrm);
            rm_tag_t *tag = &hand->rm[hand->nrm-1];
            if ( hrec->type==BCF_HL_INFO ) tag->handler = remove_info_tag;
            else if ( hrec->type==BCF_HL_FMT ) tag->handler = remove_format_tag;
            else 
            {
                tag->handler = remove_filter;
                tag->hdr_id = bcf_hdr_id2int(hand->hdr, BCF_DT_ID, hrec->vals[k]);
            }
            tag->key = strdup(hrec->vals[k]);
            bcf_hdr_remove(hand->hdr_out,hrec->type,tag->key);
        }
    }
    khash_str2int_destroy_free(keep);
    if ( !hand->nrm ) error("No matching tag in -x %s\n", hand->remove_annots);
    bcf_hdr_sync(hand->hdr_out);
}
*/

/* static void init_header_lines(hand_t *hand) */
/* { */
/*     htsFile *file = hts_open(hand->header_fname, "rb"); */
/*     if ( !file ) error("Error reading %s\n", hand->header_fname); */
/*     kstring_t str = {0,0,0}; */
/*     while ( hts_getline(file, KS_SEP_LINE, &str) > 0 ) */
/*     { */
/*         if ( bcf_hdr_append(hand->hdr_out,str.s) ) error("Could not parse %s: %s\n", hand->header_fname, str.s); */
/*         bcf_hdr_append(hand->hdr,str.s);    // the input file may not have the header line if run with -h (and nothing else) */
/*     } */
/*     hts_close(file); */
/*     free(str.s); */
/*     bcf_hdr_sync(hand->hdr_out); */
/*     bcf_hdr_sync(hand->hdr); */
/* } */

int setter_filter(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    // note: so far this works only with one filter, not a list of filters
    annot_line_t *tab = (annot_line_t*) data;
    if ( tab->cols[col->icol] && tab->cols[col->icol][0]=='.' && !tab->cols[col->icol][1] ) return 0; // don't replace with "."
    hts_expand(int,1,hand->mtmpi,hand->tmpi);
    hand->tmpi[0] = bcf_hdr_id2int(hand->hdr_out, BCF_DT_ID, tab->cols[col->icol]);
    if ( hand->tmpi[0]<0 ) error("The FILTER is not defined in the header: %s\n", tab->cols[col->icol]);
    if ( col->replace==SET_OR_APPEND ) { bcf_add_filter(hand->hdr_out,line,hand->tmpi[0]); return 0; }
    if ( col->replace!=REPLACE_MISSING )
    {
        bcf_update_filter(hand->hdr_out,line,NULL,0);
        bcf_update_filter(hand->hdr_out,line,hand->tmpi,1); 
        return 0; 
    }
    
    // only update missing FILTER
    if ( !(line->unpacked & BCF_UN_FLT) ) bcf_unpack(line, BCF_UN_FLT);
    if ( !line->d.n_flt )
        bcf_update_filter(hand->hdr_out,line,hand->tmpi,1);
    return 0;
}
int vcf_setter_filter(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    int i;
    bcf1_t *rec = (bcf1_t*) data;
    if ( !(rec->unpacked & BCF_UN_FLT) ) bcf_unpack(rec, BCF_UN_FLT);
    if ( !(line->unpacked & BCF_UN_FLT) ) bcf_unpack(line, BCF_UN_FLT);
    if ( !rec->d.n_flt ) return 0;  // don't overwrite with a missing value
    if ( col->replace==SET_OR_APPEND || col->replace==REPLACE_MISSING )
    {
        if ( col->replace==REPLACE_MISSING && line->d.n_flt ) return 0; // only update missing FILTER
        for (i=0; i<rec->d.n_flt; i++)
        {
            const char *flt = bcf_hdr_int2id(hand->files->readers[1].header, BCF_DT_ID, rec->d.flt[i]);
            bcf_add_filter(hand->hdr_out,line,bcf_hdr_id2int(hand->hdr_out, BCF_DT_ID, flt));
        }
        return 0;
    }
    hts_expand(int,rec->d.n_flt,hand->mtmpi,hand->tmpi);
    for (i=0; i<rec->d.n_flt; i++)
    {
        const char *flt = bcf_hdr_int2id(hand->files->readers[1].header, BCF_DT_ID, rec->d.flt[i]);
        hand->tmpi[i] = bcf_hdr_id2int(hand->hdr_out, BCF_DT_ID, flt);
    }
    bcf_update_filter(hand->hdr_out,line,NULL,0);
    bcf_update_filter(hand->hdr_out,line,hand->tmpi,rec->d.n_flt);
    return 0;
}
int setter_id(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    // possible cases:
    //      IN  ANNOT   OUT     ACHIEVED_BY
    //      x   y       x        -c +ID
    //      x   y       y        -c ID
    //      x   y       x,y      -c =ID
    //      x   .       x        -c +ID, ID
    //      x   .       .        -x ID
    //      .   y       y        -c +ID, -c ID
    //
    annot_line_t *tab = (annot_line_t*) data;
    if ( tab->cols[col->icol] && tab->cols[col->icol][0]=='.' && !tab->cols[col->icol][1] ) return 0; // don't replace with "."
    if ( col->replace==SET_OR_APPEND ) return bcf_add_id(hand->hdr_out,line,tab->cols[col->icol]);
    if ( col->replace!=REPLACE_MISSING ) return bcf_update_id(hand->hdr_out,line,tab->cols[col->icol]);

    // running with +ID, only update missing ids
    if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) )
        return bcf_update_id(hand->hdr_out,line,tab->cols[col->icol]);
    return 0;
}
int vcf_setter_id(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    if ( rec->d.id && rec->d.id[0]=='.' && !rec->d.id[1] ) return 0;    // don't replace with "."
    if ( col->replace==SET_OR_APPEND ) return bcf_add_id(hand->hdr_out,line,rec->d.id);
    if ( col->replace!=REPLACE_MISSING ) return bcf_update_id(hand->hdr_out,line,rec->d.id);

    // running with +ID, only update missing ids
    if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) )
        return bcf_update_id(hand->hdr_out,line,rec->d.id);
    return 0;
}
int setter_qual(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol];
    if ( str[0]=='.' && str[1]==0 ) return 0;   // empty

    if ( col->replace==REPLACE_MISSING && !bcf_float_is_missing(line->qual) ) return 0;

    line->qual = strtod(str, &str);
    if ( str == tab->cols[col->icol] )
        error("Could not parse %s at %s:%d .. [%s]\n", col->hdr_key,bcf_seqname(hand->hdr,line),line->pos+1,tab->cols[col->icol]);
    return 0;
}
int vcf_setter_qual(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    if ( bcf_float_is_missing(rec->qual) ) return 0;
    if ( col->replace==REPLACE_MISSING && !bcf_float_is_missing(line->qual) ) return 0;
    line->qual = rec->qual;
    return 0;
}
int setter_info_flag(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol];
    if ( str[0]=='.' && str[1]==0 ) return 0;

    if ( str[0]=='1' && str[1]==0 ) return bcf_update_info_flag(hand->hdr_out,line,col->hdr_key,NULL,1);
    if ( str[0]=='0' && str[1]==0 ) return bcf_update_info_flag(hand->hdr_out,line,col->hdr_key,NULL,0);
    error("Could not parse %s at %s:%d .. [%s]\n", bcf_seqname(hand->hdr,line),line->pos+1,tab->cols[col->icol]);
    return -1;
}
int vcf_setter_info_flag(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int flag = bcf_get_info_flag(hand->files->readers[1].header,rec,col->hdr_key,NULL,NULL);
    bcf_update_info_flag(hand->hdr_out,line,col->hdr_key,NULL,flag);
    return 0;
}
int setter_ARinfo_int32(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, int nals, char **als, int ntmpi)
{
    if ( col->number==BCF_VL_A && ntmpi!=nals-1 && (ntmpi!=1 || hand->tmpi[0]!=bcf_int32_missing || hand->tmpi[1]!=bcf_int32_vector_end) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpi,col->hdr_key,bcf_seqname(hand->hdr,line),line->pos+1);
    else if ( col->number==BCF_VL_R && ntmpi!=nals && (ntmpi!=1 || hand->tmpi[0]!=bcf_int32_missing || hand->tmpi[1]!=bcf_int32_vector_end) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpi,col->hdr_key,bcf_seqname(hand->hdr,line),line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    int *map = vcmp_map_ARvalues(hand->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%d\n");

    // fill in any missing values in the target VCF (or all, if not present)
    int ntmpi2 = bcf_get_info_float(hand->hdr, line, col->hdr_key, &hand->tmpi2, &hand->mtmpi2);
    if ( ntmpi2 < ndst ) hts_expand(int32_t,ndst,hand->mtmpi2,hand->tmpi2);

    int i;
    for (i=0; i<ndst; i++)
    {
        if ( map[i]<0 )
        {
            if ( ntmpi2 < ndst ) hand->tmpi2[i] = bcf_int32_missing;
            continue;
        }
        if ( ntmpi2==ndst && col->replace==REPLACE_MISSING
                && hand->tmpi2[i]!=bcf_int32_missing
                && hand->tmpi2[i]!=bcf_int32_vector_end ) continue;

        hand->tmpi2[i] = hand->tmpi[ map[i] ];
    }
    bcf_update_info_int32(hand->hdr_out,line,col->hdr_key,hand->tmpi2,ndst);
    return 0;
}
int setter_info_int(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol], *end = str;
    if ( str[0]=='.' && str[1]==0 ) return 0;

    int ntmpi = 0;
    while ( *end )
    {
        int val = strtol(str, &end, 10); 
        if ( end==str )
            error("Could not parse %s at %s:%d .. [%s]\n", bcf_seqname(hand->hdr,line),line->pos+1,tab->cols[col->icol]);
        ntmpi++;
        hts_expand(int32_t,ntmpi,hand->mtmpi,hand->tmpi);
        hand->tmpi[ntmpi-1] = val;
        str = end+1;
    }

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_int32(hand,line,col,tab->nals,tab->als,ntmpi);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_int32(hand->hdr, line, col->hdr_key, &hand->tmpi2, &hand->mtmpi2);
        if ( ret>0 && hand->tmpi2[0]!=bcf_int32_missing ) return 0;
    }

    bcf_update_info_int32(hand->hdr_out,line,col->hdr_key,hand->tmpi,ntmpi);
    return 0;
}
int vcf_setter_info_int(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int ntmpi = bcf_get_info_int32(hand->files->readers[1].header,rec,col->hdr_key,&hand->tmpi,&hand->mtmpi);
    if ( ntmpi < 0 ) return 0;    // nothing to add

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_int32(hand,line,col,rec->n_allele,rec->d.allele,ntmpi);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_int32(hand->hdr, line, col->hdr_key, &hand->tmpi2, &hand->mtmpi2);
        if ( ret>0 && hand->tmpi2[0]!=bcf_int32_missing ) return 0;
    }

    bcf_update_info_int32(hand->hdr_out,line,col->hdr_key,hand->tmpi,ntmpi);
    return 0;
}
int setter_ARinfo_real(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, int nals, char **als, int ntmpf)
{
    if ( col->number==BCF_VL_A && ntmpf!=nals-1 && (ntmpf!=1 || !bcf_float_is_missing(hand->tmpf[0]) || !bcf_float_is_vector_end(hand->tmpf[0])) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpf,col->hdr_key,bcf_seqname(hand->hdr,line),line->pos+1);
    else if ( col->number==BCF_VL_R && ntmpf!=nals && (ntmpf!=1 || !bcf_float_is_missing(hand->tmpf[0]) || !bcf_float_is_vector_end(hand->tmpf[0])) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpf,col->hdr_key,bcf_seqname(hand->hdr,line),line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    int *map = vcmp_map_ARvalues(hand->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%d\n");

    // fill in any missing values in the target VCF (or all, if not present)
    int ntmpf2 = bcf_get_info_float(hand->hdr, line, col->hdr_key, &hand->tmpf2, &hand->mtmpf2);
    if ( ntmpf2 < ndst ) hts_expand(float,ndst,hand->mtmpf2,hand->tmpf2);

    int i;
    for (i=0; i<ndst; i++)
    {
        if ( map[i]<0 )
        {
            if ( ntmpf2 < ndst ) bcf_float_set_missing(hand->tmpf2[i]);
            continue;
        }
        if ( ntmpf2==ndst && col->replace==REPLACE_MISSING
                && !bcf_float_is_missing(hand->tmpf2[i])
                && !bcf_float_is_vector_end(hand->tmpf2[i]) ) continue;

        hand->tmpf2[i] = hand->tmpf[ map[i] ];
    }
    bcf_update_info_float(hand->hdr_out,line,col->hdr_key,hand->tmpf2,ndst);
    return 0;
}
int setter_info_real(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol], *end = str;
    if ( str[0]=='.' && str[1]==0 ) return 0;

    int ntmpf = 0;
    while ( *end )
    {
        double val = strtod(str, &end);
        if ( end==str )
            error("Could not parse %s at %s:%d .. [%s]\n", bcf_seqname(hand->hdr,line),line->pos+1,tab->cols[col->icol]);
        ntmpf++;
        hts_expand(float,ntmpf,hand->mtmpf,hand->tmpf);
        hand->tmpf[ntmpf-1] = val;
        str = end+1;
    }

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_real(hand,line,col,tab->nals,tab->als,ntmpf);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_float(hand->hdr, line, col->hdr_key, &hand->tmpf2, &hand->mtmpf2);
        if ( ret>0 && !bcf_float_is_missing(hand->tmpf2[0]) ) return 0;
    }

    bcf_update_info_float(hand->hdr_out,line,col->hdr_key,hand->tmpf,ntmpf);
    return 0;
}
int vcf_setter_info_real(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int ntmpf = bcf_get_info_float(hand->files->readers[1].header,rec,col->hdr_key,&hand->tmpf,&hand->mtmpf);
    if ( ntmpf < 0 ) return 0;    // nothing to add

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_real(hand,line,col,rec->n_allele,rec->d.allele,ntmpf);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_float(hand->hdr, line, col->hdr_key, &hand->tmpf2, &hand->mtmpf2);
        if ( ret>0 && !bcf_float_is_missing(hand->tmpf2[0]) ) return 0;
    }

    bcf_update_info_float(hand->hdr_out,line,col->hdr_key,hand->tmpf,ntmpf);
    return 0;
}
int copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst); // see vcfmerge.c
int setter_ARinfo_string(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, int nals, char **als)
{
    int nsrc = 1, lsrc = 0;
    while ( hand->tmps[lsrc] )
    {
        if ( hand->tmps[lsrc]==',' ) nsrc++;
        lsrc++;
    }
    if ( col->number==BCF_VL_A && nsrc!=nals-1 && (nsrc!=1 || hand->tmps[0]!='.' || hand->tmps[1]!=0 ) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", nsrc,col->hdr_key,bcf_seqname(hand->hdr,line),line->pos+1);
    else if ( col->number==BCF_VL_R && nsrc!=nals && (nsrc!=1 || hand->tmps[0]!='.' || hand->tmps[1]!=0 ) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", nsrc,col->hdr_key,bcf_seqname(hand->hdr,line),line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    int *map = vcmp_map_ARvalues(hand->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%d\n");

    // fill in any missing values in the target VCF (or all, if not present)
    int i, empty = 0, nstr, mstr = hand->tmpks.m;
    nstr = bcf_get_info_string(hand->hdr, line, col->hdr_key, &hand->tmpks.s, &mstr); 
    hand->tmpks.m = mstr;
    if ( nstr<0 || (nstr==1 && hand->tmpks.s[0]=='.' && hand->tmpks.s[1]==0) )
    {
        empty = 0;
        hand->tmpks.l = 0;
        kputc('.',&hand->tmpks);
        for (i=1; i<ndst; i++) kputs(",.",&hand->tmpks);
    }
    else hand->tmpks.l = nstr;
    for (i=0; i<ndst; i++)
    {
        if ( map[i]<0 )
        {
            if ( empty ) copy_string_field(".",0,1,&hand->tmpks,i);
            continue;
        }
        if ( col->replace==REPLACE_MISSING )
        {
            // Do not replace filled values. The field must be looked up again because
            // of realloc in copy_string_field
            int n = 0;
            char *str = hand->tmpks.s;
            while ( *str && n<i )
            {
                if ( *str==',' ) n++;
                str++;
            }
            if ( str[0]!='.' || (str[1]!=',' && str[1]!=0) ) continue;  // value already set
        }
        int ret = copy_string_field(hand->tmps,map[i],lsrc,&hand->tmpks,i);
        assert( ret==0 );
    }
    bcf_update_info_string(hand->hdr_out,line,col->hdr_key,hand->tmpks.s);
    return 0;
}
int setter_info_str(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    int len = strlen(tab->cols[col->icol]);
    if ( !len ) return 0;
    hts_expand(char,len+1,hand->mtmps,hand->tmps);
    memcpy(hand->tmps,tab->cols[col->icol],len+1);
    if ( hand->tmps[0]=='.' && hand->tmps[1]==0 ) return 0;

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_string(hand,line,col,tab->nals,tab->als);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_string(hand->hdr, line, col->hdr_key, &hand->tmps2, &hand->mtmps2);
        if ( ret>0 && (hand->tmps2[0]!='.' || hand->tmps2[1]!=0) ) return 0;
    }

    bcf_update_info_string(hand->hdr_out,line,col->hdr_key,hand->tmps);
    return 0;
}
int vcf_setter_info_str(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int ntmps = bcf_get_info_string(hand->files->readers[1].header,rec,col->hdr_key,&hand->tmps,&hand->mtmps);
    if ( ntmps < 0 ) return 0;    // nothing to add

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_string(hand,line,col,rec->n_allele,rec->d.allele);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_string(hand->hdr, line, col->hdr_key, &hand->tmps2, &hand->mtmps2);
        if ( ret>0 && (hand->tmps2[0]!='.' || hand->tmps2[1]!=0) ) return 0;
    }

    bcf_update_info_string(hand->hdr_out,line,col->hdr_key,hand->tmps);
    return 0;
}
int vcf_setter_format_gt(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int nsrc = bcf_get_genotypes(hand->files->readers[1].header,rec,&hand->tmpi,&hand->mtmpi);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    if ( !hand->sample_map )
        return bcf_update_genotypes(hand->hdr_out,line,hand->tmpi,nsrc);

    int i, j, ndst = bcf_get_genotypes(hand->hdr,line,&hand->tmpi2,&hand->mtmpi2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(hand->hdr_out);
    nsrc /= bcf_hdr_nsamples(hand->files->readers[1].header);
    if ( ndst<=0 )  // field not present in dst file
    {
        if ( col->replace==REPLACE_EXISTING ) return 0;
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(hand->hdr_out), hand->mtmpi2, hand->tmpi2);
        for (i=0; i<bcf_hdr_nsamples(hand->hdr_out); i++)
        {
            int32_t *dst = hand->tmpi2 + nsrc*i;
            if ( hand->sample_map[i]==-1 )
            {
                dst[0] = bcf_gt_missing;
                for (j=1; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = hand->tmpi + nsrc*hand->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_genotypes(hand->hdr_out,line,hand->tmpi2,nsrc*bcf_hdr_nsamples(hand->hdr_out));
    }
    else if ( ndst >= nsrc )     
    {
        for (i=0; i<bcf_hdr_nsamples(hand->hdr_out); i++)
        {
            if ( hand->sample_map[i]==-1 ) continue;
            int32_t *src = hand->tmpi  + nsrc*hand->sample_map[i];
            int32_t *dst = hand->tmpi2 + ndst*i;
            if ( col->replace==REPLACE_EXISTING && bcf_gt_is_missing(dst[0]) ) continue;
            if ( col->replace==REPLACE_MISSING  && !bcf_gt_is_missing(dst[0]) ) continue;
            for (j=0; j<nsrc; j++) dst[j] = src[j];
            for (; j<ndst; j++) dst[j] = bcf_int32_vector_end;
        }
        return bcf_update_genotypes(hand->hdr_out,line,hand->tmpi2,ndst*bcf_hdr_nsamples(hand->hdr_out));
    }
    else    // ndst < nsrc
    {
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(hand->hdr_out), hand->mtmpi3, hand->tmpi3);
        for (i=0; i<bcf_hdr_nsamples(hand->hdr_out); i++)
        {
            int32_t *ori = hand->tmpi2 + ndst*i;
            int32_t *dst = hand->tmpi3 + nsrc*i;
            int keep_ori = 0;
            if ( hand->sample_map[i]==-1 ) keep_ori = 1;
            else if ( col->replace==REPLACE_EXISTING && bcf_gt_is_missing(ori[0]) ) keep_ori = 1;
            else if ( col->replace==REPLACE_MISSING  && !bcf_gt_is_missing(ori[0]) ) keep_ori = 1;
            if ( keep_ori )
            {
                for (j=0; j<ndst; j++) dst[j] = ori[j];
                for (; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = hand->tmpi + nsrc*hand->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_genotypes(hand->hdr_out,line,hand->tmpi3,nsrc*bcf_hdr_nsamples(hand->hdr_out));
    }
}
int count_vals(annot_line_t *tab, int icol_beg, int icol_end)
{
    int i, nmax = 0;
    for (i=icol_beg; i<icol_end; i++)
    {
        char *str = tab->cols[i], *end = str;
        if ( str[0]=='.' && !str[1] ) 
        {
            // missing value
            if ( !nmax ) nmax = 1;
            continue;
        }
        int n = 1;
        while ( *end )
        {
            if ( *end==',' ) n++;
            end++;
        }
        if ( nmax<n ) nmax = n;
    }
    return nmax;
}
int setter_format_int(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    int nsmpl = bcf_hdr_nsamples(hand->hdr_out);
    assert( col->icol+nsmpl <= tab->ncols );
    int nvals = count_vals(tab,col->icol,col->icol+nsmpl);
    assert( nvals>0 );
    hts_expand(int32_t,nvals*nsmpl,hand->mtmpi,hand->tmpi);

    int icol = col->icol, ismpl;
    for (ismpl=0; ismpl<nsmpl; ismpl++)
    {
        int32_t *ptr = hand->tmpi + ismpl*nvals;
        int ival = 0;

        char *str = tab->cols[icol];
        while ( *str )
        {
            if ( str[0]=='.' && (!str[1] || str[1]==',') )  // missing value
            {
                ptr[ival++] = bcf_int32_missing;
                str += str[1] ? 2 : 1;
                continue;
            }

            char *end = str;
            ptr[ival] = strtol(str, &end, 10); 
            if ( end==str )
                error("Could not parse %s at %s:%d .. [%s]\n", col->hdr_key,bcf_seqname(hand->hdr,line),line->pos+1,tab->cols[col->icol]);

            ival++;
            str = *end ? end+1 : end;
        }
        while ( ival<nvals ) ptr[ival++] = bcf_int32_vector_end;
        icol++;
    }
    return bcf_update_format_int32(hand->hdr_out,line,col->hdr_key,hand->tmpi,nsmpl*nvals);
}
int setter_format_real(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    int nsmpl = bcf_hdr_nsamples(hand->hdr_out);
    assert( col->icol+nsmpl <= tab->ncols );
    int nvals = count_vals(tab,col->icol,col->icol+nsmpl);
    assert( nvals>0 );
    hts_expand(float,nvals*nsmpl,hand->mtmpf,hand->tmpf);

    int icol = col->icol, ismpl;
    for (ismpl=0; ismpl<nsmpl; ismpl++)
    {
        float *ptr = hand->tmpf + ismpl*nvals;
        int ival = 0;

        char *str = tab->cols[icol];
        while ( *str )
        {
            if ( str[0]=='.' && (!str[1] || str[1]==',') )  // missing value
            {
                bcf_float_set_missing(ptr[ival]); 
                ival++;
                str += str[1] ? 2 : 1;
                continue;
            }

            char *end = str;
            ptr[ival] = strtod(str, &end); 
            if ( end==str )
                error("Could not parse %s at %s:%d .. [%s]\n", col->hdr_key,bcf_seqname(hand->hdr,line),line->pos+1,tab->cols[col->icol]);

            ival++;
            str = *end ? end+1 : end;
        }
        while ( ival<nvals ) { bcf_float_set_vector_end(ptr[ival]); ival++; }
        icol++;
    }
    return bcf_update_format_float(hand->hdr_out,line,col->hdr_key,hand->tmpf,nsmpl*nvals);
}
int setter_format_str(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    int nsmpl = bcf_hdr_nsamples(hand->hdr_out);
    assert( col->icol+nsmpl <= tab->ncols );

    int i, max_len = 0;
    for (i=col->icol; i<col->icol+nsmpl; i++)
    {
        int len = strlen(tab->cols[i]);
        if ( max_len < len ) max_len = len;
    }
    hts_expand(char,max_len*nsmpl,hand->mtmps,hand->tmps);

    int icol = col->icol, ismpl;
    for (ismpl=0; ismpl<nsmpl; ismpl++)
    {
        char *ptr = hand->tmps + ismpl*max_len;
        char *str = tab->cols[icol];
        i = 0;
        while ( str[i] )
        {
            ptr[i] = str[i];
            i++;
        }
        while ( i<max_len ) ptr[i++] = 0;
        icol++;
    }
    return bcf_update_format_char(hand->hdr_out,line,col->hdr_key,hand->tmps,nsmpl*max_len);
}
int vcf_setter_format_int(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int nsrc = bcf_get_format_int32(hand->files->readers[1].header,rec,col->hdr_key,&hand->tmpi,&hand->mtmpi);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    if ( !hand->sample_map )
        return bcf_update_format_int32(hand->hdr_out,line,col->hdr_key,hand->tmpi,nsrc);

    int i, j, ndst = bcf_get_format_int32(hand->hdr,line,col->hdr_key,&hand->tmpi2,&hand->mtmpi2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(hand->hdr_out);
    nsrc /= bcf_hdr_nsamples(hand->files->readers[1].header);
    if ( ndst<=0 )
    {
        if ( col->replace==REPLACE_EXISTING ) return 0;    // overwrite only if present
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(hand->hdr_out), hand->mtmpi2, hand->tmpi2);
        for (i=0; i<bcf_hdr_nsamples(hand->hdr_out); i++)
        {
            int32_t *dst = hand->tmpi2 + nsrc*i;
            if ( hand->sample_map[i]==-1 )
            {
                dst[0] = bcf_int32_missing;
                for (j=1; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = hand->tmpi + nsrc*hand->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_int32(hand->hdr_out,line,col->hdr_key,hand->tmpi2,nsrc*bcf_hdr_nsamples(hand->hdr_out));
    }
    else if ( ndst >= nsrc )     
    {
        for (i=0; i<bcf_hdr_nsamples(hand->hdr_out); i++)
        {
            if ( hand->sample_map[i]==-1 ) continue;
            int32_t *src = hand->tmpi  + nsrc*hand->sample_map[i];
            int32_t *dst = hand->tmpi2 + ndst*i;
            if ( col->replace==REPLACE_EXISTING && dst[0]==bcf_int32_missing ) continue;
            if ( col->replace==REPLACE_MISSING  && dst[0]!=bcf_int32_missing ) continue;
            for (j=0; j<nsrc; j++) dst[j] = src[j];
            for (; j<ndst; j++) dst[j] = bcf_int32_vector_end;
        }
        return bcf_update_format_int32(hand->hdr_out,line,col->hdr_key,hand->tmpi2,ndst*bcf_hdr_nsamples(hand->hdr_out));
    }
    else    // ndst < nsrc
    {
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(hand->hdr_out), hand->mtmpi3, hand->tmpi3);
        for (i=0; i<bcf_hdr_nsamples(hand->hdr_out); i++)
        {
            int32_t *ori = hand->tmpi2 + ndst*i;
            int32_t *dst = hand->tmpi3 + nsrc*i;
            int keep_ori = 0;
            if ( hand->sample_map[i]==-1 ) keep_ori = 1;
            else if ( col->replace==REPLACE_EXISTING && ori[0]==bcf_int32_missing ) keep_ori = 1;
            else if ( col->replace==REPLACE_MISSING  && ori[0]!=bcf_int32_missing ) keep_ori = 1;
            if ( keep_ori )
            {
                for (j=0; j<ndst; j++) dst[j] = ori[j];
                for (; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = hand->tmpi + nsrc*hand->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_int32(hand->hdr_out,line,col->hdr_key,hand->tmpi3,nsrc*bcf_hdr_nsamples(hand->hdr_out));
    }
}
int vcf_setter_format_real(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int nsrc = bcf_get_format_float(hand->files->readers[1].header,rec,col->hdr_key,&hand->tmpf,&hand->mtmpf);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    if ( !hand->sample_map )
        return bcf_update_format_float(hand->hdr_out,line,col->hdr_key,hand->tmpf,nsrc);

    int i, j, ndst = bcf_get_format_float(hand->hdr,line,col->hdr_key,&hand->tmpf2,&hand->mtmpf2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(hand->hdr_out);
    nsrc /= bcf_hdr_nsamples(hand->files->readers[1].header);
    if ( ndst<=0 )
    {
        if ( col->replace==REPLACE_EXISTING ) return 0;    // overwrite only if present
        hts_expand(float, nsrc*bcf_hdr_nsamples(hand->hdr_out), hand->mtmpf2, hand->tmpf2);
        for (i=0; i<bcf_hdr_nsamples(hand->hdr_out); i++)
        {
            float *dst = hand->tmpf2 + nsrc*i;
            if ( hand->sample_map[i]==-1 )
            {
                bcf_float_set_missing(dst[0]);
                for (j=1; j<nsrc; j++) bcf_float_set_vector_end(dst[j]);
            }
            else
            {
                float *src = hand->tmpf + nsrc*hand->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_float(hand->hdr_out,line,col->hdr_key,hand->tmpf2,nsrc*bcf_hdr_nsamples(hand->hdr_out));
    }
    else if ( ndst >= nsrc )     
    {
        for (i=0; i<bcf_hdr_nsamples(hand->hdr_out); i++)
        {
            if ( hand->sample_map[i]==-1 ) continue;
            float *src = hand->tmpf  + nsrc*hand->sample_map[i];
            float *dst = hand->tmpf2 + ndst*i;
            if ( col->replace==REPLACE_EXISTING && bcf_float_is_missing(dst[0]) ) continue;
            if ( col->replace==REPLACE_MISSING  && !bcf_float_is_missing(dst[0]) ) continue;
            for (j=0; j<nsrc; j++) dst[j] = src[j];
            for (; j<ndst; j++) bcf_float_set_vector_end(dst[j]);
        }
        return bcf_update_format_float(hand->hdr_out,line,col->hdr_key,hand->tmpf2,ndst*bcf_hdr_nsamples(hand->hdr_out));
    }
    else    // ndst < nsrc
    {
        hts_expand(float, nsrc*bcf_hdr_nsamples(hand->hdr_out), hand->mtmpf3, hand->tmpf3);
        for (i=0; i<bcf_hdr_nsamples(hand->hdr_out); i++)
        {
            float *ori = hand->tmpf2 + ndst*i;
            float *dst = hand->tmpf3 + nsrc*i;
            int keep_ori = 0;
            if ( hand->sample_map[i]==-1 ) keep_ori = 1;
            else if ( col->replace==REPLACE_EXISTING && bcf_float_is_missing(ori[0]) ) keep_ori = 1;
            else if ( col->replace==REPLACE_MISSING  && !bcf_float_is_missing(ori[0]) ) keep_ori = 1;
            if ( keep_ori )
            {
                for (j=0; j<ndst; j++) dst[j] = ori[j];
                for (; j<nsrc; j++) bcf_float_set_vector_end(dst[j]);
            }
            else
            {
                float *src = hand->tmpf + nsrc*hand->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_float(hand->hdr_out,line,col->hdr_key,hand->tmpf3,nsrc*bcf_hdr_nsamples(hand->hdr_out));
    }
}
int vcf_setter_format_str(struct anno_handler *hand, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    hand->tmpp[0] = hand->tmps;
    int ret = bcf_get_format_string(hand->files->readers[1].header,rec,col->hdr_key,&hand->tmpp,&hand->mtmps);
    hand->tmps = hand->tmpp[0]; // tmps might be realloced
    if ( ret==-3 ) return 0;    // the tag is not present
    if ( ret<=0 ) return 1;     // error

    if ( !hand->sample_map )
        return bcf_update_format_string(hand->hdr_out,line,col->hdr_key,(const char**)hand->tmpp,bcf_hdr_nsamples(hand->hdr_out));

    int i;
    hand->tmpp2[0] = hand->tmps2;
    ret = bcf_get_format_string(hand->hdr,line,col->hdr_key,&hand->tmpp2,&hand->mtmps2);
    hand->tmps2 = hand->tmpp2[0];   // tmps2 might be realloced

    if ( ret<=0 )   // not present in dst
    {
        hts_expand(char,bcf_hdr_nsamples(hand->hdr_out)*2,hand->mtmps2,hand->tmps2);
        for (i=0; i<bcf_hdr_nsamples(hand->hdr_out); i++)
        {
            hand->tmps2[2*i]   = '.';
            hand->tmps2[2*i+1] = 0;
            hand->tmpp2[i] = hand->tmps2+2*i;
        }
    }

    for (i=0; i<bcf_hdr_nsamples(hand->hdr_out); i++)
    {
        int isrc = hand->sample_map[i];
        if ( isrc==-1 ) continue;
        hand->tmpp2[i] = hand->tmpp[isrc];
    }
    return bcf_update_format_string(hand->hdr_out,line,col->hdr_key,(const char**)hand->tmpp2,bcf_hdr_nsamples(hand->hdr_out));
}

/* static void set_samples(struct anno_handler *hand, bcf_hdr_t *src, bcf_hdr_t *dst, int need_samples) */
/* { */
/*     int i; */
/*     if ( !hand->sample_names ) */
/*     { */
/*         int nmatch = 0, order_ok = 1; */
/*         for (i=0; i<bcf_hdr_nsamples(src); i++) */
/*         { */
/*             int id = bcf_hdr_id2int(dst, BCF_DT_SAMPLE, src->samples[i]); */
/*             if ( id!=-1 )  */
/*             { */
/*                 nmatch++; */
/*                 if ( i!=id ) order_ok = 0; */
/*             } */
/*         } */
/*         if ( bcf_hdr_nsamples(src)==bcf_hdr_nsamples(dst) && nmatch==bcf_hdr_nsamples(src) && order_ok && !need_samples )  */
/*             return;    // the same samples in both files */

/*         if ( !nmatch ) error("No matching samples found in the source and the destination file\n"); */
/*         if ( nmatch!=bcf_hdr_nsamples(src) || nmatch!=bcf_hdr_nsamples(dst) ) fprintf(stderr,"%d sample(s) in common\n", nmatch); */

/*         hand->nsample_map = bcf_hdr_nsamples(dst); */
/*         hand->sample_map  = (int*) malloc(sizeof(int)*hand->nsample_map); */
/*         for (i=0; i<hand->nsample_map; i++) */
/*         { */
/*             int id = bcf_hdr_id2int(src, BCF_DT_SAMPLE, dst->samples[i]); */
/*             hand->sample_map[i] = id;   // idst -> isrc, -1 if not present */
/*         } */
/*         return; */
/*     } */
/*     hand->nsample_map = bcf_hdr_nsamples(dst); */
/*     hand->sample_map  = (int*) malloc(sizeof(int)*hand->nsample_map); */
/*     for (i=0; i<hand->nsample_map; i++) hand->sample_map[i] = -1; */
/*     int nsamples = 0; */
/*     char **samples = hts_readlist(hand->sample_names, hand->sample_is_file, &nsamples); */
/*     for (i=0; i<nsamples; i++) */
/*     { */
/*         int isrc, idst; */
/*         char *ss = samples[i], *se = samples[i]; */
/*         while ( *se && !isspace(*se) ) se++; */
/*         if ( !*se )  */
/*         { */
/*             // only one sample name */
/*             isrc = bcf_hdr_id2int(src, BCF_DT_SAMPLE,ss); */
/*             if ( isrc==-1 ) error("Sample \"%s\" not found in the source file\n", ss); */
/*             idst = bcf_hdr_id2int(dst, BCF_DT_SAMPLE,ss); */
/*             if ( idst==-1 ) error("Sample \"%s\" not found in the destination file\n", ss); */
/*             hand->sample_map[idst] = isrc; */
/*             continue; */
/*         } */
/*         *se = 0; */
/*         isrc = bcf_hdr_id2int(src, BCF_DT_SAMPLE,ss); */
/*         if ( isrc==-1 ) error("Sample \"%s\" not found in the source file\n", ss); */

/*         ss = se+1; */
/*         while ( isspace(*ss) ) ss++; */
/*         se = ss; */
/*         while ( *se && !isspace(*se) ) se++; */

/*         idst = bcf_hdr_id2int(dst, BCF_DT_SAMPLE,ss); */
/*         if ( idst==-1 ) error("Sample \"%s\" not found in the destination file\n", ss); */

/*         hand->sample_map[idst] = isrc; */
/*     } */
/*     for (i=0; i<nsamples; i++) free(samples[i]); */
/*     free(samples); */
/* } */
/* static char *columns_complement(char *columns, void **skip_info, void **skip_fmt) */
/* { */
/*     kstring_t str = {0,0,0}; */
/*     char *ss = columns, *se = ss; */
/*     while ( *ss ) */
/*     { */
/*         if ( *se && *se!=',' ) { se++; continue; } */
/*         if ( *ss!='^' ) */
/*         { */
/*             if ( str.l ) kputc(',',&str); */
/*             kputsn(ss, se-ss, &str); */
/*             if ( !*se ) break; */
/*             ss = ++se; */
/*             continue; */
/*         } */

/*         if ( !strncasecmp("^INFO/",ss,6) ) */
/*         { */
/*             if ( !*skip_info ) */
/*             { */
/*                 *skip_info = khash_str2int_init(); */
/*                 if ( str.l ) kputc(',',&str); */
/*                 kputs("INFO",&str); */
/*             } */
/*             char tmp = *se; *se = 0; */
/*             khash_str2int_inc(*skip_info, strdup(ss+6)); */
/*             *se = tmp; */
/*         } */
/*         else if ( !strncasecmp("^FORMAT/",ss,8) || !strncasecmp("^FMT/",ss,5) ) */
/*         { */
/*             int n = !strncasecmp("^FMT/",ss,5) ? 5 : 8; */
/*             if ( !*skip_fmt ) */
/*             { */
/*                 *skip_fmt = khash_str2int_init(); */
/*                 if ( str.l ) kputc(',',&str); */
/*                 kputs("FORMAT",&str); */
/*             } */
/*             char tmp = *se; *se = 0; */
/*             khash_str2int_inc(*skip_fmt, strdup(ss+n)); */
/*             *se = tmp; */
/*         } */
/*         else */
/*         { */
/*             if ( !*skip_info ) */
/*             { */
/*                 *skip_info = khash_str2int_init(); */
/*                 if ( str.l ) kputc(',',&str); */
/*                 kputs("INFO",&str); */
/*             } */
/*             char tmp = *se; *se = 0; */
/*             khash_str2int_inc(*skip_info, strdup(ss+1)); */
/*             *se = tmp; */
/*         } */

/*         if ( !*se ) break; */
/*         ss = ++se; */
/*     } */
/*     free(columns); */
/*     return str.s; */
/* } */
/* static void init_columns(struct anno_handler *hand) */
/* { */
/*     void *skip_fmt = NULL, *skip_info = NULL; */
/*     if ( hand->tgts_is_vcf ) */
/*         hand->columns = columns_complement(hand->columns, &skip_info, &skip_fmt); */

/*     kstring_t str = {0,0,0}, tmp = {0,0,0}; */
/*     char *ss = hand->columns, *se = ss; */
/*     hand->ncols = 0; */
/*     int icol = -1, has_fmt_str = 0, force_samples = -1; */
/*     while ( *ss ) */
/*     { */
/*         if ( *se && *se!=',' ) { se++; continue; } */
/*         int replace = REPLACE_ALL; */
/*         if ( *ss=='+' ) { replace = REPLACE_MISSING; ss++; } */
/*         else if ( *ss=='-' ) { replace = REPLACE_EXISTING; ss++; } */
/*         else if ( *ss=='=' ) { replace = SET_OR_APPEND; ss++; } */
/*         icol++; */
/*         str.l = 0; */
/*         kputsn(ss, se-ss, &str); */
/*         if ( !str.s[0] || !strcasecmp("-",str.s) ) ; */
/*         else if ( !strcasecmp("CHROM",str.s) ) hand->chr_idx = icol; */
/*         else if ( !strcasecmp("POS",str.s) ) hand->from_idx = icol; */
/*         else if ( !strcasecmp("FROM",str.s) ) hand->from_idx = icol; */
/*         else if ( !strcasecmp("TO",str.s) ) hand->to_idx = icol; */
/*         else if ( !strcasecmp("REF",str.s) ) hand->ref_idx = icol; */
/*         else if ( !strcasecmp("ALT",str.s) ) hand->alt_idx = icol; */
/*         else if ( !strcasecmp("ID",str.s) ) */
/*         { */
/*             if ( replace==REPLACE_EXISTING ) error("Apologies, the -ID feature has not been implemented yet.\n"); */
/*             hand->ncols++; hand->cols = (annot_col_t*) realloc(hand->cols,sizeof(annot_col_t)*hand->ncols); */
/*             annot_col_t *col = &hand->cols[hand->ncols-1]; */
/*             col->icol = icol; */
/*             col->replace = replace; */
/*             col->setter = hand->tgts_is_vcf ? vcf_setter_id : setter_id; */
/*             col->hdr_key = strdup(str.s); */
/*         } */
/*         else if ( !strcasecmp("FILTER",str.s) ) */
/*         { */
/*             if ( replace==REPLACE_EXISTING ) error("Apologies, the -FILTER feature has not been implemented yet.\n"); */
/*             hand->ncols++; hand->cols = (annot_col_t*) realloc(hand->cols,sizeof(annot_col_t)*hand->ncols); */
/*             annot_col_t *col = &hand->cols[hand->ncols-1]; */
/*             col->icol = icol; */
/*             col->replace = replace; */
/*             col->setter = hand->tgts_is_vcf ? vcf_setter_filter : setter_filter; */
/*             col->hdr_key = strdup(str.s); */
/*             if ( hand->tgts_is_vcf ) */
/*             { */
/*                 bcf_hdr_t *tgts_hdr = hand->files->readers[1].header; */
/*                 int j; */
/*                 for (j=0; j<tgts_hdr->nhrec; j++) */
/*                 { */
/*                     bcf_hrec_t *hrec = tgts_hdr->hrec[j]; */
/*                     if ( hrec->type!=BCF_HL_FLT ) continue; */
/*                     int k = bcf_hrec_find_key(hrec,"ID"); */
/*                     assert( k>=0 ); // this should always be true for valid VCFs */
/*                     tmp.l = 0; */
/*                     bcf_hrec_format(hrec, &tmp); */
/*                     bcf_hdr_append(hand->hdr_out, tmp.s); */
/*                 } */
/*                 bcf_hdr_sync(hand->hdr_out); */
/*             } */
/*         } */
/*         else if ( !strcasecmp("QUAL",str.s) ) */
/*         { */
/*             if ( replace==REPLACE_EXISTING ) error("Apologies, the -QUAL feature has not been implemented yet.\n"); */
/*             if ( replace==SET_OR_APPEND ) error("Apologies, the =QUAL feature has not been implemented yet.\n"); */
/*             hand->ncols++; hand->cols = (annot_col_t*) realloc(hand->cols,sizeof(annot_col_t)*hand->ncols); */
/*             annot_col_t *col = &hand->cols[hand->ncols-1]; */
/*             col->icol = icol; */
/*             col->replace = replace; */
/*             col->setter = hand->tgts_is_vcf ? vcf_setter_qual : setter_qual; */
/*             col->hdr_key = strdup(str.s); */
/*         } */
/*         else if ( hand->tgts_is_vcf && !strcasecmp("INFO",str.s) ) // All INFO fields */
/*         { */
/*             if ( replace==REPLACE_EXISTING ) error("Apologies, the -INFO/TAG feature has not been implemented yet.\n"); */
/*             if ( replace==SET_OR_APPEND ) error("Apologies, the =INFO/TAG feature has not been implemented yet.\n"); */
/*             bcf_hdr_t *tgts_hdr = hand->files->readers[1].header; */
/*             int j; */
/*             for (j=0; j<tgts_hdr->nhrec; j++) */
/*             { */
/*                 bcf_hrec_t *hrec = tgts_hdr->hrec[j]; */
/*                 if ( hrec->type!=BCF_HL_INFO ) continue; */
/*                 int k = bcf_hrec_find_key(hrec,"ID"); */
/*                 assert( k>=0 ); // this should always be true for valid VCFs */
/*                 if ( skip_info && khash_str2int_has_key(skip_info,hrec->vals[k]) ) continue; */
/*                 tmp.l = 0; */
/*                 bcf_hrec_format(hrec, &tmp); */
/*                 bcf_hdr_append(hand->hdr_out, tmp.s); */
/*                 bcf_hdr_sync(hand->hdr_out); */
/*                 int hdr_id = bcf_hdr_id2int(hand->hdr_out, BCF_DT_ID, hrec->vals[k]); */
/*                 hand->ncols++; hand->cols = (annot_col_t*) realloc(hand->cols,sizeof(annot_col_t)*hand->ncols); */
/*                 annot_col_t *col = &hand->cols[hand->ncols-1]; */
/*                 col->icol = -1; */
/*                 col->replace = replace; */
/*                 col->hdr_key = strdup(hrec->vals[k]); */
/*                 col->number  = bcf_hdr_id2length(hand->hdr_out,BCF_HL_INFO,hdr_id); */
/*                 switch ( bcf_hdr_id2type(hand->hdr_out,BCF_HL_INFO,hdr_id) ) */
/*                 { */
/*                     case BCF_HT_FLAG:   col->setter = vcf_setter_info_flag; break; */
/*                     case BCF_HT_INT:    col->setter = vcf_setter_info_int; break; */
/*                     case BCF_HT_REAL:   col->setter = vcf_setter_info_real; break; */
/*                     case BCF_HT_STR:    col->setter = vcf_setter_info_str; break; */
/*                     default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(hand->hdr_out,BCF_HL_INFO,hdr_id)); */
/*                 } */
/*             } */
/*         } */
/*         else if ( hand->tgts_is_vcf && (!strcasecmp("FORMAT",str.s) || !strcasecmp("FMT",str.s)) ) // All FORMAT fields */
/*         { */
/*             bcf_hdr_t *tgts_hdr = hand->files->readers[1].header; */
/*             if ( force_samples<0 ) force_samples = replace; */
/*             if ( force_samples>=0 && replace!=REPLACE_ALL ) force_samples = replace; */
/*             int j; */
/*             for (j=0; j<tgts_hdr->nhrec; j++) */
/*             { */
/*                 bcf_hrec_t *hrec = tgts_hdr->hrec[j]; */
/*                 if ( hrec->type!=BCF_HL_FMT) continue; */
/*                 int k = bcf_hrec_find_key(hrec,"ID"); */
/*                 assert( k>=0 ); // this should always be true for valid VCFs */
/*                 if ( skip_fmt && khash_str2int_has_key(skip_fmt,hrec->vals[k]) ) continue; */
/*                 tmp.l = 0; */
/*                 bcf_hrec_format(hrec, &tmp); */
/*                 bcf_hdr_append(hand->hdr_out, tmp.s); */
/*                 bcf_hdr_sync(hand->hdr_out); */
/*                 int hdr_id = bcf_hdr_id2int(hand->hdr_out, BCF_DT_ID, hrec->vals[k]); */
/*                 hand->ncols++; hand->cols = (annot_col_t*) realloc(hand->cols,sizeof(annot_col_t)*hand->ncols); */
/*                 annot_col_t *col = &hand->cols[hand->ncols-1]; */
/*                 col->icol = -1; */
/*                 col->replace = replace; */
/*                 col->hdr_key = strdup(hrec->vals[k]); */
/*                 if ( !strcasecmp("GT",col->hdr_key) ) col->setter = vcf_setter_format_gt; */
/*                 else */
/*                     switch ( bcf_hdr_id2type(hand->hdr_out,BCF_HL_FMT,hdr_id) ) */
/*                     { */
/*                         case BCF_HT_INT:    col->setter = vcf_setter_format_int; break; */
/*                         case BCF_HT_REAL:   col->setter = vcf_setter_format_real; break; */
/*                         case BCF_HT_STR:    col->setter = vcf_setter_format_str; has_fmt_str = 1; break; */
/*                         default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(hand->hdr_out,BCF_HL_FMT,hdr_id)); */
/*                     } */
/*             } */
/*         } */
/*         else if ( !strncasecmp("FORMAT/",str.s, 7) || !strncasecmp("FMT/",str.s,4) ) */
/*         { */
/*             char *key = str.s + (!strncasecmp("FMT/",str.s,4) ? 4 : 7); */
/*             if ( force_samples<0 ) force_samples = replace; */
/*             if ( force_samples>=0 && replace!=REPLACE_ALL ) force_samples = replace; */
/*             if ( hand->tgts_is_vcf ) */
/*             { */
/*                 bcf_hrec_t *hrec = bcf_hdr_get_hrec(hand->files->readers[1].header, BCF_HL_FMT, "ID", key, NULL); */
/*                 tmp.l = 0; */
/*                 bcf_hrec_format(hrec, &tmp); */
/*                 bcf_hdr_append(hand->hdr_out, tmp.s); */
/*                 bcf_hdr_sync(hand->hdr_out); */
/*             } */
/*             int hdr_id = bcf_hdr_id2int(hand->hdr_out, BCF_DT_ID, key); */
/*             if ( !bcf_hdr_idinfo_exists(hand->hdr_out,BCF_HL_FMT,hdr_id) ) */
/*                 error("The tag \"%s\" is not defined in %s\n", str.s, hand->targets_fname); */
/*             hand->ncols++; hand->cols = (annot_col_t*) realloc(hand->cols,sizeof(annot_col_t)*hand->ncols); */
/*             annot_col_t *col = &hand->cols[hand->ncols-1]; */
/*             if ( !hand->tgts_is_vcf ) */
/*             { */
/*                 col->icol = icol; */
/*                 icol += bcf_hdr_nsamples(hand->hdr_out) - 1; */
/*             } */
/*             else */
/*                 col->icol = -1; */
/*             col->replace = replace; */
/*             col->hdr_key = strdup(key); */
/*             if ( !strcasecmp("GT",key) ) col->setter = vcf_setter_format_gt; */
/*             else */
/*                 switch ( bcf_hdr_id2type(hand->hdr_out,BCF_HL_FMT,hdr_id) ) */
/*                 { */
/*                     case BCF_HT_INT:    col->setter = hand->tgts_is_vcf ? vcf_setter_format_int  : setter_format_int; break; */
/*                     case BCF_HT_REAL:   col->setter = hand->tgts_is_vcf ? vcf_setter_format_real : setter_format_real; break; */
/*                     case BCF_HT_STR:    col->setter = hand->tgts_is_vcf ? vcf_setter_format_str  : setter_format_str; has_fmt_str = 1; break; */
/*                     default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(hand->hdr_out,BCF_HL_FMT,hdr_id)); */
/*                 } */
/*         } */
/*         else */
/*         { */
/*             if ( replace==REPLACE_EXISTING ) error("Apologies, the -INFO/TAG feature has not been implemented yet.\n"); */
/*             if ( replace==SET_OR_APPEND ) error("Apologies, the =INFO/TAG feature has not been implemented yet.\n"); */
/*             if ( !strncasecmp("INFO/",str.s,5) ) { memmove(str.s,str.s+5,str.l-4); } */
/*             int hdr_id = bcf_hdr_id2int(hand->hdr_out, BCF_DT_ID, str.s); */
/*             if ( !bcf_hdr_idinfo_exists(hand->hdr_out,BCF_HL_INFO,hdr_id) ) */
/*             { */
/*                 if ( hand->tgts_is_vcf ) // reading annotations from a VCF, add a new header line */
/*                 { */
/*                     bcf_hrec_t *hrec = bcf_hdr_get_hrec(hand->files->readers[1].header, BCF_HL_INFO, "ID", str.s, NULL); */
/*                     if ( !hrec ) error("The tag \"%s\" is not defined in %s\n", str.s,hand->files->readers[1].fname); */
/*                     tmp.l = 0; */
/*                     bcf_hrec_format(hrec, &tmp); */
/*                     bcf_hdr_append(hand->hdr_out, tmp.s); */
/*                     bcf_hdr_sync(hand->hdr_out); */
/*                     hdr_id = bcf_hdr_id2int(hand->hdr_out, BCF_DT_ID, str.s); */
/*                 } */
/*                 else */
/*                     error("The tag \"%s\" is not defined in %s\n", str.s, hand->targets_fname); */
/*                 assert( bcf_hdr_idinfo_exists(hand->hdr_out,BCF_HL_INFO,hdr_id) ); */
/*             } */

/*             hand->ncols++; hand->cols = (annot_col_t*) realloc(hand->cols,sizeof(annot_col_t)*hand->ncols); */
/*             annot_col_t *col = &hand->cols[hand->ncols-1]; */
/*             col->icol = icol; */
/*             col->replace = replace; */
/*             col->hdr_key = strdup(str.s); */
/*             col->number  = bcf_hdr_id2length(hand->hdr_out,BCF_HL_INFO,hdr_id); */
/*             switch ( bcf_hdr_id2type(hand->hdr_out,BCF_HL_INFO,hdr_id) ) */
/*             { */
/*                 case BCF_HT_FLAG:   col->setter = hand->tgts_is_vcf ? vcf_setter_info_flag : setter_info_flag; break; */
/*                 case BCF_HT_INT:    col->setter = hand->tgts_is_vcf ? vcf_setter_info_int  : setter_info_int; break; */
/*                 case BCF_HT_REAL:   col->setter = hand->tgts_is_vcf ? vcf_setter_info_real : setter_info_real; break; */
/*                 case BCF_HT_STR:    col->setter = hand->tgts_is_vcf ? vcf_setter_info_str  : setter_info_str; break; */
/*                 default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(hand->hdr_out,BCF_HL_INFO,hdr_id)); */
/*             } */
/*         } */
/*         if ( !*se ) break; */
/*         ss = ++se; */
/*     } */
/*     free(str.s); */
/*     free(tmp.s); */
/*     if ( hand->to_idx==-1 ) hand->to_idx = hand->from_idx; */
/*     free(hand->columns); */
/*     if ( skip_info ) khash_str2int_destroy_free(skip_info); */
/*     if ( skip_fmt ) khash_str2int_destroy_free(skip_fmt); */
/*     if ( has_fmt_str ) */
/*     { */
/*         int n = bcf_hdr_nsamples(hand->hdr_out); */
/*         if ( hand->tgts_is_vcf && n<bcf_hdr_nsamples(hand->files->readers[1].header) ) n = bcf_hdr_nsamples(hand->files->readers[1].header); */
/*         hand->tmpp  = (char**)malloc(sizeof(char*)*n); */
/*         hand->tmpp2 = (char**)malloc(sizeof(char*)*n); */
/*     } */
/*     if ( force_samples>=0 && hand->tgts_is_vcf ) */
/*         set_samples(hand, hand->files->readers[1].header, hand->hdr, force_samples==REPLACE_ALL ? 0 : 1); */
/* } */

/* static void rename_chrs(struct anno_handler *hand, char *fname) */
/* { */
/*     int n, i; */
/*     char **map = hts_readlist(fname, 1, &n); */
/*     if ( !map ) error("Could not read: %s\n", fname); */
/*     for (i=0; i<n; i++) */
/*     { */
/*         char *ss = map[i]; */
/*         while ( *ss && !isspace(*ss) ) ss++; */
/*         if ( !*ss ) error("Could not parse: %s\n", fname); */
/*         *ss = 0; */
/*         int rid = bcf_hdr_name2id(hand->hdr_out, map[i]); */
/*         bcf_hrec_t *hrec = bcf_hdr_get_hrec(hand->hdr_out, BCF_HL_CTG, "ID", map[i], NULL); */
/*         if ( !hrec ) continue;  // the sequence not present */
/*         int j = bcf_hrec_find_key(hrec, "ID"); */
/*         assert( j>=0 ); */
/*         free(hrec->vals[j]); */
/*         ss++; */
/*         while ( *ss && isspace(*ss) ) ss++; */
/*         char *se = ss; */
/*         while ( *se && !isspace(*se) ) se++; */
/*         *se = 0; */
/*         hrec->vals[j] = strdup(ss); */
/*         hand->hdr_out->id[BCF_DT_CTG][rid].key = hrec->vals[j]; */
/*     } */
/*     for (i=0; i<n; i++) free(map[i]); */
/*     free(map); */
/* } */

/* static void init_data(struct anno_handler *hand) */
/* { */
/*     hand->hdr = hand->files->readers[0].header; */
/*     hand->hdr_out = bcf_hdr_dup(hand->hdr); */

/*     if ( hand->remove_annots ) init_remove_annots(hand); */
/*     if ( hand->header_fname ) init_header_lines(hand); */
/*     if ( hand->targets_fname && hand->tgts_is_vcf ) */
/*     { */
/*         // reading annots from a VCF */
/*         if ( !bcf_sr_add_reader(hand->files, hand->targets_fname) ) */
/*             error("Failed to open %s: %s\n", hand->targets_fname,bcf_sr_strerror(hand->files->errnum)); */
/*     } */
/*     if ( hand->columns ) init_columns(hand); */
/*     if ( hand->targets_fname && !hand->tgts_is_vcf ) */
/*     { */
/*         if ( !hand->columns ) error("The -c option not given\n"); */
/*         if ( hand->chr_idx==-1 ) error("The -c CHROM option not given\n"); */
/*         if ( hand->from_idx==-1 ) error("The -c POS option not given\n"); */
/*         if ( hand->to_idx==-1 ) hand->to_idx = -hand->from_idx - 1; */

/*         hand->tgts = bcf_sr_regions_init(hand->targets_fname,1,hand->chr_idx,hand->from_idx,hand->to_idx); */
/*         if ( !hand->tgts ) error("Could not initialize the annotation file: %s\n", hand->targets_fname); */
/*         if ( !hand->tgts->tbx ) error("Expected tabix-indexed annotation file: %s\n", hand->targets_fname); */
/*     } */
/*     hand->vcmp = vcmp_init(); */

/*     if ( hand->filter_str ) */
/*         hand->filter = filter_init(hand->hdr, hand->filter_str); */

/*     if ( hand->set_ids_fmt ) */
/*     { */
/*         if ( hand->set_ids_fmt[0]=='+' ) { hand->set_ids_replace = 0; hand->set_ids_fmt++; } */
/*         hand->set_ids = convert_init(hand->hdr_out, NULL, 0, hand->set_ids_fmt); */
/*     } */

/*     if ( hand->mark_sites ) */
/*     { */
/*         if ( !hand->targets_fname ) error("The -a option not given\n"); */
/*         if ( hand->tgts_is_vcf ) error("Apologies, this has not been implemented yet: -a is a VCF\n");  // very easy to add.. */
/*         bcf_hdr_printf(hand->hdr_out,"##INFO=<ID=%s,Number=0,Type=Flag,Description=\"Sites %slisted in %s\">", */
/*             hand->mark_sites,hand->mark_sites_logic==MARK_LISTED?"":"not ",hand->mark_sites); */
/*     } */

/*      if (hand->record_cmd_line) bcf_hdr_append_version(hand->hdr_out, hand->argc, hand->argv, "bcftools_annotate"); */
/*     if ( !hand->drop_header ) */
/*     { */
/*         if ( hand->rename_chrs ) rename_chrs(hand, hand->rename_chrs); */

/*         hand->out_fh = hts_open(hand->output_fname,hts_bcf_wmode(hand->output_type)); */
/*         if ( hand->out_fh == NULL ) error("Can't write to \"%s\": %s\n", hand->output_fname, strerror(errno)); */
/*         if ( hand->n_threads ) hts_set_threads(hand->out_fh, hand->n_threads); */
/*         bcf_hdr_write(hand->out_fh, hand->hdr_out); */
/*     } */
/* } */

/* static void destroy_data(struct anno_handler *hand) */
/* { */
/*     int i; */
/*     for (i=0; i<hand->nrm; i++) free(hand->rm[i].key); */
/*     free(hand->rm); */
/*     if ( hand->hdr_out ) bcf_hdr_destroy(hand->hdr_out); */
/*     if (hand->vcmp) vcmp_destroy(hand->vcmp); */
/*     for (i=0; i<hand->ncols; i++) */
/*         free(hand->cols[i].hdr_key); */
/*     free(hand->cols); */
/*     for (i=0; i<hand->malines; i++) */
/*     { */
/*         free(hand->alines[i].cols); */
/*         free(hand->alines[i].als); */
/*         free(hand->alines[i].line.s); */
/*     } */
/*     free(hand->alines); */
/*     if ( hand->tgts ) bcf_sr_regions_destroy(hand->tgts); */
/*     free(hand->tmpks.s); */
/*     free(hand->tmpi); */
/*     free(hand->tmpf); */
/*     free(hand->tmps); */
/*     free(hand->tmpp); */
/*     free(hand->tmpi2); */
/*     free(hand->tmpf2); */
/*     free(hand->tmps2); */
/*     free(hand->tmpp2); */
/*     free(hand->tmpi3); */
/*     free(hand->tmpf3); */
/*     if ( hand->set_ids ) */
/*         convert_destroy(hand->set_ids); */
/*     if ( hand->filter ) */
/*         filter_destroy(hand->filter); */
/*     if (hand->out_fh) hts_close(hand->out_fh); */
/*     free(hand->sample_map); */
/* } */

/* static void buffer_annot_lines(struct anno_handler *hand, bcf1_t *line, int start_pos, int end_pos) */
/* { */
/*     if ( hand->nalines && hand->alines[0].rid != line->rid ) hand->nalines = 0; */

/*     int i = 0; */
/*     while ( i<hand->nalines ) */
/*     { */
/*         if ( line->pos > hand->alines[i].end ) */
/*         { */
/*             hand->nalines--; */
/*             if ( hand->nalines && i<hand->nalines ) */
/*             { */
/*                 annot_line_t tmp = hand->alines[i]; */
/*                 memmove(&hand->alines[i],&hand->alines[i+1],(hand->nalines-i)*sizeof(annot_line_t)); */
/*                 hand->alines[hand->nalines] = tmp; */
/*             } */
/*         } */
/*         else i++; */
/*     } */

/*     if ( hand->ref_idx==-1 && hand->nalines ) return; */

/*     while ( !bcf_sr_regions_overlap(hand->tgts, bcf_seqname(hand->hdr,line), start_pos,end_pos) ) */
/*     { */
/*         hand->nalines++; */
/*         hts_expand0(annot_line_t,hand->nalines,hand->malines,hand->alines); */
/*         annot_line_t *tmp = &hand->alines[hand->nalines-1]; */
/*         tmp->rid   = line->rid; */
/*         tmp->start = hand->tgts->start; */
/*         tmp->end   = hand->tgts->end; */
/*         tmp->line.l = 0; */
/*         kputs(hand->tgts->line.s, &tmp->line); */
/*         char *s = tmp->line.s; */
/*         tmp->ncols = 1; */
/*         hts_expand(char*,tmp->ncols,tmp->mcols,tmp->cols); */
/*         tmp->cols[0] = s; */
/*         while ( *s ) */
/*         { */
/*             if ( *s=='\t' ) */
/*             { */
/*                 tmp->ncols++; */
/*                 hts_expand(char*,tmp->ncols,tmp->mcols,tmp->cols); */
/*                 tmp->cols[tmp->ncols-1] = s+1; */
/*                 *s = 0; */
/*             } */
/*             s++; */
/*         } */
/*         if ( hand->ref_idx != -1 ) */
/*         { */
/*             if ( hand->ref_idx >= tmp->ncols )  */
/*                 error("Could not parse the line, expected %d+ columns, found %d:\n\t%s\n",hand->ref_idx+1,tmp->ncols,hand->tgts->line.s); */
/*             if ( hand->alt_idx >= tmp->ncols ) */
/*                 error("Could not parse the line, expected %d+ columns, found %d:\n\t%s\n",hand->alt_idx+1,tmp->ncols,hand->tgts->line.s); */
/*             tmp->nals = 2; */
/*             hts_expand(char*,tmp->nals,tmp->mals,tmp->als); */
/*             tmp->als[0] = tmp->cols[hand->ref_idx]; */
/*             tmp->als[1] = s = tmp->cols[hand->alt_idx]; */
/*             while ( *s ) */
/*             { */
/*                 if ( *s==',' ) */
/*                 { */
/*                     tmp->nals++; */
/*                     hts_expand(char*,tmp->nals,tmp->mals,tmp->als); */
/*                     tmp->als[tmp->nals-1] = s+1; */
/*                     *s = 0; */
/*                 } */
/*                 s++; */
/*             } */
/*             int iseq = hand->tgts->iseq; */
/*             if ( bcf_sr_regions_next(hand->tgts)<0 || hand->tgts->iseq!=iseq ) break; */
/*         } */
/*         else break; */
/*     } */
/* } */

/* static void annotate(struct anno_handler *hand, bcf1_t *line) */
/* { */
/*     int i, j; */
/*     for (i=0; i<hand->nrm; i++) */
/*         hand->rm[i].handler(hand, line, &hand->rm[i]); */

/*     if ( hand->tgts ) */
/*     { */
/*         // Buffer annotation lines. When multiple ALT alleles are present in the */
/*         // annotation file, at least one must match one of the VCF alleles. */
/*         int len = 0; */
/*         bcf_get_variant_types(line); */
/*         for (i=1; i<line->n_allele; i++) */
/*             if ( len > line->d.var[i].n ) len = line->d.var[i].n; */
/*         int end_pos = len<0 ? line->pos - len : line->pos; */
/*         buffer_annot_lines(hand, line, line->pos, end_pos); */
/*         for (i=0; i<hand->nalines; i++) */
/*         { */
/*             if ( line->pos > hand->alines[i].end || end_pos < hand->alines[i].start ) continue; */
/*             if ( hand->ref_idx != -1 ) */
/*             { */
/*                 if ( vcmp_set_ref(hand->vcmp, line->d.allele[0], hand->alines[i].als[0]) < 0 ) continue;   // refs not compatible */
/*                 for (j=1; j<hand->alines[i].nals; j++) */
/*                 { */
/*                     if ( line->n_allele==1 && hand->alines[i].als[j][0]=='.' && hand->alines[i].als[j][1]==0 ) break;   // no ALT allele in VCF and annot file has "." */
/*                     if ( vcmp_find_allele(hand->vcmp, line->d.allele+1, line->n_allele - 1, hand->alines[i].als[j]) >= 0 ) break; */
/*                 } */
/*                 if ( j==hand->alines[i].nals ) continue;    // none of the annot alleles present in VCF's ALT */
/*             } */
/*             break; */
/*         } */

/*         if ( i<hand->nalines ) */
/*         { */
/*             // there is a matching line */
/*             for (j=0; j<hand->ncols; j++) */
/*                 if ( hand->cols[j].setter(hand,line,&hand->cols[j],&hand->alines[i]) ) */
/*                     error("fixme: Could not set %s at %s:%d\n", hand->cols[j].hdr_key,bcf_seqname(hand->hdr,line),line->pos+1); */

/*         } */

/*         if ( hand->mark_sites ) */
/*         { */
/*             // ideally, we'd like to be far more general than this in future, see https://github.com/samtools/bcftools/issues/87 */
/*             if ( hand->mark_sites_logic==MARK_LISTED ) */
/*                 bcf_update_info_flag(hand->hdr_out,line,hand->mark_sites,NULL,i<hand->nalines?1:0); */
/*             else */
/*                 bcf_update_info_flag(hand->hdr_out,line,hand->mark_sites,NULL,i<hand->nalines?0:1); */
/*         } */
/*     } */
/*     else if ( hand->files->nreaders == 2 && bcf_sr_has_line(hand->files,1) ) */
/*     { */
/*         bcf1_t *aline = bcf_sr_get_line(hand->files,1); */
/*         for (j=0; j<hand->ncols; j++) */
/*             if ( hand->cols[j].setter(hand,line,&hand->cols[j],aline) ) */
/*                 error("fixme: Could not set %s at %s:%d\n", hand->cols[j].hdr_key,bcf_seqname(hand->hdr,line),line->pos+1); */
/*     } */
/*     if ( hand->set_ids ) */
/*     { */
/*         hand->tmpks.l = 0; */
/*         convert_line(hand->set_ids, line, &hand->tmpks); */
/*         if ( hand->tmpks.l ) */
/*         { */
/*             int replace = 0; */
/*             if ( hand->set_ids_replace ) replace = 1; */
/*             else if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) ) replace = 1; */
/*             if ( replace ) */
/*                 bcf_update_id(hand->hdr_out,line,hand->tmpks.s); */
/*         } */
/*     } */
/* } */

/* static void usage(struct anno_handler *hand) */
/* { */
/*     fprintf(stderr, "\n"); */
/*     fprintf(stderr, "About:   Annotate and edit VCF/BCF files.\n"); */
/*     fprintf(stderr, "Usage:   bcftools annotate [options] <in.vcf.gz>\n"); */
/*     fprintf(stderr, "\n"); */
/*     fprintf(stderr, "Options:\n"); */
/*     fprintf(stderr, "   -a, --annotations <file>       VCF file or tabix-indexed file with annotations: CHR\\tPOS[\\tVALUE]+\n"); */
/*     fprintf(stderr, "   -c, --columns <list>           list of columns in the annotation file, e.g. CHROM,POS,REF,ALT,-,INFO/TAG. See man page for details\n"); */
/*     fprintf(stderr, "   -e, --exclude <expr>           exclude sites for which the expression is true (see man page for details)\n"); */
/*     fprintf(stderr, "   -h, --header-lines <file>      lines which should be appended to the VCF header\n"); */
/*     fprintf(stderr, "   -I, --set-id [+]<format>       set ID column, see man page for details\n"); */
/*     fprintf(stderr, "   -i, --include <expr>           select sites for which the expression is true (see man page for details)\n"); */
/*     fprintf(stderr, "   -m, --mark-sites [+-]<tag>     add INFO/tag flag to sites which are (\"+\") or are not (\"-\") listed in the -a file\n"); */
/*     fprintf(stderr, "       --no-version               do not append version and command line to the header\n"); */
/*     fprintf(stderr, "   -o, --output <file>            write output to a file [standard output]\n"); */
/*     fprintf(stderr, "   -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n"); */
/*     fprintf(stderr, "   -r, --regions <region>         restrict to comma-separated list of regions\n"); */
/*     fprintf(stderr, "   -R, --regions-file <file>      restrict to regions listed in a file\n"); */
/*     fprintf(stderr, "       --rename-chrs <file>       rename sequences according to map file: from\\tto\n"); */
/*     fprintf(stderr, "   -s, --samples [^]<list>        comma separated list of samples to annotate (or exclude with \"^\" prefix)\n"); */
/*     fprintf(stderr, "   -S, --samples-file [^]<file>   file of samples to annotate (or exclude with \"^\" prefix)\n"); */
/*     fprintf(stderr, "   -x, --remove <list>            list of annotations to remove (e.g. ID,INFO/DP,FORMAT/DP,FILTER). See man page for details\n"); */
/*     fprintf(stderr, "       --threads <int>            number of extra output compression threads [0]\n"); */
/*     fprintf(stderr, "\n"); */
/*     exit(1); */
/* } */

/* int main_vcfannotate(int argc, char *argv[]) */
/* { */
/*     int c; */
/*     struct anno_handler *hand  = (struct anno_handler*) calloc(1,sizeof(struct anno_handler)); */
/*     hand->argc    = argc; hand->argv = argv; */
/*     hand->files   = bcf_sr_init(); */
/*     hand->output_fname = "-"; */
/*     hand->output_type = FT_VCF; */
/*     hand->n_threads = 0; */
/*     hand->record_cmd_line = 1; */
/*     hand->ref_idx = hand->alt_idx = hand->chr_idx = hand->from_idx = hand->to_idx = -1; */
/*     hand->set_ids_replace = 1; */
/*     int regions_is_file = 0; */

/*     static struct option loptions[] = */
/*     { */
/*         {"mark-sites",required_argument,NULL,'m'}, */
/*         {"set-id",required_argument,NULL,'I'}, */
/*         {"output",required_argument,NULL,'o'}, */
/*         {"output-type",required_argument,NULL,'O'}, */
/*         {"threads",required_argument,NULL,9}, */
/*         {"annotations",required_argument,NULL,'a'}, */
/*         {"include",required_argument,NULL,'i'}, */
/*         {"exclude",required_argument,NULL,'e'}, */
/*         {"regions",required_argument,NULL,'r'}, */
/*         {"regions-file",required_argument,NULL,'R'}, */
/*         {"remove",required_argument,NULL,'x'}, */
/*         {"columns",required_argument,NULL,'c'}, */
/*         {"rename-chrs",required_argument,NULL,1}, */
/*         {"header-lines",required_argument,NULL,'h'}, */
/*         {"samples",required_argument,NULL,'s'}, */
/*         {"samples-file",required_argument,NULL,'S'}, */
/*         {"no-version",no_argument,NULL,8}, */
/*         {NULL,0,NULL,0} */
/*     }; */
/*     while ((c = getopt_long(argc, argv, "h:?o:O:r:R:a:x:c:i:e:S:s:I:m:",loptions,NULL)) >= 0) */
/*     { */
/*         switch (c) { */
/*             case 'm':  */
/*                 hand->mark_sites_logic = MARK_LISTED; */
/*                 if ( optarg[0]=='+' ) hand->mark_sites = optarg+1; */
/*                 else if ( optarg[0]=='-' ) { hand->mark_sites = optarg+1; hand->mark_sites_logic = MARK_UNLISTED; } */
/*                 else hand->mark_sites = optarg;  */
/*                 break; */
/*             case 'I': hand->set_ids_fmt = optarg; break; */
/*             case 's': hand->sample_names = optarg; break; */
/*             case 'S': hand->sample_names = optarg; hand->sample_is_file = 1; break; */
/*             case 'c': hand->columns = strdup(optarg); break; */
/*             case 'o': hand->output_fname = optarg; break; */
/*             case 'O': */
/*                 switch (optarg[0]) { */
/*                     case 'b': hand->output_type = FT_BCF_GZ; break; */
/*                     case 'u': hand->output_type = FT_BCF; break; */
/*                     case 'z': hand->output_type = FT_VCF_GZ; break; */
/*                     case 'v': hand->output_type = FT_VCF; break; */
/*                     default: error("The output type \"%s\" not recognised\n", optarg); */
/*                 }; */
/*                 break; */
/*             case 'e': hand->filter_str = optarg; hand->filter_logic |= FLT_EXCLUDE; break; */
/*             case 'i': hand->filter_str = optarg; hand->filter_logic |= FLT_INCLUDE; break; */
/*             case 'x': hand->remove_annots = optarg; break; */
/*             case 'a': hand->targets_fname = optarg; break; */
/*             case 'r': hand->regions_list = optarg; break; */
/*             case 'R': hand->regions_list = optarg; regions_is_file = 1; break; */
/*             case 'h': hand->header_fname = optarg; break; */
/*             case  1 : hand->rename_chrs = optarg; break; */
/*             case  9 : hand->n_threads = strtol(optarg, 0, 0); break; */
/*             case  8 : hand->record_cmd_line = 0; break; */
/*             case '?': usage(hand); break; */
/*             default: error("Unknown argument: %s\n", optarg); */
/*         } */
/*     } */

/*     char *fname = NULL; */
/*     if ( optind>=argc ) */
/*     { */
/*         if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin */
/*         else usage(hand); */
/*     } */
/*     else fname = argv[optind]; */

/*     if ( hand->regions_list ) */
/*     { */
/*         if ( bcf_sr_set_regions(hand->files, hand->regions_list, regions_is_file)<0 ) */
/*             error("Failed to read the regions: %s\n", hand->regions_list); */
/*     } */
/*     if ( hand->targets_fname ) */
/*     { */
/*         htsFile *fp = hts_open(hand->targets_fname,"r");  */
/*         htsFormat type = *hts_get_format(fp); */
/*         hts_close(fp); */

/*         if ( type.format==vcf || type.format==bcf ) */
/*         { */
/*             hand->tgts_is_vcf = 1; */
/*             hand->files->require_index = 1; */
/*             hand->files->collapse |= COLLAPSE_SOME; */
/*         } */
/*     } */
/*     if ( !bcf_sr_add_reader(hand->files, fname) ) error("Failed to open %s: %s\n", fname,bcf_sr_strerror(hand->files->errnum)); */

/*     init_data(hand); */
/*     while ( bcf_sr_next_line(hand->files) ) */
/*     { */
/*         if ( !bcf_sr_has_line(hand->files,0) ) continue; */
/*         bcf1_t *line = bcf_sr_get_line(hand->files,0); */
/*         if ( line->errcode ) error("Encountered error, cannot proceed. Please check the error output above.\n"); */
/*         if ( hand->filter ) */
/*         { */
/*             int pass = filter_test(hand->filter, line, NULL); */
/*             if ( hand->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1; */
/*             if ( !pass ) continue; */
/*         } */
/*         annotate(hand, line); */
/*         bcf_write1(hand->out_fh, hand->hdr_out, line); */
/*     } */
/*     destroy_data(hand); */
/*     bcf_sr_destroy(hand->files); */
/*     free(hand); */
/*     return 0; */

