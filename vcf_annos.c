// vcf_annos.c - define annotate functions for each tag type
//
#include "utils.h"
#include "anno.h"
#include "plugin.h"
#include "vcmp.h"

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define MARK_LISTED   1
#define MARK_UNLISTED 2

struct anno_aux *anno_handler_init (void)
{
    struct anno_aux *hand = (struct anno_aux*)malloc(sizeof(struct anno_aux));
    hand->hdr = NULL;
    hand->hdr_out = NULL;
    
}

int setter_filter(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    // note: so far this works only with one filter, not a list of filters
    struct annot_line *tab = (struct annot_line*) data;
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

int vcf_setter_filter(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    int i;
    bcf1_t *rec = (bcf1_t*) data;
    bcf_hdr_t *hdr = hand->files->readers[hand->ti].header;
    if ( !(rec->unpacked & BCF_UN_FLT) ) bcf_unpack(rec, BCF_UN_FLT);
    if ( !(line->unpacked & BCF_UN_FLT) ) bcf_unpack(line, BCF_UN_FLT);
    if ( !rec->d.n_flt ) return 0;  // don't overwrite with a missing value
    if ( col->replace==SET_OR_APPEND || col->replace==REPLACE_MISSING )
    {
        if ( col->replace==REPLACE_MISSING && line->d.n_flt ) return 0; // only update missing FILTER
        for (i=0; i<rec->d.n_flt; i++)
        {
            const char *flt = bcf_hdr_int2id(hdr, BCF_DT_ID, rec->d.flt[i]);
            bcf_add_filter(hand->hdr_out,line,bcf_hdr_id2int(hand->hdr_out, BCF_DT_ID, flt));
        }
        return 0;
    }
    hts_expand(int,rec->d.n_flt,hand->mtmpi,hand->tmpi);
    for (i=0; i<rec->d.n_flt; i++)
    {
        const char *flt = bcf_hdr_int2id(hdr, BCF_DT_ID, rec->d.flt[i]);
        hand->tmpi[i] = bcf_hdr_id2int(hand->hdr_out, BCF_DT_ID, flt);
    }
    bcf_update_filter(hand->hdr_out,line,NULL,0);
    bcf_update_filter(hand->hdr_out,line,hand->tmpi,rec->d.n_flt);
    return 0;
}
int setter_id(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
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
    struct annot_line *tab = (struct annot_line*) data;
    if ( tab->cols[col->icol] && tab->cols[col->icol][0]=='.' && !tab->cols[col->icol][1] ) return 0; // don't replace with "."
    if ( col->replace==SET_OR_APPEND ) return bcf_add_id(hand->hdr_out,line,tab->cols[col->icol]);
    if ( col->replace!=REPLACE_MISSING ) return bcf_update_id(hand->hdr_out,line,tab->cols[col->icol]);

    // running with +ID, only update missing ids
    if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) )
        return bcf_update_id(hand->hdr_out,line,tab->cols[col->icol]);
    return 0;
}
int vcf_setter_id(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
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
int setter_qual(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    struct annot_line *tab = (struct annot_line*) data;
    char *str = tab->cols[col->icol];
    if ( str[0]=='.' && str[1]==0 ) return 0;   // empty

    if ( col->replace==REPLACE_MISSING && !bcf_float_is_missing(line->qual) ) return 0;

    line->qual = strtod(str, &str);
    if ( str == tab->cols[col->icol] )
        error("Could not parse %s at %s:%d [%s]\n", col->hdr_key, bcf_seqname(hand->hdr,line),line->pos+1,tab->cols[col->icol]);
    return 0;
}
int vcf_setter_qual(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    if ( bcf_float_is_missing(rec->qual) ) return 0;
    if ( col->replace==REPLACE_MISSING && !bcf_float_is_missing(line->qual) ) return 0;
    line->qual = rec->qual;
    return 0;
}
int setter_info_flag(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    struct annot_line *tab = (struct annot_line*) data;
    char *str = tab->cols[col->icol];
    if ( str[0]=='.' && str[1]==0 ) return 0;

    if ( str[0]=='1' && str[1]==0 ) return bcf_update_info_flag(hand->hdr_out,line,col->hdr_key,NULL,1);
    if ( str[0]=='0' && str[1]==0 ) return bcf_update_info_flag(hand->hdr_out,line,col->hdr_key,NULL,0);
    error("Could not parse %s at %s:%d .. [%s]",col->hdr_key, bcf_seqname(hand->hdr,line), line->pos+1, tab->cols[col->icol]);
    return -1;
}
int vcf_setter_info_flag(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    bcf_hdr_t *hdr = hand->files->readers[hand->ti].header;
    int flag = bcf_get_info_flag(hdr,rec,col->hdr_key,NULL,NULL);
    bcf_update_info_flag(hand->hdr_out,line,col->hdr_key,NULL,flag);
    return 0;
}
int setter_ARinfo_int32(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, int nals, char **als, int ntmpi)
{
    if ( col->number==BCF_VL_A && ntmpi!=nals-1 && (ntmpi!=1 || hand->tmpi[0]!=bcf_int32_missing || hand->tmpi[1]!=bcf_int32_vector_end) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpi, col->hdr_key, bcf_seqname(hand->hdr,line), line->pos+1);
    else if ( col->number==BCF_VL_R && ntmpi!=nals && (ntmpi!=1 || hand->tmpi[0]!=bcf_int32_missing || hand->tmpi[1]!=bcf_int32_vector_end) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpi, col->hdr_key, bcf_seqname(hand->hdr,line), line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    int *map = vcmp_map_ARvalues(hand->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%d\n", bcf_seqname(hand->hdr, line), line->pos +1);

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
int setter_info_int(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    struct annot_line *tab = (struct annot_line*) data;
    char *str = tab->cols[col->icol], *end = str;
    if ( str[0]=='.' && str[1]==0 ) return 0;

    int ntmpi = 0;
    while ( *end )
    {
        int val = strtol(str, &end, 10); 
        if ( end==str )
            error("Could not parse %s at %s:%d .. [%s]\n", col->hdr_key, bcf_seqname(hand->hdr,line),line->pos+1,tab->cols[col->icol]);
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
int vcf_setter_info_int(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    bcf_hdr_t *hdr = hand->files->readers[hand->ti].header;
    int ntmpi = bcf_get_info_int32(hdr, rec, col->hdr_key, &hand->tmpi, &hand->mtmpi);
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
int setter_ARinfo_real(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, int nals, char **als, int ntmpf)
{
    if ( col->number==BCF_VL_A && ntmpf!=nals-1 && (ntmpf!=1 || !bcf_float_is_missing(hand->tmpf[0]) || !bcf_float_is_vector_end(hand->tmpf[0])) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpf,col->hdr_key,bcf_seqname(hand->hdr,line),line->pos+1);
    else if ( col->number==BCF_VL_R && ntmpf!=nals && (ntmpf!=1 || !bcf_float_is_missing(hand->tmpf[0]) || !bcf_float_is_vector_end(hand->tmpf[0])) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpf,col->hdr_key,bcf_seqname(hand->hdr,line),line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    int *map = vcmp_map_ARvalues(hand->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%d\n", bcf_seqname(hand->hdr, line), line->pos +1);

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
int setter_info_real(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    struct annot_line *tab = (struct annot_line*) data;
    char *str = tab->cols[col->icol], *end = str;
    if ( str[0]=='.' && str[1]==0 ) return 0;

    int ntmpf = 0;
    while ( *end )
    {
        double val = strtod(str, &end);
        if ( end==str )
            error("Could not parse %s at %s:%d .. [%s]\n",col->hdr_key, bcf_seqname(hand->hdr,line),line->pos+1,tab->cols[col->icol]);
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
int vcf_setter_info_real(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    bcf_hdr_t *hdr = hand->files->readers[hand->ti].header;
    int ntmpf = bcf_get_info_float(hdr,rec,col->hdr_key,&hand->tmpf,&hand->mtmpf);
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

int copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst)
{
    int ith_src = 0, start_src = 0;    // i-th field in src string
    while ( ith_src<isrc && start_src<src_len )
    {
        if ( src[start_src]==',' ) { ith_src++; }
        start_src++;
    }
    if ( ith_src!=isrc ) return -1; // requested field not found
    int end_src = start_src;
    while ( end_src<src_len && src[end_src] && src[end_src]!=',' ) end_src++;

    int nsrc_cpy = end_src - start_src;
    if ( nsrc_cpy==1 && src[start_src]=='.' ) return 0;   // don't write missing values, dst is already initialized

    int ith_dst = 0, start_dst = 0;
    while ( ith_dst<idst && start_dst<dst->l )
    {
        if ( dst->s[start_dst]==',' ) { ith_dst++; }
        start_dst++;
    }
    if ( ith_dst!=idst ) return -2;
    int end_dst = start_dst;
    while ( end_dst<dst->l && dst->s[end_dst]!=',' ) end_dst++;

    if ( end_dst - start_dst>1 || dst->s[start_dst]!='.' ) return 0;   // do not overwrite non-empty values

    // Now start_dst and end_dst are indexes to the destination memory area
    // which needs to be replaced with nsrc_cpy
    // source bytes, end_dst points just after.
    int ndst_shift = nsrc_cpy - (end_dst - start_dst);
    int ndst_move  = dst->l - end_dst + 1;  // how many bytes must be moved (including \0)
    if ( ndst_shift )
    {
        ks_resize(dst, dst->l + ndst_shift + 1);    // plus \0
        memmove(dst->s+end_dst+ndst_shift, dst->s+end_dst, ndst_move);
    }
    memcpy(dst->s+start_dst, src+start_src, nsrc_cpy);
    dst->l += ndst_shift;
    return 0;
}

int setter_ARinfo_string(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, int nals, char **als)
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
    if ( !map ) error("REF alleles not compatible at %s:%d\n", bcf_seqname(hand->hdr, line), line->pos+1);

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
int setter_info_str(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    struct annot_line *tab = (struct annot_line*) data;
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
int vcf_setter_info_str(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    bcf_hdr_t *hdr = hand->files->readers[hand->ti].header;
    int ntmps = bcf_get_info_string(hdr,rec,col->hdr_key,&hand->tmps,&hand->mtmps);
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
int vcf_setter_format_gt(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    assert(hand->ti > 0);
    bcf_hdr_t *header = hand->files->readers[hand->ti].header;
    int nsrc = bcf_get_genotypes(header, rec, &hand->tmpi, &hand->mtmpi);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    if ( !hand->sample_map )
        return bcf_update_genotypes(hand->hdr_out,line,hand->tmpi,nsrc);

    int i, j, ndst = bcf_get_genotypes(hand->hdr,line,&hand->tmpi2,&hand->mtmpi2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(hand->hdr_out);
    nsrc /= bcf_hdr_nsamples(header);
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
int count_vals(struct annot_line *tab, int icol_beg, int icol_end)
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
int setter_format_int(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    struct annot_line *tab = (struct annot_line*) data;
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

int setter_format_real(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    struct annot_line *tab = (struct annot_line*) data;
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
int setter_format_str(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    struct annot_line *tab = (struct annot_line*) data;
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

int vcf_setter_format_int(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    assert(hand->ti > 0);
    bcf_hdr_t *header = hand->files->readers[hand->ti].header;
    int nsrc = bcf_get_format_int32(header, rec, col->hdr_key, &hand->tmpi, &hand->mtmpi);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    if ( !hand->sample_map )
        return bcf_update_format_int32(hand->hdr_out,line,col->hdr_key,hand->tmpi,nsrc);

    int i, j, ndst = bcf_get_format_int32(hand->hdr,line,col->hdr_key,&hand->tmpi2,&hand->mtmpi2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(hand->hdr_out);
    nsrc /= bcf_hdr_nsamples(header);
    if ( ndst<=0 ) {
        if ( col->replace==REPLACE_EXISTING ) return 0;    // overwrite only if present
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(hand->hdr_out), hand->mtmpi2, hand->tmpi2);
        for (i=0; i<bcf_hdr_nsamples(hand->hdr_out); i++) {
            int32_t *dst = hand->tmpi2 + nsrc*i;
            if ( hand->sample_map[i]==-1 ) {
                dst[0] = bcf_int32_missing;
                for (j=1; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            } else {
                int32_t *src = hand->tmpi + nsrc*hand->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_int32(hand->hdr_out,line,col->hdr_key,hand->tmpi2,nsrc*bcf_hdr_nsamples(hand->hdr_out));
    } else if ( ndst >= nsrc ) {
        for (i=0; i<bcf_hdr_nsamples(hand->hdr_out); i++) {
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
int vcf_setter_format_real(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    bcf_hdr_t *header = hand->files->readers[hand->ti].header;
    int nsrc = bcf_get_format_float(header, rec, col->hdr_key, &hand->tmpf, &hand->mtmpf);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    if ( !hand->sample_map )
        return bcf_update_format_float(hand->hdr_out,line,col->hdr_key,hand->tmpf,nsrc);

    int i, j, ndst = bcf_get_format_float(hand->hdr,line,col->hdr_key,&hand->tmpf2,&hand->mtmpf2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(hand->hdr_out);
    nsrc /= bcf_hdr_nsamples(header);
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
int vcf_setter_format_str(struct anno_aux *hand, bcf1_t *line, struct annot_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    hand->tmpp[0] = hand->tmps;
    bcf_hdr_t *header = hand->files->readers[hand->ti].header;
    int ret = bcf_get_format_string(header,rec,col->hdr_key,&hand->tmpp,&hand->mtmps);
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
