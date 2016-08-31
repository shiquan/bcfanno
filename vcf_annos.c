// vcf_annos.c - define annotate functions for each tag type
//
#include "utils.h"
#include "anno.h"
#include "anno_vcf.h"
#include "vcmp.h"
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/bgzf.h>

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define MARK_LISTED   1
#define MARK_UNLISTED 2

int setter_filter(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    // note: so far this works only with one filter, not a list of filters
    struct anno_line *tab = (struct anno_line*) data;
    if ( tab->cols[col->icol] && tab->cols[col->icol][0]=='.' && !tab->cols[col->icol][1] ) return 0; // don't replace with "."
    hts_expand(int,1,opts->mtmpi,opts->tmpi);
    opts->tmpi[0] = bcf_hdr_id2int(opts->hdr_out, BCF_DT_ID, tab->cols[col->icol]);
    if ( opts->tmpi[0]<0 ) error("The FILTER is not defined in the header: %s\n", tab->cols[col->icol]);
    if ( col->replace==SET_OR_APPEND ) { bcf_add_filter(opts->hdr_out,line,opts->tmpi[0]); return 0; }
    if ( col->replace!=REPLACE_MISSING )
    {
        bcf_update_filter(opts->hdr_out,line,NULL,0);
        bcf_update_filter(opts->hdr_out,line,opts->tmpi,1); 
        return 0; 
    }
    
    // only update missing FILTER
    if ( !(line->unpacked & BCF_UN_FLT) ) bcf_unpack(line, BCF_UN_FLT);
    if ( !line->d.n_flt )
        bcf_update_filter(opts->hdr_out,line,opts->tmpi,1);
    return 0;
}

int vcf_setter_filter(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    int i;
    bcf1_t *rec = (bcf1_t*) data;
    bcf_hdr_t *hdr = opts->files[col->ifile].hdr;
    if ( !(rec->unpacked & BCF_UN_FLT) ) bcf_unpack(rec, BCF_UN_FLT);
    if ( !(line->unpacked & BCF_UN_FLT) ) bcf_unpack(line, BCF_UN_FLT);
    if ( !rec->d.n_flt ) return 0;  // don't overwrite with a missing value
    if ( col->replace==SET_OR_APPEND || col->replace==REPLACE_MISSING )
    {
        if ( col->replace==REPLACE_MISSING && line->d.n_flt ) return 0; // only update missing FILTER
        for (i=0; i<rec->d.n_flt; i++)
        {
            const char *flt = bcf_hdr_int2id(hdr, BCF_DT_ID, rec->d.flt[i]);
            bcf_add_filter(opts->hdr_out,line,bcf_hdr_id2int(opts->hdr_out, BCF_DT_ID, flt));
        }
        return 0;
    }
    hts_expand(int,rec->d.n_flt,opts->mtmpi,opts->tmpi);
    for (i=0; i<rec->d.n_flt; i++)
    {
        const char *flt = bcf_hdr_int2id(hdr, BCF_DT_ID, rec->d.flt[i]);
        opts->tmpi[i] = bcf_hdr_id2int(opts->hdr_out, BCF_DT_ID, flt);
    }
    bcf_update_filter(opts->hdr_out,line,NULL,0);
    bcf_update_filter(opts->hdr_out,line,opts->tmpi,rec->d.n_flt);
    return 0;
}
int setter_id(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
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
    struct anno_line *tab = (struct anno_line*) data;
    if ( tab->cols[col->icol] && tab->cols[col->icol][0]=='.' && !tab->cols[col->icol][1] ) return 0; // don't replace with "."
    if ( col->replace==SET_OR_APPEND ) return bcf_add_id(opts->hdr_out,line,tab->cols[col->icol]);
    if ( col->replace!=REPLACE_MISSING ) return bcf_update_id(opts->hdr_out,line,tab->cols[col->icol]);

    // running with +ID, only update missing ids
    if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) )
        return bcf_update_id(opts->hdr_out,line,tab->cols[col->icol]);
    return 0;
}
int vcf_setter_id(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    if ( rec->d.id && rec->d.id[0]=='.' && !rec->d.id[1] ) return 0;    // don't replace with "."
    if ( col->replace==SET_OR_APPEND ) return bcf_add_id(opts->hdr_out,line,rec->d.id);
    if ( col->replace!=REPLACE_MISSING ) return bcf_update_id(opts->hdr_out,line,rec->d.id);

   // running with +ID, only update missing ids
    if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) )
        return bcf_update_id(opts->hdr_out,line,rec->d.id);
    return 0;
}
int setter_qual(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    struct anno_line *tab = (struct anno_line*) data;
    char *str = tab->cols[col->icol];
    if ( str[0]=='.' && str[1]==0 ) return 0;   // empty

    if ( col->replace==REPLACE_MISSING && !bcf_float_is_missing(line->qual) ) return 0;

    line->qual = strtod(str, &str);
    if ( str == tab->cols[col->icol] )
        error("Could not parse %s at %s:%d [%s]\n", col->hdr_key, bcf_seqname(opts->hdr_out,line),line->pos+1,tab->cols[col->icol]);
    return 0;
}
int vcf_setter_qual(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    if ( bcf_float_is_missing(rec->qual) ) return 0;
    if ( col->replace==REPLACE_MISSING && !bcf_float_is_missing(line->qual) ) return 0;
    line->qual = rec->qual;
    return 0;
}
int setter_info_flag(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    struct anno_line *tab = (struct anno_line*) data;
    char *str = tab->cols[col->icol];
    if ( str[0]=='.' && str[1]==0 ) return 0;

    if ( str[0]=='1' && str[1]==0 ) return bcf_update_info_flag(opts->hdr_out,line,col->hdr_key,NULL,1);
    if ( str[0]=='0' && str[1]==0 ) return bcf_update_info_flag(opts->hdr_out,line,col->hdr_key,NULL,0);
    error("Could not parse %s at %s:%d .. [%s]",col->hdr_key, bcf_seqname(opts->hdr_out,line), line->pos+1, tab->cols[col->icol]);
    return -1;
}
int vcf_setter_info_flag(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    bcf_hdr_t *hdr = opts->files[col->ifile].hdr;
    int flag = bcf_get_info_flag(hdr,rec,col->hdr_key,NULL,NULL);
    bcf_update_info_flag(opts->hdr_out,line,col->hdr_key,NULL,flag);
    return 0;
}
int setter_ARinfo_int32(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, int nals, char **als, int ntmpi)
{
    if ( col->number==BCF_VL_A && ntmpi!=nals-1 && (ntmpi!=1 || opts->tmpi[0]!=bcf_int32_missing || opts->tmpi[1]!=bcf_int32_vector_end) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpi, col->hdr_key, bcf_seqname(opts->hdr_out,line), line->pos+1);
    else if ( col->number==BCF_VL_R && ntmpi!=nals && (ntmpi!=1 || opts->tmpi[0]!=bcf_int32_missing || opts->tmpi[1]!=bcf_int32_vector_end) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpi, col->hdr_key, bcf_seqname(opts->hdr_out,line), line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    int *map = vcmp_map_ARvalues(opts->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%d\n", bcf_seqname(opts->hdr_out, line), line->pos +1);

    // fill in any missing values in the target VCF (or all, if not present)
    int ntmpi2 = bcf_get_info_float(opts->hdr_out, line, col->hdr_key, &opts->tmpi2, &opts->mtmpi2);
    if ( ntmpi2 < ndst ) hts_expand(int32_t,ndst,opts->mtmpi2,opts->tmpi2);

    int i;
    for (i=0; i<ndst; i++)
    {
        if ( map[i]<0 )
        {
            if ( ntmpi2 < ndst ) opts->tmpi2[i] = bcf_int32_missing;
            continue;
        }
        if ( ntmpi2==ndst && col->replace==REPLACE_MISSING
                && opts->tmpi2[i]!=bcf_int32_missing
                && opts->tmpi2[i]!=bcf_int32_vector_end ) continue;

        opts->tmpi2[i] = opts->tmpi[ map[i] ];
    }
    bcf_update_info_int32(opts->hdr_out,line,col->hdr_key,opts->tmpi2,ndst);
    return 0;
}
int setter_info_int(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    struct anno_line *tab = (struct anno_line*) data;
    char *str = tab->cols[col->icol], *end = str;
    if ( str[0]=='.' && str[1]==0 ) return 0;

    int ntmpi = 0;
    while ( *end )
    {
        int val = strtol(str, &end, 10); 
        if ( end==str )
            error("Could not parse %s at %s:%d .. [%s]\n", col->hdr_key, bcf_seqname(opts->hdr_out,line),line->pos+1,tab->cols[col->icol]);
        ntmpi++;
        hts_expand(int32_t,ntmpi,opts->mtmpi,opts->tmpi);
        opts->tmpi[ntmpi-1] = val;
        str = end+1;
    }

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_int32(opts,line,col,tab->nals,tab->als,ntmpi);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_int32(opts->hdr_out, line, col->hdr_key, &opts->tmpi2, &opts->mtmpi2);
        if ( ret>0 && opts->tmpi2[0]!=bcf_int32_missing ) return 0;
    }

    bcf_update_info_int32(opts->hdr_out,line,col->hdr_key,opts->tmpi,ntmpi);
    return 0;
}
int vcf_setter_info_int(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    bcf_hdr_t *hdr = opts->files[col->ifile].hdr;
    int ntmpi = bcf_get_info_int32(hdr, rec, col->hdr_key, &opts->tmpi, &opts->mtmpi);
    // debug_print("ntmpi : %d", ntmpi);
    if ( ntmpi < 0 ) return 0;    // nothing to add

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_int32(opts,line,col,rec->n_allele,rec->d.allele,ntmpi);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_int32(opts->hdr_out, line, col->hdr_key, &opts->tmpi2, &opts->mtmpi2);
        if ( ret>0 && opts->tmpi2[0]!=bcf_int32_missing ) return 0;
    }

    bcf_update_info_int32(opts->hdr_out,line,col->hdr_key,opts->tmpi,ntmpi);
    return 0;
}
int setter_ARinfo_real(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, int nals, char **als, int ntmpf)
{
    if ( col->number==BCF_VL_A && ntmpf!=nals-1 && (ntmpf!=1 || !bcf_float_is_missing(opts->tmpf[0]) || !bcf_float_is_vector_end(opts->tmpf[0])) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpf,col->hdr_key,bcf_seqname(opts->hdr_out,line),line->pos+1);
    else if ( col->number==BCF_VL_R && ntmpf!=nals && (ntmpf!=1 || !bcf_float_is_missing(opts->tmpf[0]) || !bcf_float_is_vector_end(opts->tmpf[0])) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpf,col->hdr_key,bcf_seqname(opts->hdr_out,line),line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    // debug_print("ndst : %d, nals : %d", ndst, nals);
    int *map = vcmp_map_ARvalues(opts->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%d\n", bcf_seqname(opts->hdr_out, line), line->pos +1);

    // fill in any missing values in the target VCF (or all, if not present)
    int ntmpf2 = bcf_get_info_float(opts->hdr_out, line, col->hdr_key, &opts->tmpf2, &opts->mtmpf2);
    if ( ntmpf2 < ndst ) hts_expand(float,ndst,opts->mtmpf2,opts->tmpf2);

    int i;
    for (i=0; i<ndst; i++)
    {
        if ( map[i]<0 )
        {
            if ( ntmpf2 < ndst ) bcf_float_set_missing(opts->tmpf2[i]);
            continue;
        }
        if ( ntmpf2==ndst && col->replace==REPLACE_MISSING
                && !bcf_float_is_missing(opts->tmpf2[i])
                && !bcf_float_is_vector_end(opts->tmpf2[i]) ) continue;

        opts->tmpf2[i] = opts->tmpf[ map[i] ];
    }
    bcf_update_info_float(opts->hdr_out,line,col->hdr_key,opts->tmpf2,ndst);
    return 0;
}
int setter_info_real(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    struct anno_line *tab = (struct anno_line*) data;
    char *str = tab->cols[col->icol], *end = str;
    if ( str[0]=='.' && str[1]==0 ) return 0;

    int ntmpf = 0;
    while ( *end )
    {
        double val = strtod(str, &end);
        if ( end==str )
            error("Could not parse %s at %s:%d .. [%s]\n",col->hdr_key, bcf_seqname(opts->hdr_out,line),line->pos+1,tab->cols[col->icol]);
        ntmpf++;
        hts_expand(float,ntmpf,opts->mtmpf,opts->tmpf);
        opts->tmpf[ntmpf-1] = val;
        str = end+1;
    }

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_real(opts,line,col,tab->nals,tab->als,ntmpf);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_float(opts->hdr_out, line, col->hdr_key, &opts->tmpf2, &opts->mtmpf2);
        if ( ret>0 && !bcf_float_is_missing(opts->tmpf2[0]) ) return 0;
    }

    bcf_update_info_float(opts->hdr_out,line,col->hdr_key,opts->tmpf,ntmpf);
    return 0;
}
int vcf_setter_info_real(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    bcf_hdr_t *hdr = opts->files[col->ifile].hdr;
    int ntmpf = bcf_get_info_float(hdr,rec,col->hdr_key,&opts->tmpf,&opts->mtmpf);
    if ( ntmpf < 0 ) return 0;    // nothing to add

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_real(opts,line,col,rec->n_allele,rec->d.allele,ntmpf);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_float(opts->hdr_out, line, col->hdr_key, &opts->tmpf2, &opts->mtmpf2);
        if ( ret>0 && !bcf_float_is_missing(opts->tmpf2[0]) ) return 0;
    }

    bcf_update_info_float(opts->hdr_out,line,col->hdr_key,opts->tmpf,ntmpf);
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

int setter_ARinfo_string(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, int nals, char **als)
{
    int nsrc = 1, lsrc = 0;
    while ( opts->tmps[lsrc] )
    {
        if ( opts->tmps[lsrc]==',' ) nsrc++;
        lsrc++;
    }
    if ( col->number==BCF_VL_A && nsrc!=nals-1 && (nsrc!=1 || opts->tmps[0]!='.' || opts->tmps[1]!=0 ) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", nsrc,col->hdr_key,bcf_seqname(opts->hdr_out,line),line->pos+1);
    else if ( col->number==BCF_VL_R && nsrc!=nals && (nsrc!=1 || opts->tmps[0]!='.' || opts->tmps[1]!=0 ) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", nsrc,col->hdr_key,bcf_seqname(opts->hdr_out,line),line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    int *map = vcmp_map_ARvalues(opts->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%d\n", bcf_seqname(opts->hdr_out, line), line->pos+1);

    // fill in any missing values in the target VCF (or all, if not present)
    int i, empty = 0, nstr, mstr = opts->tmpks.m;
    nstr = bcf_get_info_string(opts->hdr_out, line, col->hdr_key, &opts->tmpks.s, &mstr); 
    opts->tmpks.m = mstr;
    if ( nstr<0 || (nstr==1 && opts->tmpks.s[0]=='.' && opts->tmpks.s[1]==0) )
    {
        empty = 0;
        opts->tmpks.l = 0;
        kputc('.',&opts->tmpks);
        for (i=1; i<ndst; i++) kputs(",.",&opts->tmpks);
    }
    else opts->tmpks.l = nstr;
    for (i=0; i<ndst; i++)
    {
        if ( map[i]<0 )
        {
            if ( empty ) copy_string_field(".",0,1,&opts->tmpks,i);
            continue;
        }
        if ( col->replace==REPLACE_MISSING )
        {
            // Do not replace filled values. The field must be looked up again because
            // of realloc in copy_string_field
            int n = 0;
            char *str = opts->tmpks.s;
            while ( *str && n<i )
            {
                if ( *str==',' ) n++;
                str++;
            }
            if ( str[0]!='.' || (str[1]!=',' && str[1]!=0) ) continue;  // value already set
        }
        int ret = copy_string_field(opts->tmps,map[i],lsrc,&opts->tmpks,i);
        assert( ret==0 );
    }
    bcf_update_info_string(opts->hdr_out,line,col->hdr_key,opts->tmpks.s);
    return 0;
}
int setter_info_str(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    struct anno_line *tab = (struct anno_line*) data;
    int len = strlen(tab->cols[col->icol]);
    if ( !len ) return 0;
    hts_expand(char,len+1,opts->mtmps,opts->tmps);
    memcpy(opts->tmps,tab->cols[col->icol],len+1);
    if ( opts->tmps[0]=='.' && opts->tmps[1]==0 ) return 0;

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_string(opts,line,col,tab->nals,tab->als);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_string(opts->hdr_out, line, col->hdr_key, &opts->tmps2, &opts->mtmps2);
        if ( ret>0 && (opts->tmps2[0]!='.' || opts->tmps2[1]!=0) ) return 0;
    }

    bcf_update_info_string(opts->hdr_out,line,col->hdr_key,opts->tmps);
    return 0;
}
int vcf_setter_info_str(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    bcf_hdr_t *hdr = opts->files[col->ifile].hdr;
    int ntmps = bcf_get_info_string(hdr,rec,col->hdr_key,&opts->tmps,&opts->mtmps);
    if ( ntmps < 0 ) return 0;    // nothing to add

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_string(opts,line,col,rec->n_allele,rec->d.allele);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_string(opts->hdr_out, line, col->hdr_key, &opts->tmps2, &opts->mtmps2);
        if ( ret>0 && (opts->tmps2[0]!='.' || opts->tmps2[1]!=0) ) return 0;
    }

    bcf_update_info_string(opts->hdr_out,line,col->hdr_key,opts->tmps);
    return 0;
}
int vcf_setter_format_gt(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;

    bcf_hdr_t *header = opts->hdr_out;
    int nsrc = bcf_get_genotypes(header, rec, &opts->tmpi, &opts->mtmpi);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    if ( !opts->sample_map )
        return bcf_update_genotypes(opts->hdr_out,line,opts->tmpi,nsrc);

    int i, j, ndst = bcf_get_genotypes(opts->hdr_out,line,&opts->tmpi2,&opts->mtmpi2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(opts->hdr_out);
    nsrc /= bcf_hdr_nsamples(header);
    if ( ndst<=0 )  // field not present in dst file
    {
        if ( col->replace==REPLACE_EXISTING ) return 0;
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(opts->hdr_out), opts->mtmpi2, opts->tmpi2);
        for (i=0; i<bcf_hdr_nsamples(opts->hdr_out); i++)
        {
            int32_t *dst = opts->tmpi2 + nsrc*i;
            if ( opts->sample_map[i]==-1 )
            {
                dst[0] = bcf_gt_missing;
                for (j=1; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = opts->tmpi + nsrc*opts->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_genotypes(opts->hdr_out,line,opts->tmpi2,nsrc*bcf_hdr_nsamples(opts->hdr_out));
    }
    else if ( ndst >= nsrc )     
    {
        for (i=0; i<bcf_hdr_nsamples(opts->hdr_out); i++)
        {
            if ( opts->sample_map[i]==-1 ) continue;
            int32_t *src = opts->tmpi  + nsrc*opts->sample_map[i];
            int32_t *dst = opts->tmpi2 + ndst*i;
            if ( col->replace==REPLACE_EXISTING && bcf_gt_is_missing(dst[0]) ) continue;
            if ( col->replace==REPLACE_MISSING  && !bcf_gt_is_missing(dst[0]) ) continue;
            for (j=0; j<nsrc; j++) dst[j] = src[j];
            for (; j<ndst; j++) dst[j] = bcf_int32_vector_end;
        }
        return bcf_update_genotypes(opts->hdr_out,line,opts->tmpi2,ndst*bcf_hdr_nsamples(opts->hdr_out));
    }
    else    // ndst < nsrc
    {
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(opts->hdr_out), opts->mtmpi3, opts->tmpi3);
        for (i=0; i<bcf_hdr_nsamples(opts->hdr_out); i++)
        {
            int32_t *ori = opts->tmpi2 + ndst*i;
            int32_t *dst = opts->tmpi3 + nsrc*i;
            int keep_ori = 0;
            if ( opts->sample_map[i]==-1 ) keep_ori = 1;
            else if ( col->replace==REPLACE_EXISTING && bcf_gt_is_missing(ori[0]) ) keep_ori = 1;
            else if ( col->replace==REPLACE_MISSING  && !bcf_gt_is_missing(ori[0]) ) keep_ori = 1;
            if ( keep_ori )
            {
                for (j=0; j<ndst; j++) dst[j] = ori[j];
                for (; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = opts->tmpi + nsrc*opts->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_genotypes(opts->hdr_out,line,opts->tmpi3,nsrc*bcf_hdr_nsamples(opts->hdr_out));
    }
}
int count_vals(struct anno_line *tab, int icol_beg, int icol_end)
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
int setter_format_int(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    struct anno_line *tab = (struct anno_line*) data;
    int nsmpl = bcf_hdr_nsamples(opts->hdr_out);
    assert( col->icol+nsmpl <= tab->ncols );
    int nvals = count_vals(tab,col->icol,col->icol+nsmpl);
    assert( nvals>0 );
    hts_expand(int32_t,nvals*nsmpl,opts->mtmpi,opts->tmpi);

    int icol = col->icol, ismpl;
    for (ismpl=0; ismpl<nsmpl; ismpl++)
    {
        int32_t *ptr = opts->tmpi + ismpl*nvals;
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
                error("Could not parse %s at %s:%d .. [%s]\n", col->hdr_key,bcf_seqname(opts->hdr_out,line),line->pos+1,tab->cols[col->icol]);

            ival++;
            str = *end ? end+1 : end;
        }
        while ( ival<nvals ) ptr[ival++] = bcf_int32_vector_end;
        icol++;
    }
    return bcf_update_format_int32(opts->hdr_out,line,col->hdr_key,opts->tmpi,nsmpl*nvals);
}

int setter_format_real(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    struct anno_line *tab = (struct anno_line*) data;
    int nsmpl = bcf_hdr_nsamples(opts->hdr_out);
    assert( col->icol+nsmpl <= tab->ncols );
    int nvals = count_vals(tab,col->icol,col->icol+nsmpl);
    assert( nvals>0 );
    hts_expand(float,nvals*nsmpl,opts->mtmpf,opts->tmpf);

    int icol = col->icol, ismpl;
    for (ismpl=0; ismpl<nsmpl; ismpl++)
    {
        float *ptr = opts->tmpf + ismpl*nvals;
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
                error("Could not parse %s at %s:%d .. [%s]\n", col->hdr_key,bcf_seqname(opts->hdr_out,line),line->pos+1,tab->cols[col->icol]);

            ival++;
            str = *end ? end+1 : end;
        }
        while ( ival<nvals ) { bcf_float_set_vector_end(ptr[ival]); ival++; }
        icol++;
    }
    return bcf_update_format_float(opts->hdr_out,line,col->hdr_key,opts->tmpf,nsmpl*nvals);
}
int setter_format_str(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    struct anno_line *tab = (struct anno_line*) data;
    int nsmpl = bcf_hdr_nsamples(opts->hdr_out);
    assert( col->icol+nsmpl <= tab->ncols );

    int i, max_len = 0;
    for (i=col->icol; i<col->icol+nsmpl; i++)
    {
        int len = strlen(tab->cols[i]);
        if ( max_len < len ) max_len = len;
    }
    hts_expand(char,max_len*nsmpl,opts->mtmps,opts->tmps);

    int icol = col->icol, ismpl;
    for (ismpl=0; ismpl<nsmpl; ismpl++)
    {
        char *ptr = opts->tmps + ismpl*max_len;
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
    return bcf_update_format_char(opts->hdr_out,line,col->hdr_key,opts->tmps,nsmpl*max_len);
}

int vcf_setter_format_int(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;

    bcf_hdr_t *header = opts->hdr_out;
    int nsrc = bcf_get_format_int32(header, rec, col->hdr_key, &opts->tmpi, &opts->mtmpi);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    if ( !opts->sample_map )
        return bcf_update_format_int32(opts->hdr_out,line,col->hdr_key,opts->tmpi,nsrc);

    int i, j, ndst = bcf_get_format_int32(opts->hdr_out,line,col->hdr_key,&opts->tmpi2,&opts->mtmpi2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(opts->hdr_out);
    nsrc /= bcf_hdr_nsamples(header);
    if ( ndst<=0 ) {
        if ( col->replace==REPLACE_EXISTING ) return 0;    // overwrite only if present
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(opts->hdr_out), opts->mtmpi2, opts->tmpi2);
        for (i=0; i<bcf_hdr_nsamples(opts->hdr_out); i++) {
            int32_t *dst = opts->tmpi2 + nsrc*i;
            if ( opts->sample_map[i]==-1 ) {
                dst[0] = bcf_int32_missing;
                for (j=1; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            } else {
                int32_t *src = opts->tmpi + nsrc*opts->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_int32(opts->hdr_out,line,col->hdr_key,opts->tmpi2,nsrc*bcf_hdr_nsamples(opts->hdr_out));
    } else if ( ndst >= nsrc ) {
        for (i=0; i<bcf_hdr_nsamples(opts->hdr_out); i++) {
            if ( opts->sample_map[i]==-1 ) continue;
            int32_t *src = opts->tmpi  + nsrc*opts->sample_map[i];
            int32_t *dst = opts->tmpi2 + ndst*i;
            if ( col->replace==REPLACE_EXISTING && dst[0]==bcf_int32_missing ) continue;
            if ( col->replace==REPLACE_MISSING  && dst[0]!=bcf_int32_missing ) continue;
            for (j=0; j<nsrc; j++) dst[j] = src[j];
            for (; j<ndst; j++) dst[j] = bcf_int32_vector_end;
        }
        return bcf_update_format_int32(opts->hdr_out,line,col->hdr_key,opts->tmpi2,ndst*bcf_hdr_nsamples(opts->hdr_out));
    }
    else    // ndst < nsrc
    {
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(opts->hdr_out), opts->mtmpi3, opts->tmpi3);
        for (i=0; i<bcf_hdr_nsamples(opts->hdr_out); i++)
        {
            int32_t *ori = opts->tmpi2 + ndst*i;
            int32_t *dst = opts->tmpi3 + nsrc*i;
            int keep_ori = 0;
            if ( opts->sample_map[i]==-1 ) keep_ori = 1;
            else if ( col->replace==REPLACE_EXISTING && ori[0]==bcf_int32_missing ) keep_ori = 1;
            else if ( col->replace==REPLACE_MISSING  && ori[0]!=bcf_int32_missing ) keep_ori = 1;
            if ( keep_ori )
            {
                for (j=0; j<ndst; j++) dst[j] = ori[j];
                for (; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = opts->tmpi + nsrc*opts->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_int32(opts->hdr_out,line,col->hdr_key,opts->tmpi3,nsrc*bcf_hdr_nsamples(opts->hdr_out));
    }
}
int vcf_setter_format_real(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    bcf_hdr_t *header = opts->hdr_out;
    int nsrc = bcf_get_format_float(header, rec, col->hdr_key, &opts->tmpf, &opts->mtmpf);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    if ( !opts->sample_map )
        return bcf_update_format_float(opts->hdr_out,line,col->hdr_key,opts->tmpf,nsrc);

    int i, j, ndst = bcf_get_format_float(opts->hdr_out,line,col->hdr_key,&opts->tmpf2,&opts->mtmpf2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(opts->hdr_out);
    nsrc /= bcf_hdr_nsamples(header);
    if ( ndst<=0 )
    {
        if ( col->replace==REPLACE_EXISTING ) return 0;    // overwrite only if present
        hts_expand(float, nsrc*bcf_hdr_nsamples(opts->hdr_out), opts->mtmpf2, opts->tmpf2);
        for (i=0; i<bcf_hdr_nsamples(opts->hdr_out); i++)
        {
            float *dst = opts->tmpf2 + nsrc*i;
            if ( opts->sample_map[i]==-1 )
            {
                bcf_float_set_missing(dst[0]);
                for (j=1; j<nsrc; j++) bcf_float_set_vector_end(dst[j]);
            }
            else
            {
                float *src = opts->tmpf + nsrc*opts->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_float(opts->hdr_out,line,col->hdr_key,opts->tmpf2,nsrc*bcf_hdr_nsamples(opts->hdr_out));
    }
    else if ( ndst >= nsrc )     
    {
        for (i=0; i<bcf_hdr_nsamples(opts->hdr_out); i++)
        {
            if ( opts->sample_map[i]==-1 ) continue;
            float *src = opts->tmpf  + nsrc*opts->sample_map[i];
            float *dst = opts->tmpf2 + ndst*i;
            if ( col->replace==REPLACE_EXISTING && bcf_float_is_missing(dst[0]) ) continue;
            if ( col->replace==REPLACE_MISSING  && !bcf_float_is_missing(dst[0]) ) continue;
            for (j=0; j<nsrc; j++) dst[j] = src[j];
            for (; j<ndst; j++) bcf_float_set_vector_end(dst[j]);
        }
        return bcf_update_format_float(opts->hdr_out,line,col->hdr_key,opts->tmpf2,ndst*bcf_hdr_nsamples(opts->hdr_out));
    }
    else    // ndst < nsrc
    {
        hts_expand(float, nsrc*bcf_hdr_nsamples(opts->hdr_out), opts->mtmpf3, opts->tmpf3);
        for (i=0; i<bcf_hdr_nsamples(opts->hdr_out); i++)
        {
            float *ori = opts->tmpf2 + ndst*i;
            float *dst = opts->tmpf3 + nsrc*i;
            int keep_ori = 0;
            if ( opts->sample_map[i]==-1 ) keep_ori = 1;
            else if ( col->replace==REPLACE_EXISTING && bcf_float_is_missing(ori[0]) ) keep_ori = 1;
            else if ( col->replace==REPLACE_MISSING  && !bcf_float_is_missing(ori[0]) ) keep_ori = 1;
            if ( keep_ori )
            {
                for (j=0; j<ndst; j++) dst[j] = ori[j];
                for (; j<nsrc; j++) bcf_float_set_vector_end(dst[j]);
            }
            else
            {
                float *src = opts->tmpf + nsrc*opts->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_float(opts->hdr_out,line,col->hdr_key,opts->tmpf3,nsrc*bcf_hdr_nsamples(opts->hdr_out));
    }
}
int vcf_setter_format_str(struct vcfs_options *opts, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    opts->tmpp[0] = opts->tmps;
    bcf_hdr_t *header = opts->hdr_out;
    int ret = bcf_get_format_string(header,rec,col->hdr_key,&opts->tmpp,&opts->mtmps);
    opts->tmps = opts->tmpp[0]; // tmps might be realloced
    if ( ret==-3 ) return 0;    // the tag is not present
    if ( ret<=0 ) return 1;     // error

    if ( !opts->sample_map )
        return bcf_update_format_string(opts->hdr_out,line,col->hdr_key,(const char**)opts->tmpp,bcf_hdr_nsamples(opts->hdr_out));

    int i;
    opts->tmpp2[0] = opts->tmps2;
    ret = bcf_get_format_string(opts->hdr_out,line,col->hdr_key,&opts->tmpp2,&opts->mtmps2);
    opts->tmps2 = opts->tmpp2[0];   // tmps2 might be realloced

    if ( ret<=0 )   // not present in dst
    {
        hts_expand(char,bcf_hdr_nsamples(opts->hdr_out)*2,opts->mtmps2,opts->tmps2);
        for (i=0; i<bcf_hdr_nsamples(opts->hdr_out); i++)
        {
            opts->tmps2[2*i]   = '.';
            opts->tmps2[2*i+1] = 0;
            opts->tmpp2[i] = opts->tmps2+2*i;
        }
    }

    for (i=0; i<bcf_hdr_nsamples(opts->hdr_out); i++)
    {
        int isrc = opts->sample_map[i];
        if ( isrc==-1 ) continue;
        opts->tmpp2[i] = opts->tmpp[isrc];
    }
    return bcf_update_format_string(opts->hdr_out,line,col->hdr_key,(const char**)opts->tmpp2,bcf_hdr_nsamples(opts->hdr_out));
}
int vcfs_options_init(struct vcfs_options *opts)
{
    memset(opts, 0, sizeof(struct vcfs_options));
    opts->vcmp = vcmp_init();
    opts->vcfs_is_inited = 1;
    return 0;    
}
int vcf_file_destroy(struct anno_vcf_file *file)
{
    hts_close(file->fp);
    bcf_hdr_destroy(file->hdr);
    if ( file->bcf_idx )
	hts_idx_destroy(file->bcf_idx);
    if ( file->tbx_idx )
	tbx_destroy(file->tbx_idx);
    free(file->fname);
    int i;
    for ( i = 0; i < file->max; ++i ) 
	bcf_destroy1(file->buffer[i]);
    for ( i = 0; i < file->ncols; ++i )
	free(file->cols[i].hdr_key);
    if ( file->ncols )
	free(file->cols);
    return 0;
}
int vcfs_options_destroy(struct vcfs_options *opts)
{
    if ( opts->vcfs_is_inited == 0 )
	return 1;
    int i;
    for ( i = 0; i < opts->n_files; ++i ) 
	vcf_file_destroy(&opts->files[i]);
    free(opts->files);
    vcmp_destroy(opts->vcmp);
    free(opts->tmpks.s);
    free(opts->tmpi);
    free(opts->tmpi2);
    free(opts->tmpi3);
    free(opts->tmpf);
    free(opts->tmpf2);
    free(opts->tmpf3);
    free(opts->tmps);
    free(opts->tmps2);
    free(opts->sample_map);
    return 0;
}
int vcfs_columns_init(struct anno_vcf_file *file, bcf_hdr_t *hdr_out, char *columns)
{
    kstring_t string = KSTRING_INIT;
    kputs(columns, &string);
    int nfields = 0;
    int *splits = ksplit(&string, ',', &nfields);
    if (nfields == 0)
	return 1;
    file->cols = (struct anno_col*)malloc(nfields * sizeof(struct anno_col));
    int i;
    kstring_t temp = KSTRING_INIT;
    file->ncols = 0;
    for ( i = 0; i < nfields; ++i ) {
	struct anno_col *col = &file->cols[file->ncols];
	char *ss = string.s + splits[i];
	col->replace = REPLACE_MISSING;
	if ( *ss == '+') {
	    col->replace = REPLACE_MISSING;
	    ss++;
	} else if ( *ss == '-' ) {
	    col->replace = REPLACE_EXISTING;
	    ss++;
	}
	if ( ss == NULL || *ss == '\0' || *ss == ' ' ) {
	    warnings("Empty tag.");
	    continue;
	}
	if ( strcasecmp("CHROM", ss) == 0 || strcasecmp("POS", ss) == 0 || strcasecmp("FEF", ss) == 0 ||
	     strcasecmp("ALT", ss) == 0 || strcasecmp("FILTER", ss) == 0 || strcasecmp("QUAL", ss) == 0 ) {
	    warnings("Skip %s", ss);
	    continue;
	}	
	if ( strcasecmp("INFO", ss) == 0 || strcasecmp("FORMAT", ss) == 0 ) {
	    warnings("Do not support annotate all INFO or FORMAT fields. todo INFO/TAGS instead.");
	    continue;
	}
	if ( strcasecmp("ID", ss) == 0 ) {	    
	    col->setter.vcf = vcf_setter_id;	    
	} else if ( strncasecmp("FORMAT/", ss, 7) == 0 || strncasecmp("FMT/", ss, 4) == 0) {
	    ss += (strncasecmp("FORMAT/", ss, 7) == 0 ? 7 : 4);
	    if ( strcasecmp("GT", ss) == 0 ) {
		warnings("It is not allowed to change GT tag.");
		continue;
	    }
	    int hdr_id = bcf_hdr_id2int(hdr_out, BCF_DT_ID, ss);
	    if ( !bcf_hdr_idinfo_exists(hdr_out, BCF_HL_FMT, hdr_id) ) {
		bcf_hrec_t *hrec = bcf_hdr_get_hrec(file->hdr, BCF_HL_FMT, "ID", ss, NULL);
		if ( hrec == NULL ) {
		    warnings("The tag \"%s\" is not defined in header. Skip..", ss);
		    continue;
		}
		temp.l = 0;
		bcf_hrec_format(hrec, &temp);
		bcf_hdr_append(hdr_out, temp.s);
		bcf_hdr_sync(hdr_out);
		hdr_id = bcf_hdr_id2int(hdr_out, BCF_DT_ID, ss);
		assert ( bcf_hdr_idinfo_exists(hdr_out, BCF_HL_FMT, hdr_id) );
	    }
	    switch ( bcf_hdr_id2type(hdr_out, BCF_HL_FMT, hdr_id) ) {		
		case BCF_HT_INT :
		    col->setter.vcf = vcf_setter_format_int;
		    break;
		case BCF_HT_REAL:
		    col->setter.vcf = vcf_setter_format_real;
		    break;
		case BCF_HT_STR:
		    col->setter.vcf = vcf_setter_format_str;
		    break;
		default:
		    error("The type of %s not recongized (%d).", ss, bcf_hdr_id2type(hdr_out, BCF_HL_FMT, hdr_id));
	    }
	    col->number = bcf_hdr_id2length(hdr_out, BCF_HL_FMT, hdr_id);
	    // col->icol = hdr_id;
	} else {
	    if ( strncasecmp("INFO/", ss, 5) == 0 )
		ss += 5;
	    int hdr_id = bcf_hdr_id2int(hdr_out, BCF_DT_ID, ss);
	    if ( !bcf_hdr_idinfo_exists(hdr_out, BCF_HL_INFO, hdr_id) ) {
		bcf_hrec_t *hrec = bcf_hdr_get_hrec(file->hdr, BCF_HL_INFO, "ID", ss, NULL);
		if ( hrec == NULL ) {
		    warnings("The tag \"%s\" is not defined in header. Skip.", ss);
		    continue;
		}		    
		temp.l = 0;
		bcf_hrec_format(hrec, &temp);
		bcf_hdr_append(hdr_out, temp.s);
		bcf_hdr_sync(hdr_out);
		hdr_id = bcf_hdr_id2int(hdr_out, BCF_DT_ID, ss);
		assert( bcf_hdr_idinfo_exists(hdr_out, BCF_HL_INFO, hdr_id));
	    }

	    switch ( bcf_hdr_id2type(hdr_out, BCF_HL_INFO, hdr_id) ) {
		case BCF_HT_FLAG:
		    col->setter.vcf = vcf_setter_info_flag;
		    break;
		case BCF_HT_INT:
		    col->setter.vcf = vcf_setter_info_int;
		    break;

		case BCF_HT_REAL:
		    col->setter.vcf = vcf_setter_info_real;
		    break;

		case BCF_HT_STR:
		    col->setter.vcf = vcf_setter_info_str;
		    break;

		default:
		    error("The type of %s not recongized (%d).", ss, bcf_hdr_id2type(hdr_out, BCF_HL_INFO, hdr_id));
	    }
	    col->number = bcf_hdr_id2length(hdr_out, BCF_HL_INFO, hdr_id);
	    // col->icol = hdr_id;
	}
	col->hdr_key = strdup(ss);
	// set the file id to col ifile
	col->ifile = file->id;
	file->ncols++;	
    }
    free(splits);    
    if ( string.m )
	free(string.s);    
    if ( file->ncols == 0 )
	return 1;
    return 0;    
}
int vcfs_database_add(struct vcfs_options *opts, const char *fname, char *columns)
{
    if ( opts->n_files == opts->m_files ) {
	opts->m_files = opts->m_files + 4;
	opts->files = (struct anno_vcf_file *)realloc(opts->files, opts->m_files*sizeof(struct anno_vcf_file));
    }
    struct anno_vcf_file *file = &opts->files[opts->n_files];
    file->fp = hts_open(fname, "r");
    if ( file->fp == NULL ) {
	warnings("Failed to open %s : %s.", fname, strerror(errno));
	return 1;
    }
    htsFormat type = *hts_get_format(file->fp);
    if ( type.format != vcf && type.format != bcf)
	error("Unsupported input format; only accept pre-indexed BCF/VCF file. %s", fname);

    if ( file->fp->format.compression != bgzf )
	error("This file is not compressed by bgzip. %s", fname);
    BGZF *bgzf = hts_get_bgzfp(file->fp);
    if (bgzf && bgzf_check_EOF(bgzf) == 0) {
	warnings("no BGZF EOF marker; file may be truncated. %s", fname);
    }
    file->bcf_idx = NULL;
    file->tbx_idx = NULL;
    if ( type.format == bcf ) {
	file->bcf_idx = bcf_index_load(fname);
	if (file->bcf_idx == NULL)
	    error("Failed to load bcf index of %s", fname);
    } else {
	file->tbx_idx = tbx_index_load(fname);
	if (file->tbx_idx == NULL)
	    error("Failed to load tabix index of %s", fname);
    }
    file->id= opts->n_files;
    file->fname = strdup(fname);
    file->hdr = bcf_hdr_read(file->fp);
    file->itr = NULL;
    file->buffer = NULL;
    file->cached = file->max = 0;
    if ( vcfs_columns_init(file, opts->hdr_out, columns) == 1) {
	bcf_hdr_destroy(file->hdr);
	hts_close(file->fp);
	if ( file->tbx_idx )
	    tbx_destroy(file->tbx_idx);
	if ( file->bcf_idx )
	    hts_idx_destroy(file->bcf_idx);
	return 1;
    }

    opts->n_files++;
    if (opts->vcfs_is_inited == 0)
	opts->vcfs_is_inited = 1;
    return 0;
}
// find the bcf line from database of same positions and allele
static int vcf_fill_buffer(struct anno_vcf_file *file, bcf_hdr_t *hdr_out, bcf1_t *line)
{
    if ( file->itr ) {
	hts_itr_destroy(file->itr);
	file->itr = NULL;
    }
    
    if ( file->tbx_idx ) {
	int tid = tbx_name2id(file->tbx_idx, bcf_seqname(hdr_out, line));
	if ( tid == -1) {
	    warnings("no chromosome %s found in databases %s.", bcf_seqname(hdr_out, line), file->fname);
	    return 1;
	}
	file->itr = tbx_itr_queryi(file->tbx_idx, tid, line->pos, line->pos + line->rlen);
    } else if ( file->bcf_idx ) {
	int tid = bcf_hdr_name2id(file->hdr, bcf_seqname(hdr_out, line));
	if ( tid == -1) {
	    warnings("no chromosome %s found in databases %s.", bcf_seqname(hdr_out, line), file->fname);
	    return 1;
	}
	file->itr = bcf_itr_queryi(file->bcf_idx, tid, line->pos, line->pos + line->rlen);	
    } else {
	error("no index");
    }
    if ( file->itr == NULL )
	return 1;
    file->cached = 0;
    // debug_print("tbx : %p, bcf : %p, cached : %d, max : %d", file->tbx_idx, file->bcf_idx, file->cached, file->max);
    while (1) {
	if ( file->cached == file->max ) {
	    file->max += 8;
	    file->buffer = (bcf1_t**)realloc(file->buffer, sizeof(bcf1_t)*file->max);
	    int i;
	    for (i = 8; i > 0; i--) {
		file->buffer[file->max-i] = bcf_init();
		//file->buffer[file->max-i]->max_unpack = file->max_unpack;
		file->buffer[file->max-i]->pos = -1;
	    } 
	}
	if ( file->tbx_idx ) {
	    kstring_t tmps = KSTRING_INIT;		
	    if ( tbx_itr_next(file->fp, file->tbx_idx, file->itr, &tmps) < 0) 
		break;
	    //debug_print("%s", tmps.s);
	    vcf_parse1(&tmps, file->hdr, file->buffer[file->cached]);
	    free(tmps.s);	    
	} else if ( file->bcf_idx ) {
	    if ( bcf_itr_next(file->fp, file->itr, file->buffer[file->cached]) < 0)
		break;
	} else {
	    error("no index");
	}
	file->cached++;
    } 
    return file->cached ? 0 : 1;      
}
   
bcf1_t *anno_vcfs_core(struct vcfs_options *opts, bcf1_t *line)
{
    // just skip the next steps if vcfs databases not inited
    if ( opts->vcfs_is_inited == 0 )
	return line;
    assert(opts->hdr_out);    
    int i, j;
    int k;    
    // int ret = 0;
    // debug_print("n_files : %d", opts->n_files);
    for ( i = 0; i < opts->n_files; ++i ) {
	struct anno_vcf_file *file = &opts->files[i];
	vcf_fill_buffer(file, opts->hdr_out, line);
	// debug_print("%d, file: %d", file->cached, i);
	for ( j = 0; j < file->cached; ++j ) {
	    bcf1_t *dat = file->buffer[j];
	    // todo : for deletion, rough search
	    if ( dat->rlen != line->rlen)
		continue;
	    if (dat->pos != line->pos)
		continue;
	    if (dat->pos > line->pos)
		break;
	    // debug_print("dat->pos : %d, line->pos : %d", dat->pos,line->pos);
	    for ( k = 0; k < file->ncols; ++k ) {
		struct anno_col *col = &file->cols[k];
		// col->replace = REPLACE_MISSING;
		// debug_print("ifile : %d, hdr_key : %s", col->ifile, col->hdr_key);
		col->setter.vcf(opts, line, col, dat);
	    }
	}
    }
    return line;
}
#ifdef _VCF_ANNOS_MAIN

#include "config.h"
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

const char *json_fname = 0;
const char *input_fname = 0;
const char *output_fname_type = 0;
const char *output_fname = 0;

int main(int argc, char **argv)
{
    if (argc == 1)
	error("Usage : vcf_annos -c config.json -O z -o output.vcf.gz input.vcf.gz");
    int i;
    for ( i = 1; i < argc; ) {
	const char *a = argv[i++];
	const char **var = 0;
	if ( strcmp(a, "-c") == 0 )
	    var = &json_fname;
	else if ( strcmp(a, "-O") == 0 )
	    var = &output_fname_type;
	else if ( strcmp(a, "-o") == 0 )
	    var = &output_fname;

	if ( var != 0 ) {
	    if ( i == argc )
		error("Missing an argument after %s", a);
	    *var = argv[i++];
	    continue;
	}

	if ( input_fname == 0 ) {
	    input_fname = a;
	    continue;
	}

	error("Unknown argument : %s.", a);
    }

    struct vcfanno_config *con = vcfanno_config_init();
    if ( vcfanno_load_config(con, json_fname) != 0 )
	error("Failed to load configure file. %s : %s", json_fname, strerror(errno));
    vcfanno_config_debug(con);
    if ( con->vcfs.n_vcfs == 0)
	error("No BCF/VCF database specified.");
    if (input_fname == 0 && (!isatty(fileno(stdin))))
	input_fname = "-";
    if (input_fname == 0)
	error("No input file.");

    int out_type = FT_VCF;
    if (output_fname_type != 0) {
	switch (output_fname_type[0]) {
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

    htsFile *fp = hts_open(input_fname, "r");
    if ( fp == NULL )
	error("Failed to open %s : %s.", input_fname, strerror(errno));
    htsFormat type = *hts_get_format(fp);
    if ( type.format != vcf && type.format != bcf )
	error("Unsupported input format. %s", input_fname);
    
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if ( hdr == NULL )
	error("Failed to parse header.");	
    bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);    
    htsFile *fout = output_fname == 0 ? hts_open("-", hts_bcf_wmode(out_type)) : hts_open(output_fname, hts_bcf_wmode(out_type));
    struct vcfs_options opts = { .vcfs_is_inited = 0,};
    vcfs_options_init(&opts);
    opts.hdr_out = hdr_out;
    
    for ( i = 0; i < con->vcfs.n_vcfs; ++i ) {
	vcfs_database_add(&opts, con->vcfs.files[i].fname, con->vcfs.files[i].columns);
    }

    bcf_hdr_write(fout, hdr_out);
    bcf1_t *line = bcf_init();
    while ( bcf_read(fp, hdr, line) == 0 ) {
	anno_vcfs_core(&opts, line);
	bcf_write1(fout, hdr_out, line);
    }
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    bcf_hdr_destroy(hdr_out);
    vcfs_options_destroy(&opts);
    hts_close(fp);
    hts_close(fout);
    return 0;
}

#endif
