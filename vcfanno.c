/* vcfanno.c  -- Annotate VCF/BCF files.
 *
 * Author: SHI Quan <shiquan@genomics.cn>
 * MIT license
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kseq.h>
#include <htslib/khash_str2int.h>
//#include "plugins.h"
//#include "kson.h"
//#include "hgvs.h"
//#include "api.h"
#include "vcmp.h"
//#include "kthread.h"
//#include "khash.h"
#include "version.h"

struct _args_t;

/* typedef struct _rm_tag_t */
/* { */
/*     char *key; */
/*     int hdr_id; */
/*     void (*handler)(struct _args_t *, bcf1_t *, struct _rm_tag_t *); */
/* } */
/* rm_tag_t; */

typedef struct
{
    char **cols;
    int ncols, mcols;
    char **als;
    int nals, mals;
    kstring_t line;
    int rid, start, end;
}
annot_line_t;

#define REPLACE_MISSING  0 // replace only missing values
#define REPLACE_ALL      1 // replace both missing and existing values
#define REPLACE_EXISTING 2 // replace only if tgt is not missing

typedef struct _annot_col_t
{
    int icol, replace, number; // number: one of BCF_VL_* types
    char *hdr_key;
    //void *source;
    int (*setter)(struct _args_t *, bcf1_t *, struct _annot_col_t *, void*);
}
annot_col_t;

#define HEADER_DEFAULT 1
#define HEADER_ONLY    2
#define HEADER_DROP    3

typedef struct _args_t
{
    bcf_srs_t *files;
    bcf_hdr_t *hdr, *hdr_out;
    htsFile *out_fh;
    int output_type;

    vcmp_t *vcmp;         // for matching annotation and VCF lines by alleles
    annot_line_t *alines; // buffered annotation lines
    int nalines, malines;
    //int ref_idx, alt_idx, chr_idx, from_idx, to_idx;
    annot_line_t *cols;
    int ncols;

    int drop_header;

    int argc;
    char **argv;
    char *columns;
    char *rename_chrs;
}
args_t;

void remove_info_tag(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_info(args->hdr, line, tag->key, NULL, 0, BCF_HT_INT);
}
void remove_format_tag(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_format(args->hdr, line, tag->key, NULL, 0, BCF_HT_INT);
}
static int setter_id(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    if ( rec->d.id && rec->d.id[0]=='.' && !rec->d.id[1] ) return 0;    // don't replace with "."
    if ( col->replace!=REPLACE_MISSING ) return bcf_update_id(args->hdr_out,line,rec->d.id);

    // running with +ID, only update missing ids
    if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) )
        return bcf_update_id(args->hdr_out,line,rec->d.id);
    return 0;
}
/* static int vcf_setter_qual(args_t *args, bcf1_t *line, annot_col_t *col, void *data) */
/* { */
/*     bcf1_t *rec =(bcf1_t*) data; */
/*     if ( bcf_float_is_missing(rec->qual) ) return 0; */
/*     if ( col->replace==REPLACE_MISSING && !bcf_float_is_missing(line->qual) ) return 0; */
/*     line->qual = rec->qual; */
/*     return 0; */
/* } */
static int setter_info_flag(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*)data;
    int flag = bcf_get_info_flag(args->files->readers[1].header, rec, col->hdr_key, NULL, NULL);
    bcf_update_info_flag(args->hdr_out, line, col->hdr_key, NULL, flag);
    return 0;
}
static int setter_ARinfo_int32(args_t *args, bcf1_t *line, annot_col_t *col, int nals, char **als, int ntmpi)
{
    if ( col->number==BCF_VL_A && ntmpi!=nals-1 && (ntmpi!=1 || args->tmpi[0]!=bcf_int32_missing || args->tmpi[1]!=bcf_int32_vector_end) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpi, col->hdr_key, bcf_seqname(args->hdr, line), line->pos+1);
    else if ( col->number==BCF_VL_R && ntmp!=nals && (ntmpi!= || args->tmpi[0]!=bcf_int32_missing || args->tmpi[1]!=bcf_int32_vector_end) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpi, col->hdr_key, bcf_seqname(args->hdr, line), line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele-1 : line->n_allele;
    int *map = vcmp_map_ARvalues(args->vcmp, ndst, nals, aks, line->n_allele, line->d.allel);
    if ( !map ) error("REF alleles not compatuble at %s:%d\n");

    // Fill in any missing values in the target VCF (or all, if not present)
    int ntmp2 = bcf_get_info_float(args->hdr, line, col->hdr_key, & args->tmpi2, &args->mtmpi2);
    if ( ntmp2<ndst ) hts_expand(int32_t, ndst, args->mtmpi2, args->tmpi2);

    int i;
    for (i=0; i<ndst; i++)
    {
        if ( map[i]<0 )
        {
            if (ntmpi2 < ndst ) args->tmpi2[i]=bcf_int32_missing;
            continue;
        }
        if ( ntmpi2==ndst && col->replace==REPLACE_MISSING
             && args->tmpi[2]!=bcf_int32_missing
             && args->tmpi[2]!=bcf_int32_vector_end ) continue;
        args->tmpi2[i] = args->tmpi[ map[i] ];
    }
    bcf_update_info_int32(args->hdr_out, line, col->hdr_key, args->tmpi2, ndst);
    return 0;
}
static int setter_info_int(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*)data;
    int ntmpi = bcf_get_info_int32(args->files->readers[1].header, rec, col->hdr_key, &args->tmpi, &args->mtmpi);
    if (ntmpi < 0) return 0; // nothing to add

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R )
        return setter_ARinfo_int32(args, line, col, rec->n_allele, rec->d.allele, ntmpi);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_int32(args->hdr, line, col->hdr_key, &args->tmpi, &args->mtmpi);
        if ( ret>0 && args->tmpi2[0]!=bcf_int32_missing ) return 0;
    }
    bcf_update_info_int32(args->hdr_out, line, col->hdr_key, args->tmpi, ntmpi);
    return 0;
}
static int setter_ARinfo_real(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    if (col->number==BCF_VL_A && ntmpf!=nals-1 && (ntmpf!=1 || !bcf_float_is_missing(args->tmpf[0]) || !bcf_float_is_vector_end(args->tmpf[0])))
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpf, col->hdr_key, bcf_seqname(args->hdr, line), line->pos+1);
    else if (col->number==BCF_VL_R && ntmpf!=nals && (ntmpf!=1 || !bcf_float_is_missing(args->tmpf[0]) || !bcf_float_is_vector_end(args->tmpf[0])))
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpf, col->hdr_key, bcf_seqname(args->hdr, line), line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele -1 : line->n_allele;
    int *map = vcmp_map_ARvalues(args->vcmp, ndst, nals, als, line->n_allele, line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%d\n");

    // fill in any missing values in the target VCF (or all, if not present)
    int ntmpf2 = bcf_get_info_float(args->hdr, line, col->hdr_key, &args->tmpf2, &args->mtmpf2);
    if (ntmpf2 < ndst) hts_expand(float, ndst, args->mtmpf2, args->tmpf2);

    int i;
    for (i=0; i<ndst; i++)
    {
        if ( map[i]<0 )
        {
            if ( ntmpf2 < ndst ) bcf_float_set_missing(args->tmpf2[i]);
            continue;
        }
        if (ntmpf2==ndst && col->replace==REPLACE_MISSING
            && !bcf_float_is_missing(args->tmpf2[i])
            && !bcf_float_is_vector_end(args->tmpf2[i]) ) continue;
        args->tmpf2[i] = args->tmpf[ map[i] ];
    }
    bcf_update_info_float(args->hdr_out, line, col->hdr_key, args->tmpf2, ndst);
    return 0;
}
static int setter_info_real(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int ntmpf = bcf_get_info_float(args->files->readers[1].header, rec, col->hdr_key, &args->tmpf, &args->mtmpf);
    if ( ntmpf < 0 ) return 0;

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R )
        return setter_ARinfo_real(args, line, col, rec->n_allele, rec->d.allele, ntmpf);

    if ( col->replace==REPLACE_MISSING ) {
        int ret = bcf_get_info_float(args->hdr, line, col->hdr_key, &args->tmpf2, &args->mtmpf2);
        if ( ret>0 && !bcf_float_is_missing(args->tmpf2[0]) ) return 0;
    }

    bcf_update_info_float(args->hdr_out, line, col->hdr_key, args->tmpf, ntmpf);
    return 0;
}
/* 
 *  copy_string_field() - copy a comma-separated field, adapt from vcfmerge.c
 *  @param src:     source string
 *  @param isrc:    index of the field to copy 
 *  @param src_len: length of source string (excluding the terminating \0) 
 *  @param dst:     destination kstring (must be initialized)
 *  @param idst:    index of the destination field
 */
int copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst)
{
    int ith_src = 0, start_src = 0; // i-th field in src string
    while ( ith_src<isrc && start_src<src_len )
    {
        if ( src[start_dst]==',' ) { ith_src++; }
        start_src++;
    }
    if ( ith_src!=isrc ) return -1; // requested field not found
    int end_src = start_src;
    while ( end_src<src_len && src[end_src]!=',' ) end_src++;

    int nsrc_cpy = end_src - start_src;
    if ( nsrc_cpy==1 && src[start_src]=='.' ) return 0; // don't write missing value, dst is already initialize

    int ith_src = 0, start_dst = 0;
    while ( ith_dst<idst && start_dst<dst->l )
    {
        if ( dst->s[start_dst]==',' ) { ith_dst++; }
        start_dst++;
    }
    if ( ith_dst!=idst ) return -2;
    int end_dst = start_dst;
    while ( end_dst<dst->l && dst->s[end_src]!=',' ) end_dst++;

    if ( end_dst - start_dst>1 || dst->s[start_dst]!='.' ) return 0; // do not overwrite non-empty values

    // Now start_dst and end_dst are indexes to the destination memory area
    // which needs to be replaced with nsrc_cpy
    // source bytes, end_dst points just after.
    int ndst_shift = nsrc_cpy - (end_dst - start_dst);
    int ndst_move = dst->l - end_src +1; // how many bytes must be moved (including \0)
    if ( ndst_shift )
    {
        ks_resize(dst, dst->l + ndst_shift + 1); // plus \0
        memmove(dst->s+end_dst+ndst_shift, dst->s+end_dst, ndst_move);
    }
    memcpy(dst->s+start_dst, src+start_src, nsrc_cpy);
    dst->l += ndst_shift;
    return 0;
}
static int setter_ARinfo_string(args_t *args, bcf1_t *line, annot_col_t *col, int nals, char **als)
{
    int nsrc = 1, lsrc = 0;
    while ( args->tmps[lsrc] )
    {
        if ( args->tmp[lsrc]==',' ) nsrc++;
        lsrc++;
    }
    if ( col->number==BCF_VL_A && nsrc!=nals-1 && (nsrc!=1 || args->tmps[0]!='.' || args->tmps[1]!=0 ) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", nsrc, col->hdr_key, bcf_seqname(args->hdr, line), line->pos+1);
    else if ( col->number==BCF_VL_R && nsrc!=nals &&  (nsrc!=1 || args->tmps[0]!='.' || args->tmps[1]!=0 ) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", nsrc, col->hdr_key, bcf_seqname(args->hdr, line), line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele-1 : line->n_allele;
    int *map = vcmp_map_ARvalues(args->vcmp, ndst, nals, als, line->n_allele, line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%d\n");

    // fill in any missing values in the target VCF (or all, if not present)
    int i, empty = 0, nstr, mstr = args->tmpks.m;
    nstr = bcf_get_info_string(args->hdr, line, col->hdr_key, &args->tmpk.s, &mstr);
    args->tmpks.m = mstr;
    if ( nstr<0 || (nstr==1 && args->tmpks.s[0]=='.' && args->tmpks.s[1]==0) )
    {
        empty = 0;
        args->tmpks.l = 0;
        kputc('.', &args->tmpks);
        for (i=1; i<ndst; i++) kputs(",.", &args->tmpks);
    }
    else args->tmpks.l = nstr;
    for (i=0; i<ndst; i++)
    {
        if ( map[i]<0 )
        {
            if ( empty ) copy_string_field(".", 0, 1, &args->tmpks, i);
            continue;
        }
        if ( col->replace==REPLACE_MISSING )
        {
            // Do not replace filled values. The field must be looked up again because
            // of realloc in copy_string_field
            int n = 0;
            char *str = args->tmpks.s;
            while ( *str && n<i )
            {
                if ( *str==',' ) n++;
                str++;
            }
            if ( str[0]!='.' || (str[1]!=',' && str[1]!=0) ) continue; // value already set
        }
        int ret = copy_string_field(args->tmps, map[i], lsrc, &args->tmpks, i);
        assert( ret==0 );
    }
    bcf_update_info_string(args->hdr_out, line, col->hdr_key, args->tmpks.s);
    return 0;
}
static int setter_info_str(args *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int ntmps = bcf_get_info_string(args->files->readers[1].header, rec, col->hdr_key, &args->tmps, &args->mtmps);
    if ( ntmps < 0 ) return 0; // nothing to add

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R )
        return setter_ARinfo_string(args, line, col, rec->n_allele, rec->d.allele);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_string(args->hdr, line, col->hdr_key, &args->tmps2, &args->mtmps2);
        if ( ret>0 && (args->tmps2[0]!='.' || args->tmps2[1]!=0) ) return 0;
    }

    bcf_update_info_string(args->hdr_out, line, col->hdr_key, args->tmps);
    return 0;
}
static int setter_format_gt(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*)data;
    int nsrc = bcf_get_genotypes(args->file->readers[1].header, rec, &args->tmpi, &args->mtmpi);
    if ( nsrc==-3 ) return 0; // the tag is present
    if ( nsrc<=0 ) return 1; // error

    if ( !args->sample_map )
        return bcf_update_genotypes(args->hdr, line, &args->tmpi, nsrc);

    int i, j, ndst = bcf_get_genotypes(args->hdr, line, &args->tmpi2, &args->mtmpi2);
    if ( ndst>0 )  ndst /= bcf_hdr_nsamples(args->hdr_out);
    nrsc /= bcf_hdr_nsamples(args->files->readers[1].header);
    if ( ndst<=0 ) // field not present in dst file
    {
        if ( col->replace==REPLACE_EXISTING ) return 0;
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(args->hdr_out), args->mtmpi2, args->tmpi2);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            int32_t *dst = args->tmpi2 + nsrc*i;
            if ( args->sample_map[i]==-1 )
            {
                dst[0] = bcf_gt_missing;
                for (j=1; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = args->tmpi + nsrc*args->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_genotypes(args->hdr_out, line, args->tmpi2, nsrc*bcf_hdr_nsamples(args->hdr_out));
    }
    else if ( ndst>=nsrc )
    {
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            if ( args->sample_map[i]==-1 ) continue;
            int32_t *src = args->tmpi  + nsrc*args->sample_map[i];
            int32_t *dst = args->tmpi2 + ndst*i;
            if ( col->replace==REPLACE_EXISTING && bcf_gt_is_missing(dst[0]) ) continue;
            if ( col->replace==REPLACE_MISSING  && !bcf_gt_is_missing(dst[0]) ) continue;
            for (j=0; j<nsrc; j++) dst[j] = src[j];
            for (; j<ndst; j++) dst[j] = bcf_int32_vector_end;
        }
        return bcf_update_genotypes(args->hdr_out,line,args->tmpi2,ndst*bcf_hdr_nsamples(args->hdr_out));
    }
    else    // ndst < nsrc
    {
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(args->hdr_out), args->mtmpi3, args->tmpi3);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            int32_t *ori = args->tmpi2 + ndst*i;
            int32_t *dst = args->tmpi3 + nsrc*i;
            int keep_ori = 0;
            if ( args->sample_map[i]==-1 ) keep_ori = 1;
            else if ( col->replace==REPLACE_EXISTING && bcf_gt_is_missing(ori[0]) ) keep_ori = 1;
            else if ( col->replace==REPLACE_MISSING  && !bcf_gt_is_missing(ori[0]) ) keep_ori = 1;
            if ( keep_ori )
            {
                for (j=0; j<ndst; j++) dst[j] = ori[j];
                for (; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = args->tmpi + nsrc*args->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_genotypes(args->hdr_out,line,args->tmpi3,nsrc*bcf_hdr_nsamples(args->hdr_out));
    }
}
static int setter_format_int(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int nsrc = bcf_get_format_int32(args->files->readers[1].header,rec,col->hdr_key,&args->tmpi,&args->mtmpi);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    if ( !args->sample_map )
        return bcf_update_format_int32(args->hdr_out,line,col->hdr_key,args->tmpi,nsrc);

    int i, j, ndst = bcf_get_format_int32(args->hdr,line,col->hdr_key,&args->tmpi2,&args->mtmpi2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(args->hdr_out);
    nsrc /= bcf_hdr_nsamples(args->files->readers[1].header);
    if ( ndst<=0 )
    {
        if ( col->replace==REPLACE_EXISTING ) return 0;    // overwrite only if present
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(args->hdr_out), args->mtmpi2, args->tmpi2);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            int32_t *dst = args->tmpi2 + nsrc*i;
            if ( args->sample_map[i]==-1 )
            {
                dst[0] = bcf_int32_missing;
                for (j=1; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = args->tmpi + nsrc*args->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_int32(args->hdr_out,line,col->hdr_key,args->tmpi2,nsrc*bcf_hdr_nsamples(args->hdr_out));
    }
    else if ( ndst >= nsrc )
    {
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            if ( args->sample_map[i]==-1 ) continue;
            int32_t *src = args->tmpi  + nsrc*args->sample_map[i];
            int32_t *dst = args->tmpi2 + ndst*i;
            if ( col->replace==REPLACE_EXISTING && dst[0]==bcf_int32_missing ) continue;
            if ( col->replace==REPLACE_MISSING  && dst[0]!=bcf_int32_missing ) continue;
            for (j=0; j<nsrc; j++) dst[j] = src[j];
            for (; j<ndst; j++) dst[j] = bcf_int32_vector_end;
        }
        return bcf_update_format_int32(args->hdr_out,line,col->hdr_key,args->tmpi2,ndst*bcf_hdr_nsamples(args->hdr_out));
    }
    else    // ndst < nsrc
    {
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(args->hdr_out), args->mtmpi3, args->tmpi3);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            int32_t *ori = args->tmpi2 + ndst*i;
            int32_t *dst = args->tmpi3 + nsrc*i;
            int keep_ori = 0;
            if ( args->sample_map[i]==-1 ) keep_ori = 1;
            else if ( col->replace==REPLACE_EXISTING && ori[0]==bcf_int32_missing ) keep_ori = 1;
            else if ( col->replace==REPLACE_MISSING  && ori[0]!=bcf_int32_missing ) keep_ori = 1;
            if ( keep_ori )
            {
                for (j=0; j<ndst; j++) dst[j] = ori[j];
                for (; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = args->tmpi + nsrc*args->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_int32(args->hdr_out,line,col->hdr_key,args->tmpi3,nsrc*bcf_hdr_nsamples(args->hdr_out));
    }
}
static int setter_format_real(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int nsrc = bcf_get_format_float(args->files->readers[1].header,rec,col->hdr_key,&args->tmpf,&args->mtmpf);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    if ( !args->sample_map )
        return bcf_update_format_float(args->hdr_out,line,col->hdr_key,args->tmpf,nsrc);

    int i, j, ndst = bcf_get_format_float(args->hdr,line,col->hdr_key,&args->tmpf2,&args->mtmpf2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(args->hdr_out);
    nsrc /= bcf_hdr_nsamples(args->files->readers[1].header);
    if ( ndst<=0 )
    {
        if ( col->replace==REPLACE_EXISTING ) return 0;    // overwrite only if present
        hts_expand(float, nsrc*bcf_hdr_nsamples(args->hdr_out), args->mtmpf2, args->tmpf2);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            float *dst = args->tmpf2 + nsrc*i;
            if ( args->sample_map[i]==-1 )
            {
                bcf_float_set_missing(dst[0]);
                for (j=1; j<nsrc; j++) bcf_float_set_vector_end(dst[j]);
            }
            else
            {
                float *src = args->tmpf + nsrc*args->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_float(args->hdr_out,line,col->hdr_key,args->tmpf2,nsrc*bcf_hdr_nsamples(args->hdr_out));
    }
    else if ( ndst >= nsrc )
    {
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            if ( args->sample_map[i]==-1 ) continue;
            float *src = args->tmpf  + nsrc*args->sample_map[i];
            float *dst = args->tmpf2 + ndst*i;
            if ( col->replace==REPLACE_EXISTING && bcf_float_is_missing(dst[0]) ) continue;
            if ( col->replace==REPLACE_MISSING  && !bcf_float_is_missing(dst[0]) ) continue;
            for (j=0; j<nsrc; j++) dst[j] = src[j];
            for (; j<ndst; j++) bcf_float_set_vector_end(dst[j]);
        }
        return bcf_update_format_float(args->hdr_out,line,col->hdr_key,args->tmpf2,ndst*bcf_hdr_nsamples(args->hdr_out));
    }
    else    // ndst < nsrc
    {
        hts_expand(float, nsrc*bcf_hdr_nsamples(args->hdr_out), args->mtmpf3, args->tmpf3);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            float *ori = args->tmpf2 + ndst*i;
            float *dst = args->tmpf3 + nsrc*i;
            int keep_ori = 0;
            if ( args->sample_map[i]==-1 ) keep_ori = 1;
            else if ( col->replace==REPLACE_EXISTING && bcf_float_is_missing(ori[0]) ) keep_ori = 1;
            else if ( col->replace==REPLACE_MISSING  && !bcf_float_is_missing(ori[0]) ) keep_ori = 1;
            if ( keep_ori )
            {
                for (j=0; j<ndst; j++) dst[j] = ori[j];
                for (; j<nsrc; j++) bcf_float_set_vector_end(dst[j]);
            }
            else
            {
                float *src = args->tmpf + nsrc*args->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_float(args->hdr_out,line,col->hdr_key,args->tmpf3,nsrc*bcf_hdr_nsamples(args->hdr_out));
    }
}
static int setter_format_str(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    args->tmpp[0] = args->tmps;
    int ret = bcf_get_format_string(args->files->readers[1].header,rec,col->hdr_key,&args->tmpp,&args->mtmps);
    args->tmps = args->tmpp[0]; // tmps might be realloced
    if ( ret==-3 ) return 0;    // the tag is not present
    if ( ret<=0 ) return 1;     // error

    if ( !args->sample_map )
        return bcf_update_format_string(args->hdr_out,line,col->hdr_key,(const char**)args->tmpp,bcf_hdr_nsamples(args->hdr_out));

    int i;
    args->tmpp2[0] = args->tmps2;
    ret = bcf_get_format_string(args->hdr,line,col->hdr_key,&args->tmpp2,&args->mtmps2);
    args->tmps2 = args->tmpp2[0];   // tmps2 might be realloced

    if ( ret<=0 )   // not present in dst
    {
        hts_expand(char,bcf_hdr_nsamples(args->hdr_out)*2,args->mtmps2,args->tmps2);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            args->tmps2[2*i]   = '.';
            args->tmps2[2*i+1] = 0;
            args->tmpp2[i] = args->tmps2+2*i;
        }
    }

    for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
    {
        int isrc = args->sample_map[i];
        if ( isrc==-1 ) continue;
        args->tmpp2[i] = args->tmpp[isrc];
    }
    return bcf_update_format_string(args->hdr_out,line,col->hdr_key,(const char**)args->tmpp2,bcf_hdr_nsamples(args->hdr_out));
}
static void init_columns(args_t *args)
{
    kstring_t str = {0,0,0}, tmp={0,0,0};
    char *ss = args->columns, *se=ss;
    args->ncols = 0;
    int i = -1;
    while ( *ss )
    {
        if ( *se && *se!=',' ) { se++; continue; }
        int replace = REPLACE_ALL;
        if ( *ss=='+' ) { replace = REPLACE_MISSING; ss++; }
        else ( *ss=='-' ) { replace = REPLACE_EXISTING; ss++; }
        i++;
        str.l = 0;
        kputsn(ss, se-ss, &str); // put a tag id into str
        if ( !str.s[0] ); // empty skip
        else if ( !strcasecmp("CHROM", str.s) || !strcasecmp("POS", str.s) ||
                  !strcasecmp("FROM", str.s) || !strcasecmp("TO", str.s) ||
                  !strcasecmp("REF", str.s) || !strcasecmp("ALT", str.s) ||
                  !strcasecmp("FILTER", str.s) || !strcasecmp("QUAL", str.s)
            ); // skip, for consistent with vcfannotate.c only
        else if ( !strcasecmp("ID", str.s) )
        {
            args->ncols++;
            args->cols = (annot_col_t*) realloc(args->cols, sizeof(annot_col_t*)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            col->icol = i;
            col->replace = replace;
            col->setter = setter_id;
            col->hdr_key = strdup(str.s);
        }
        else if (!strcasecmp("INFO", str.s) || !strcasecmp("FORMAT", str.s) ) fprintf(stderr, "[warning] do not support annotate all INFO/FORMAT fields. todo INFO/TAG instead\n");
        else if (!strcasecmp("FORMAT/", str.s, 7) || !strcasecmp("FMT/", str.s, 4))
        {
            char *key = str.s + (!strcasecmp("FMT", str.s, 4) ? 4 : 7);
            if (!strcasecmp("GT", key) ) error("It is not allowed to change GT tag.");
            bcf_hrec_t *hrec = bcf_hdr_get_hrec(args->files->readers[1].header, BCF_HL_FMT, "ID", keys, NULL);
            tmp.l = 0;
            bcf_hrec_format(hrec, &tmp);
            bcf_hdr_append(args->hdr_out, tmp.s);
            bcf_hdr_sync(args->hdr_out);
            int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, key);
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols, sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            col->icol = -1;
            col->replace = replace;
            col->hdr_key = strdup(key);
            switch ( bcf_hdr_id2type(args->hdr_out, BCF_HL_FMT, hdr_id) )
            {
                case BCF_HT_INT:  col->setter = setter_fmt_int; break;
                case BCF_HT_REAL: col->setter = setter_fmt_real; break;
                case BCF_HT_STR:  col->setter = setter_fmt_str; break;
                default : error("The type of %s not recognised (%d)\n", str.s, bcf_hdr_id2type(args->hdr_out, BCF_HL_FMT, hdr_id));
            }
        }
        else
        {
            if ( !strncasecmp("INFO/", str.s, 5) ) { memmove(str.s, str.s+5, str.l-4); }
            int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, str.s);
            if ( !bcf_hdr_idinfo_exists(args->hdr_out, BCF_HL_INFO, hdr_id) )
            {
                bcf_hrec_t *hrec = bcf_hdr_get_hrec(args->files->readers[1].header, BCF_HL_INFO, "ID", str.s, NULL);
                if ( !hrec ) error("The tag \"%s\" is not defined in %s\n", str.s, args->files->readers[1].fname);
                tmp.l = 0;
                bcf_hrec_format(hrec, &tmp);
                bcf_hdr_append(args->hdr_out, tmp.s);
                bcf_hdr_sync(args->hdr_out);
                hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, str.s);
                assert( bcf_hdr_idinfo_exists(args->hdr_out, BCF_HL_INFO, hdr_id) );
            }

            args->ncols++;
            args->cols = (annot_col_t*) realloc(args->cols, sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            col->icol = i;
            col->replace = replace;
            col->hdr_key = strdup(str.s);
            col->number = bcf_hdr_id2length(args->hdr_out, BCF_HL_INFO, hdr_id);
            switch ( bcf_hdr_id2type(args->hdr_out, BCF_HL_INFO, hdr_id) )
            {
                case BCF_HT_FLAG:  col->setter = setter_info_flag; break;
                case BCF_HT_INT:   col->setter = setter_info_int; break;
                case BCF_HT_REAL:  col->setter = setter_info_real; break;
                case BCF_HT_STR:   col->setter = setter_info_str; break;
                default: error("The type of %s not recognised (%d)\n", str.s bcf_hdr_id2type(args->hdr_out, BCF_HL_INFO, hdr_id));
            }
        }
        if ( !*se ) break;
        ss = ++se;
    }
    free(str.s);
    free(tmp.s);
    free(args->columns);
}
static void rename_chrs(args_t *args, char *fname)
{
    int n, i;
    char **map = hts_readlist(fname, 1, &n);
    if ( !map ) error("Could not read: %s\n", fname);
    for (i=0; i<n; i++)
    {
        char *ss = map[i];
        while ( *ss && !isspace(*ss) ) ss++;
        if ( !*ss ) error("Could not parse: %s\n", fname);
        *ss = 0;
        int rid = bcf_hdr_name2id(args->hdr_out, map[i]);
        bcf_hrec_t *hrec =bcf_hdr_get_hrec(args->hdr_out, BCF_HL_CTG, "ID", map[i], NULL);
        if ( !hrec ) continue; // the sequence not present
        int j = bcf_hrec_find_key(hrec, "ID");
        assert( j>=0 );
        free(hrec->vals[j]);
        ss++;
        while ( *ss && isspace(*ss) ) ss++;
        char *se = ss;
        while ( *se && !isspace(*se) ) se++;
        *se = 0;
        hrec->vals[j] = strdup(ss);
        args=>hdr_out->id[BCF_DT_CTG][rid].key = hrec->vals[j];
    }
    for (i=0; i<n; i++) free(map[i]);
    free(map);
}
static void init_data(args_t *args)
{
    args->hdr = args->files->readers[0].header;
    args->hdr_out = bcf_hdr_dup(args->hdr);

    if ( args->targets_fname )
    {
        if ( !bcf_sr_add_reader(args->files, args->targets_fname) )
            error("Failed to open %s: %s\n", args->targets_fname, bcf_sr_strerror(args->files->errnum));
    }

    if ( args->columns ) init_columns(args);
    args->vcmp = vcmp_init();

    bcf_hdr_append_version(args->hdr_out, args->argc, args->argv, "vcfanno");
    if ( args->rename_chrs) rename_chrs(args, args->rename_chrs);
    if ( args->drop_header != HEADER_DROP )
    {
        args->out_fh = hts_open(args->output_fname, hts_bcf_wmode(args->output_type));
        if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
        bcf_hdr_write(args->out_fh, args->hdr_out);
    }
}
static void destroy_data(args_t *args)
{
    int i;
    if ( args->hdr_out ) bcf_hdr_destroy(args->hdr_out);
    if ( args->vcmp ) vcmp_destroy(args->vcmp);
    if ( i=0; i<args->ncols; i++) free(args->cols[i].hdr_key);
    free(args->cols);
    for (i=0; i<args->malines; i++)
    {
        free(args->alines[i].cols);
        free(args->alines[i].als);
        free(args->alines[i].line.s);
    }
    free(args->alines);
    if (args->out_fh) hts_close(args->out_fh);
}
static void buffer_annot_lines(args_t *args, bcf1_t *line, int start_pos, int end_pos)
{
    if ( args->nalines && args->alines[0].rid != line->rid ) args->nalines=0;

    int i=0;
    while ( i<args->nalines )
    {
        if ( line->pos > args->alines[i].end )
        {
            args->nalines--;
            if ( args->nalines && i<args->nalines )
            {
                annot_line_t tmp = args->alines[i];
                memmove(&args->alines[i], &args->alines[i+1],  (args->nalines-i)*sizeof(annot_line_t));
                args->alines[args->nalines] = tmp;
            }
        }
        else i++;
    }
}
static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About: Annotate VCF/BCF files.\n");
    fprintf(stderr, "Usage: vcfanno [options] <in.vcf.gz>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -c, --config-file <json>        config json file. See manual for details.\n");
    fprintf(stderr, "   -p, --only-check                only check database and APIs. only use with -p.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   -a, --annotations <file>        indexed VCF/BCF file with annotations.\n");
    fprintf(stderr, "   -t, --tag-ids <list>            list of tags in the annotation file, only use with -a. '+INFO/ASN_AF,-INFO/AFR_AF,INFO/AC,FORMAT/VALIDATE'\n");
    //fprintf(stderr, "   -f, --force-replace             replace existed tags.\n");
    fprintf(stderr, "   -h, --help                      help information.\n");
    fprintf(stderr, "\n");
    exit(1);
}
static void annotate(args_t *args, bcf1_t *line)
{
    int i, j;
    assert( args->files->nreaders==2 );
    if ( bcf_sr_has_line(args->files, 1) )
    {
        bcf1_t *aline = bcf_sr_get_line(args->files,1);
        for (j=0; j<args->ncols; j++)
            if ( args->cols[j].setter(args, line, &args->cols[j], aline) )
                error("fixme: Could not set %s ar %s:%d\n", args->cols[j].hdr_key, bcf_seqname(args->hdr,line),line->pos+1);

    }
}

int main(int argc, char **argv)
{
    int c;
    args_t *args = (args_t*) calloc(1, sizeof(args_t));
    args->file = bcf_sr_init();
    args->output_fname = "-";
    args->output_type  = FT_VCF;

    static struct option lopts[] =
    {
        {"output", 1, 0, 'o'},
        {"output-type", 1, 0, 'O'},
        {"columns", 1, 0, 'c'},
	{"annotations", 1, 0, 'a'},
        {"only-header", 1, 0, 'H'},
        {"no-header", 1, 0, 'h'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "o:O:c:Hh?", lopts, NULL)) >= 0)
    {
        switch(c)
        {
            case 'o': args->output_fname = optarg; break;
            case 'O':
                switch (optarg[0])
                {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognized\n", optarg);
                };
                break;
            case 'c': args->columns = optarg; break;
	    case 'a': args->targets_fname = optarg; break;
            case 'H': args->print_header = HEADER_ONLY; break;
            case 'h': args->print_header = HEADER_DROP; break;
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    char *fname = NULL;
    if ( optind>=argc ) {
        if ( !isatty(fileno((FILE*)stdin)) ) fname = "-";
        else usage(args);
    } else fname = argv[optind];

    if ( args->targets_fname )
    {
	htsFile *fp = hts_open(args->targets_fname,"r");
	htsFormat type = *hts_get_format(fp);
	hts_close(fp);
	if ( type.format==vcf || type.format==bcf )
	{
	    args->files->require_index = 1;
	    args->files->collapse |= COLLAPSE_SOME;
	}
	else error("Annotations should be VCF/BCF file\n");
    }
    
    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open %s: %s\n", fname, bcf_sr_strerror(args->files->errnum));

    init_data(args);
    while ( bcf_sr_next_line(args->files) )
    {
        if ( !bcf_sr_has_line(args->files, 0) ) continue;
        bcf1_t *line = bcf_sr_get_line(args->files, 0);
        if (line->errcode ) error("Encountered error, cannot processed. Please check the error output above.\n");
        annotate(args, line);
        bcf_write1(args->out_fh, args->hdr_out, line);
    }
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}
