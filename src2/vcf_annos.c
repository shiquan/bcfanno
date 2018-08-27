// Please note this code is from BCFtools
// All credit for this code goes to BCFtools's authors.
// Original copyright,
// Copyright (C) 2013-2016 Genome Research Ltd.
// Author : Petr Danecek <pd3@sanger.ac.uk>
//
#include "utils.h"
#include "anno_col.h"
#include "anno_vcf.h"
#include "vcmp.h"
#include "htslib/vcf.h"
#include "htslib/kstring.h"

int vcf_setter_filter(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, void *data)
{
    int i;
    bcf1_t *rec = (bcf1_t*) data;
    struct anno_vcf_buffer *b = f->buffer;
    if ( !(rec->unpacked & BCF_UN_FLT) ) bcf_unpack(rec, BCF_UN_FLT);
    if ( !(line->unpacked & BCF_UN_FLT) ) bcf_unpack(line, BCF_UN_FLT);
    if ( !rec->d.n_flt ) return 0;  // don't overwrite with a missing value

    if ( col->replace==SET_OR_APPEND || col->replace==REPLACE_MISSING ) {
        if ( col->replace==REPLACE_MISSING && line->d.n_flt ) return 0; // only update missing FILTER
        for ( i = 0; i < rec->d.n_flt; i++) {
            const char *flt = bcf_hdr_int2id(f->hdr, BCF_DT_ID, rec->d.flt[i]);
            bcf_add_filter(hdr, line, bcf_hdr_id2int(hdr, BCF_DT_ID, flt));
        }
        return 0;
    }
    hts_expand(int,rec->d.n_flt,b->mtmpi,b->tmpi);
    for ( i = 0; i < rec->d.n_flt; i++) {
        const char *flt = bcf_hdr_int2id(f->hdr, BCF_DT_ID, rec->d.flt[i]);
        b->tmpi[i] = bcf_hdr_id2int(hdr, BCF_DT_ID, flt);
    }
    bcf_update_filter(hdr,line,NULL,0);
    bcf_update_filter(hdr,line,b->tmpi,rec->d.n_flt);
    return 0;
}

int vcf_setter_id(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    //struct anno_vcf_buffer *b = f->buffer;
    // don't replace with "."
    if ( rec->d.id && rec->d.id[0]=='.' && !rec->d.id[1] )
        return 0;   
    if ( col->replace==SET_OR_APPEND )
        return bcf_add_id(hdr,line,rec->d.id);
    if ( col->replace!=REPLACE_MISSING )
        return bcf_update_id(hdr,line,rec->d.id);

   // running with +ID, only update missing ids
    if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) )
        return bcf_update_id(hdr,line,rec->d.id);
    return 0;
}
int vcf_setter_qual(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    //struct anno_vcf_buffer *b = f->buffer;
    if ( bcf_float_is_missing(rec->qual) ) return 0;
    if ( col->replace==REPLACE_MISSING && !bcf_float_is_missing(line->qual) ) return 0;
    line->qual = rec->qual;
    return 0;
}
int vcf_setter_info_flag(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    //struct anno_vcf_buffer *b = f->buffer;

    int flag = bcf_get_info_flag(f->hdr,rec,col->hdr_key,NULL,NULL);
    assert(flag >= 0);
    int ret = bcf_get_info_flag(hdr,line,col->hdr_key,NULL,NULL);
    if ( ret == -3 ) {
        bcf_update_info_flag(hdr,line,col->hdr_key,NULL,flag);
    } 
    return 0;
}
static int setter_ARinfo_int32(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, int nals, char **als, int ntmpi)
{
    struct anno_vcf_buffer *b = f->buffer;
    if ( col->number==BCF_VL_A && ntmpi!=nals-1 && (ntmpi!=1 || b->tmpi[0]!=bcf_int32_missing || b->tmpi[1]!=bcf_int32_vector_end) ) {
        error_return("Incorrect number of values (%d) for the %s tag at %s:%d", ntmpi, col->hdr_key, bcf_seqname(hdr,line), line->pos+1);
        return 1;
    } else if ( col->number==BCF_VL_R && ntmpi!=nals && (ntmpi!=1 || b->tmpi[0]!=bcf_int32_missing || b->tmpi[1]!=bcf_int32_vector_end) ) {
        error_return("Incorrect number of values (%d) for the %s tag at %s:%d", ntmpi, col->hdr_key, bcf_seqname(hdr,line), line->pos+1);
        return 1;
    }

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    int *map = vcmp_map_ARvalues(b->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) {
        error_return("REF alleles not compatible at %s:%d", bcf_seqname(hdr, line), line->pos +1);
        return 1;
    }
    
    // fill in any missing values in the target VCF (or all, if not present)
    int ntmpi2 = bcf_get_info_float(hdr, line, col->hdr_key, &b->tmpi2, &b->mtmpi2);

    if ( ntmpi2 < ndst )
        hts_expand(int32_t,ndst,b->mtmpi2,b->tmpi2);

    int i;
    for ( i = 0; i < ndst; i++) {
        if ( map[i] < 0 ) {
            if ( ntmpi2 < ndst ) b->tmpi2[i] = bcf_int32_missing;
            continue;
        }
        if ( ntmpi2==ndst && col->replace==REPLACE_MISSING && b->tmpi2[i]!=bcf_int32_missing && b->tmpi2[i]!=bcf_int32_vector_end )
            continue;

        b->tmpi2[i] = b->tmpi[ map[i] ];
    }
    return bcf_update_info_int32_fixed(hdr,line,col->hdr_key,b->tmpi2,ndst);
}

int vcf_setter_info_int(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    struct anno_vcf_buffer *b = f->buffer;
    if ( !(line->unpacked & BCF_UN_INFO) )
        bcf_unpack(line, BCF_UN_INFO);
    if ( !(rec->unpacked & BCF_UN_INFO) )
        bcf_unpack(rec, BCF_UN_INFO);
    
    int ntmpi = bcf_get_info_int32(f->hdr, rec, col->hdr_key, &b->tmpi, &b->mtmpi);

    if ( ntmpi < 0 ) return 0;    // nothing to add

    // check missing tag come first, changed by shiquan, 2018/01/30
    if ( col->replace==REPLACE_MISSING ) {    
        int ret = bcf_get_info_int32(hdr, line, col->hdr_key, &b->tmpi2, &b->mtmpi2);
        if ( ret>0 && b->tmpi2[0]!=bcf_int32_missing ) return 0;
    }
    
    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_int32(f,hdr,line,col,rec->n_allele,rec->d.allele,ntmpi);
   
    return bcf_update_info_int32_fixed(hdr,line,col->hdr_key,b->tmpi,ntmpi);
}
static int setter_ARinfo_real(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, int nals, char **als, int ntmpf)
{
    struct anno_vcf_buffer *b = f->buffer;
    if ( col->number==BCF_VL_A && ntmpf!=nals-1 && (ntmpf!=1 || !bcf_float_is_missing(b->tmpf[0]) || !bcf_float_is_vector_end(b->tmpf[0])) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpf,col->hdr_key,bcf_seqname(hdr,line),line->pos+1);
    else if ( col->number==BCF_VL_R && ntmpf!=nals && (ntmpf!=1 || !bcf_float_is_missing(b->tmpf[0]) || !bcf_float_is_vector_end(b->tmpf[0])) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpf,col->hdr_key,bcf_seqname(hdr,line),line->pos+1);

    
    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;

    int *map = vcmp_map_ARvalues(b->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%d\n", bcf_seqname(hdr, line), line->pos +1);

    // fill in any missing values in the target VCF (or all, if not present)
    int ntmpf2 = bcf_get_info_float(hdr, line, col->hdr_key, &b->tmpf2, &b->mtmpf2);
    if ( ntmpf2 < ndst )
        hts_expand(float,ndst,b->mtmpf2,b->tmpf2);

    int i;
    for (i=0; i<ndst; i++) {
        if ( map[i]<0 ) {
            if ( ntmpf2 < ndst )
                bcf_float_set_missing(b->tmpf2[i]);
            continue;
        }
        if ( ntmpf2==ndst && col->replace==REPLACE_MISSING
                && !bcf_float_is_missing(b->tmpf2[i])
                && !bcf_float_is_vector_end(b->tmpf2[i]) )
            continue;

        b->tmpf2[i] = b->tmpf[ map[i] ];
    }
    return bcf_update_info_float_fixed(hdr,line,col->hdr_key,b->tmpf2,ndst);
}
int vcf_setter_info_real(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    struct anno_vcf_buffer *b = f->buffer;

    if ( !(line->unpacked & BCF_UN_INFO) )
        bcf_unpack(line, BCF_UN_INFO);
    if ( !(rec->unpacked & BCF_UN_INFO) )
        bcf_unpack(rec, BCF_UN_INFO);
    
    int ntmpf = bcf_get_info_float(f->hdr,rec,col->hdr_key,&b->tmpf,&b->mtmpf);    
    if ( ntmpf < 0 ) return 0;    // nothing to add

    // check missing tag come first, changed by shiquan, 2018/01/30
    if ( col->replace==REPLACE_MISSING ) {
        int ret = bcf_get_info_float(hdr, line, col->hdr_key, &b->tmpf2, &b->mtmpf2);
        if ( ret>0 && !bcf_float_is_missing(b->tmpf2[0]) ) return 0;
    }

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_real(f,hdr,line,col,rec->n_allele,rec->d.allele,ntmpf);

    return bcf_update_info_float_fixed(hdr,line,col->hdr_key,b->tmpf,ntmpf);
}

static int copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst)
{
    int ith_src = 0, start_src = 0;    // i-th field in src string
    while ( ith_src<isrc && start_src<src_len ) {
        if ( src[start_src]==',' ) { ith_src++; }
        start_src++;
    }
    if ( ith_src!=isrc ) return -1; // requested field not found
    int end_src = start_src;
    while ( end_src<src_len && src[end_src] && src[end_src]!=',' ) end_src++;

    int nsrc_cpy = end_src - start_src;
    if ( nsrc_cpy==1 && src[start_src]=='.' )
        return 0;   // don't write missing values, dst is already initialized

    int ith_dst = 0, start_dst = 0;
    while ( ith_dst<idst && start_dst<dst->l ) {
        if ( dst->s[start_dst]==',' ) { ith_dst++; }
        start_dst++;
    }
    if ( ith_dst!=idst ) return -2;
    int end_dst = start_dst;
    while ( end_dst<dst->l && dst->s[end_dst]!=',' )
        end_dst++;

    if ( end_dst - start_dst>1 || dst->s[start_dst]!='.' )
        return 0;   // do not overwrite non-empty values

    // Now start_dst and end_dst are indexes to the destination memory area
    // which needs to be replaced with nsrc_cpy
    // source bytes, end_dst points just after.
    int ndst_shift = nsrc_cpy - (end_dst - start_dst);
    int ndst_move  = dst->l - end_dst + 1;  // how many bytes must be moved (including \0)
    if ( ndst_shift ) {
        ks_resize(dst, dst->l + ndst_shift + 1);    // plus \0
        memmove(dst->s+end_dst+ndst_shift, dst->s+end_dst, ndst_move);
    }
    memcpy(dst->s+start_dst, src+start_src, nsrc_cpy);
    dst->l += ndst_shift;
    return 0;
}

static int setter_ARinfo_string(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, int nals, char **als)
{
    int nsrc = 1, lsrc = 0;
    struct anno_vcf_buffer *b = f->buffer;
    while ( b->tmps[lsrc] ) {
        if ( b->tmps[lsrc]==',' ) nsrc++;
        lsrc++;
    }
    if ( col->number==BCF_VL_A && nsrc!=nals-1 && (nsrc!=1 || b->tmps[0]!='.' || b->tmps[1]!=0 ) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", nsrc,col->hdr_key,bcf_seqname(hdr,line),line->pos+1);
    else if ( col->number==BCF_VL_R && nsrc!=nals && (nsrc!=1 || b->tmps[0]!='.' || b->tmps[1]!=0 ) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", nsrc,col->hdr_key,bcf_seqname(hdr,line),line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    int *map = vcmp_map_ARvalues(b->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) {
        error_return("REF alleles not compatible at %s:%d", bcf_seqname(hdr, line), line->pos+1);
        return 1;
    }

    // fill in any missing values in the target VCF (or all, if not present)
    int i, empty = 0, nstr, mstr = b->tmpks.m;
    nstr = bcf_get_info_string(hdr, line, col->hdr_key, &b->tmpks.s, &mstr); 
    b->tmpks.m = mstr;
    if ( nstr<0 || (nstr==1 && b->tmpks.s[0]=='.' && b->tmpks.s[1]==0) ) {
        empty = 0;
        b->tmpks.l = 0;
        kputc('.',&b->tmpks);
        for ( i = 1; i < ndst; i++) kputs(",.",&b->tmpks);
    }
    else b->tmpks.l = nstr;
    for (i=0; i<ndst; i++) {
        if ( map[i]<0 ) {
            if ( empty ) copy_string_field(".",0,1,&b->tmpks,i);
            continue;
        }
        if ( col->replace==REPLACE_MISSING ) {
            // Do not replace filled values. The field must be looked up again because
            // of realloc in copy_string_field
            int n = 0;
            char *str = b->tmpks.s;
            while ( *str && n<i ) {
                if ( *str==',' ) n++;
                str++;
            }
            if ( str[0]!='.' || (str[1]!=',' && str[1]!=0) ) continue;  // value already set
        }
        int ret = copy_string_field(b->tmps,map[i],lsrc,&b->tmpks,i);
        assert( ret==0 );
    }
    return bcf_update_info_string_fixed(hdr,line,col->hdr_key,b->tmpks.s);
}
int vcf_setter_info_str(struct anno_vcf_file *f, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    struct anno_vcf_buffer *b = f->buffer;
    if ( !(line->unpacked & BCF_UN_INFO) )
        bcf_unpack(line, BCF_UN_INFO);
    if ( !(rec->unpacked & BCF_UN_INFO) )
        bcf_unpack(rec, BCF_UN_INFO);
    int ntmps = bcf_get_info_string(f->hdr,rec,col->hdr_key,&b->tmps,&b->mtmps);
    if ( ntmps < 0 ) return 0;    // nothing to add

    // check missing tag come first, changed by shiquan, 2018/01/30 
    if ( col->replace==REPLACE_MISSING ) {
        int ret = bcf_get_info_string(hdr, line, col->hdr_key, &b->tmps2, &b->mtmps2);
        if ( ret>0 && (b->tmps2[0]!='.' || b->tmps2[1]!=0) ) return 0;
    }
   
    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_string(f,hdr,line,col,rec->n_allele,rec->d.allele);
    
    return bcf_update_info_string(hdr,line,col->hdr_key,b->tmps);
}
