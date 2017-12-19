#include "anno_col.h"
#include "htslib/vcf.h"

void anno_col_clean(struct anno_col *c) {
    if ( c->hdr_key ) free(c->hdr_key);
    memset(c, 0, sizeof(*c));
}

void anno_col_copy(struct anno_col *src, struct anno_col *dest) {
    //anno_col_clean(dest);    
    dest->icol = src->icol;
    dest->replace = src->replace;
    dest->number = src->number;
    dest->hdr_key = strdup(src->hdr_key);
    dest->func = src->func;
    dest->curr_name = src->curr_name;
    dest->curr_line = src->curr_line;
}

//
static inline uint8_t *bcf_unpack_info_core1(uint8_t *ptr, bcf_info_t *info)
{
    uint8_t *ptr_start = ptr;
    info->key = bcf_dec_typed_int1(ptr, &ptr);
    info->len = bcf_dec_size(ptr, &ptr, &info->type);
    info->vptr = ptr;
    info->vptr_off  = ptr - ptr_start;
    info->vptr_free = 0;
    info->v1.i = 0;
    if (info->len == 1) {
        if (info->type == BCF_BT_INT8 || info->type == BCF_BT_CHAR) info->v1.i = *(int8_t*)ptr;
        else if (info->type == BCF_BT_INT32) info->v1.i = le_to_i32(ptr);
        else if (info->type == BCF_BT_FLOAT) info->v1.f = le_to_float(ptr);
        else if (info->type == BCF_BT_INT16) info->v1.i = le_to_i16(ptr);
    }
    ptr += info->len << bcf_type_shift[info->type];
    info->vptr_len = ptr - info->vptr;
    return ptr;
}

int bcf_update_info_fixed(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type)
{
    // Is the field already present?
    int i, inf_id = bcf_hdr_id2int(hdr,BCF_DT_ID,key);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,inf_id) ) return -1;    // No such INFO field in the header
    if ( !(line->unpacked & BCF_UN_INFO) ) bcf_unpack(line, BCF_UN_INFO);

    for (i=0; i<line->n_info; i++)
        if ( inf_id==line->d.info[i].key ) break;
    bcf_info_t *inf = i==line->n_info ? NULL : &line->d.info[i];

    if ( !n || (type==BCF_HT_STR && !values) )
    {
        if ( n==0 && !strcmp("END",key) )
            line->rlen = line->n_allele ? strlen(line->d.allele[0]) : 0;
        if ( inf )
        {
            // Mark the tag for removal, free existing memory if necessary
            if ( inf->vptr_free )
            {
                free(inf->vptr - inf->vptr_off);
                inf->vptr_free = 0;
            }
            line->d.shared_dirty |= BCF1_DIRTY_INF;
            inf->vptr = NULL;
            inf->vptr_off = inf->vptr_len = 0;
        }
        return 0;
    }

    // Encode the values and determine the size required to accommodate the values
    kstring_t str = {0,0,0};
    bcf_enc_int1(&str, inf_id);
    if ( type==BCF_HT_INT )
        bcf_enc_vint(&str, n, (int32_t*)values, -1);
    else if ( type==BCF_HT_REAL )
        bcf_enc_vfloat(&str, n, (float*)values);
    else if ( type==BCF_HT_FLAG || type==BCF_HT_STR )
    {
        if ( values==NULL )
            bcf_enc_size(&str, 0, BCF_BT_NULL);
        else
            bcf_enc_vchar(&str, strlen((char*)values), (char*)values);
    }
    else
    {
        error("The type %d not implemented yet", type);
    }

    // Is the INFO tag already present
    if ( inf )
    {
        // Is it big enough to accommodate new block?
        if ( str.l <= inf->vptr_len + inf->vptr_off )
        {
            if ( str.l != inf->vptr_len + inf->vptr_off ) line->d.shared_dirty |= BCF1_DIRTY_INF;
            uint8_t *ptr = inf->vptr - inf->vptr_off;
            memcpy(ptr, str.s, str.l);
            free(str.s);
            int vptr_free = inf->vptr_free;
            bcf_unpack_info_core1(ptr, inf);
            inf->vptr_free = vptr_free;
        }
        else
        {
            //assert( !inf->vptr_free );  // fix the caller or improve here: this has been modified before

            bcf_unpack_info_core1((uint8_t*)str.s, inf);
            inf->vptr_free = 1;
            line->d.shared_dirty |= BCF1_DIRTY_INF;
        }
    }
    else
    {
        // The tag is not present, create new one
        line->n_info++;
        hts_expand0(bcf_info_t, line->n_info, line->d.m_info , line->d.info);
        inf = &line->d.info[line->n_info-1];
        bcf_unpack_info_core1((uint8_t*)str.s, inf);
        inf->vptr_free = 1;
        line->d.shared_dirty |= BCF1_DIRTY_INF;
    }
    line->unpacked |= BCF_UN_INFO;

    if ( n==1 && !strcmp("END",key) ) line->rlen = ((int32_t*)values)[0] - line->pos;
    return 0;
}

