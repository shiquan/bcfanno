// anno_seqon.c -- Sequence Ontology and Human Genome Variantion Society nomenclature annotation
#include "anno_seqon.h"
#include "name_list.h"
#include "gea.h"
#include "sort_list.h"
#include "stack_lite.h"

static char *safe_duplicate_string(char *str)
{
    if (str == NULL ) return NULL;
    else return strdup(str);
}

#define MAX_GAP_GENE_DISTANCE 10000
#define MAX_TRANSCRIPT_LENGTH 100000
#define MAX_NORMLIZE_LENGTH   1000 // over 1000? cann't take care it anyway

// molecular consequence
struct MolecularConsequenceTerms {
    const char *lname;
    const char *sname;

    // rank to interpret the variant types
    //
    // to do : rank 
    // int rank;
};

const struct MolecularConsequenceTerms MCT[] = {
    { "unknown", "Unknown",},
    // intergenic variants
    { "intergenic_1KB_variant",                     "Intergenic", },   
    { "upstream_1K_of_gene",                        "Upstream1K", },
    { "upstream_10K_of_gene",                       "Upstream10K", },
    { "downstream_1K_of_gene",                      "Downstream1K", },
    { "downstream_10K_of_gene",                     "Downstream10K", },
    { "intergenic_variant",                         "Intergenic", },
    { "TFBS_ablation",                              "TFBS", },
    { "TFBS_variant",                               "TFBS", },

    // intragenoc variant, variant located inside gene but donot overlapped with any transcript
    { "intragenic_variant",                         "Intragenic", },
    // exon loss
    { "exon_loss",                                  "ExonLoss", },
    // intron variant
    { "noncoding_transcript_intron_variant",        "Intron", },
    { "coding_transcript_intron_variant",           "Intron", },
    { "3_prime_UTR_intron_variant",                 "Intron", },
    { "5_prime_UTR_intron_variant",                 "Intron", },
    // splice variant
    { "exonic_splice_region_variant",               "SpliceSites", },
    { "splice_donor_variant",                       "SpliceDonor", },
    { "splice_acceptor_variant",                    "SpliceAcceptor", },
    { "splice_donor_region_variant",                "SpliceDonor", },
    { "splice_polypyrimidine_tract_variant",        "SpliceAcceptor", },
    // noncoding variants
    { "noncoding_transcript_exon_variant",          "Noncoding", },
    { "noncoding_transcript_splice_region_variant", "SpliceSites", },
    // UTR variant
    { "5_prime_UTR_exon_variant",                   "UTR5", },
    { "3_prime_UTR_exon_variant",                   "UTR3", },
    { "5_prime_premature_start_codon_gain_varaint", "UTR5", },
    // coding variant
    { "frameshift_truncate",                        "Frameshift", },
    { "frameshift_elongation",                      "Frameshift", },
    { "inframe_indel",                              "InframeIndel", },
    { "disruptive_inframe_deletion",                "InframeDeletion", },
    { "disruptive_inframe_insertion",               "InframeInsertion", },
    { "conservative_inframe_deletion",              "InframeDeletion", },
    { "conservative_inframe_insertion",             "InframeInsertion", },
    { "start_loss",                                 "StartLoss", },
    { "start_retained_variant",                     "StartRetained", },
    { "stop_loss",                                  "StopLoss", },
    { "stop_retained_variant",                      "StopRetained", },
    { "stop_gained",                                "Nonsense", },
    { "missense_variant",                           "Missense", },
    { "synonymous_variant",                         "Synonymous", },
    { "nocall_variant",                             "NoCall", },
    { "transcript_ablation",                        "TranscriptAblation", },
    { "coding_sequence_variant",                    "CodingVariant",},
};

struct {
    int range; // splice site range
    int acceptor_region;
    //int splice_acceptor_range;
    int donor_region;
} splice_site_options = {  // splice site options
    .range = 3,
    .acceptor_region = 8,
    .donor_region = 3,
};

int file_is_GEA(const char *fn)
{
    return gea_check_format(fn);
}

extern struct mc_handler *mc_handler_init(const char *rna_fname, const char *data_fname, const char *reference_fname, const char *name_list);
extern struct mc_handler *mc_handler_duplicate(struct mc_handler *h);
extern void mc_handler_destroy(struct mc_handler *h, int l);

extern struct mc *mc_init(const char *chrom, int start, int end, char *ref, char *alt);
extern int mc_destroy(struct mc *h);
extern int mc_anno_trans(struct mc *n, struct mc_handler *h);
extern int mc_handler_fill_buffer_chunk(struct mc_handler *h, char* name, int start, int end);
extern int mc_anno_trans_chunk(struct mc *n, struct mc_handler *h);

struct mc_handler *mc_handler_init(const char *rna_fname, const char *data_fname, const char *reference_fname, const char *name_list)
{
    struct mc_handler *h = malloc(sizeof(*h));
    memset(h, 0, sizeof(*h));
    h->reference_fname = reference_fname;
    h->rna_fname = rna_fname;
    h->data_fname = data_fname;
    h->rna_fai = fai_load(rna_fname);
    
    if ( h->rna_fai == NULL ) error("Failed to load index of %s : %s.", rna_fname, strerror(errno));
    if ( gea_check_format(data_fname) ) error("Unsupported data format. Please try GenomeElementAnnotation format. %s.", data_fname);
    
    h->fp_idx = hts_open(data_fname, "r");
    if ( h->fp_idx == NULL ) error("%s : %s.", data_fname, strerror(errno));
    
    h->idx = tbx_index_load(data_fname);
    if ( h->idx == NULL ) error("Failed to load index of %s : %s.", data_fname, strerror(errno));

    h->hdr = gea_hdr_read(h->fp_idx);
    if ( h->hdr == NULL ) error("Failed to read header of %s.", data_fname);
    
    if ( name_list ) h->name_hash = name_hash_init(name_list);
    
    return h;
}

struct mc_handler *mc_handler_duplicate(struct mc_handler *h)
{
    struct mc_handler *d = malloc(sizeof(*d));
    memset(d, 0, sizeof(*d));
    d->rna_fname = h->rna_fname;
    d->data_fname = h->data_fname;
    d->reference_fname = h->reference_fname;
    
    d->rna_fai = fai_load(d->rna_fname);
    d->idx = tbx_index_load(d->data_fname);
    d->fp_idx = hts_open(d->data_fname, "r");

    d->hdr = h->hdr;
    d->name_hash = h->name_hash;    
    return d;
}
void mc_handler_destroy(struct mc_handler *h, int l)
{
    fai_destroy(h->rna_fai);
    tbx_destroy(h->idx);
    hts_close(h->fp_idx);
    if (l==0) gea_hdr_destroy(h->hdr);
    int i;
    for ( i = 0; i < h->n_record; ++i ) gea_destroy((struct gea_record*)h->records[i]);
    if ( h->n_record) free(h->records);
    free(h);
}

static int check_mitochondrial_chrom(const char *chrom)
{
    const char *mito = getenv("BCFANNO_MITOCHR");
    if ( mito == NULL ) 
        error("Failed to get mitochondrail chrom from ENVIR");
    return strcmp(chrom, mito) == 0;
}
/*
  Init variant.
*/
struct mc *mc_init(const char *chrom, int start, int end, char *ref, char *alt)
{
    struct mc *h = malloc(sizeof(*h));
    memset(h, 0, sizeof(*h));
   
    h->chr = chrom;
    h->start = start;
    h->end   = end;
    h->is_mito = check_mitochondrial_chrom(chrom);
    
    // if not consider alleles, only convert locations
    if ( ref == NULL && alt == NULL ) return h;
    // trim capped and tails
    int lr = strlen(ref);
    int la = strlen(alt);

    if ( ref && alt ) {
        // for GATK users, treat <NON_REF> as reference
        if ( strcmp(alt, "<NON_REF>") == 0 ) { h->type = var_type_nonref; return h;}

        for ( ;; ) {
            if ( ref == NULL || alt == NULL ) break;
            if (*ref == *alt) {
                ref++; alt++; lr--; la--;
                h->start++;
            }
            else if ( *(ref+lr-1) == *(alt+la-1) && lr > 0 && la > 0) {
                lr--, la--;
                h->end--;
            }
            else break;
        }

        if ( lr > 0) {
            h->ref = safe_duplicate_string(ref);
            h->ref[lr] = '\0';
        }
        if ( la > 0 ) {
            h->alt = safe_duplicate_string(alt);
            h->alt[la] = '\0';
        }
        if ( la == 0 ) {
            if ( lr == 0 ) h->type = var_type_ref;
            else h->type = var_type_del;
        } 
        else if ( lr == 0 ) h->type = var_type_ins;
        else {
            if ( la == 1 && lr == 1 ) h->type = var_type_snp;
            else h->type = var_type_delins;
        }
    }
    else if ( ref == NULL ) {
        h->type = var_type_ins;

        h->alt = safe_duplicate_string(alt);
    }
    else if ( alt == NULL ) {
        h->type = var_type_del;
        h->ref = safe_duplicate_string(ref);
    }
    if (*alt == 'N') {
        if ( h->ref ) free(h->ref);
        if ( h->alt ) free(h->alt);
        free(h);
        return NULL;
    }
    // for insert, end smaller than start.
    if ( lr == 0 && h->end < h->start ) {
        h->end = h->start;
        h->start--; // point start to capped base

        // there are several points should be notice here:
        // 1, mc::pos point to capped base, but ref point to a NULL address;
        // 2, capped base do NOT consider minus strand transcript, for minus insertion capped position will update to mc::end_pos
    }
    return h;
}

int mc_destroy(struct mc *h)
{
    if ( h == NULL ) return 0;
    if ( h->ref ) free(h->ref);
    if ( h->alt ) free(h->alt);
    int i;
    for ( i = 0; i < h->n_tran; ++i) {
        struct mc_inf *inf = &h->trans[i].inf;
        struct mc_type *type = &h->trans[i].type;
        if ( inf->transcript ) free(inf->transcript);
        if ( inf->gene ) free(inf->gene);
        if ( inf->ref ) free(inf->ref);
        if ( inf->alt ) free(inf->alt);
        if ( type->aminos ) free(type->aminos);        
    }
    if ( h->inter.gene ) free(h->inter.gene);
    free(h->trans);
    free(h);

    return 0;
}

static void clean_buffer_chunk(struct mc_handler *h)
{
    int i;
    for ( i = 0; i < h->n_record; ++i ) gea_destroy((struct gea_record*)h->records[i]);
    h->n_record = 0;
    h->i_record = 0;
    h->end_pos_for_skip = 0;
    h->last_gene = NULL;
    h->next_gene = NULL;
}

struct list_buffer {
    struct list_buffer *next;
    struct gea_record  *data;
};
struct list_buffer *list_buffer_create()
{
    struct list_buffer *l = calloc(1, sizeof(*l));
    return l;
}

// return length of list
static int retrieve_gea_records_from_region(struct mc_handler *h, int id, int start, int end, struct list_buffer **header, struct list_buffer **tail, int *tail_edge)
{
    // retrieve annotation records from database
    hts_itr_t *itr = tbx_itr_queryi(h->idx, id, start, end+1);
    kstring_t string = {0,0,0};
    int l = 0;
    *tail_edge = 0;
    
    while ( tbx_itr_next(h->fp_idx, h->idx, itr, &string) >= 0 ) {
        struct gea_record *r = gea_init();
        if ( gea_parse(&string, h->hdr, r) ) continue;
        const char *type = h->hdr->id[GEA_DT_BIOTYPE][r->biotype].key;
        if ( h->name_hash) {
            if (strcmp(type, "mRNA") == 0 || strcmp(type, "ncRNA") == 0 ) {
                if (!name_hash_key_exists(h->name_hash, r->name)) continue;
            }
            else if ( strcmp(type, "gene") == 0 ) {
                if (!name_hash_key_exists(h->name_hash, r->geneName)) continue;
            }
        }
              
        gea_unpack(h->hdr, r, GEA_UN_TRANS|GEA_UN_CIGAR);
        if ( r->chromEnd > *tail_edge ) *tail_edge = r->chromEnd;
        
        struct list_buffer *b = list_buffer_create();
        b->data = r;
        
        if ( *header == NULL ) *header = b;
        if ( *tail == NULL ) *tail = *header;
        else {
            (*tail)->next = b;  *tail = b;
        }
        l++;
    }
    free(string.s);
    tbx_itr_destroy(itr);
    
    return l;
}
static int is_gene(const struct gea_hdr *h, struct gea_record *v)
{
    char *bt = (char*)h->id[GEA_DT_BIOTYPE][v->biotype].key;
    if ( strcmp(bt, "Gene") == 0 || strcmp(bt, "mRNA") == 0 || strcmp(bt, "ncRNA") == 0 ) return 1;
    return 0;
}
// Update buffer for each chunk.
int mc_handler_fill_buffer_chunk(struct mc_handler *h, char* name, int start, int end)
{
    extern int name_hash_key_exists(void *hash, char *key); 
    
    // clean buffer
    clean_buffer_chunk(h);

    //debug_print("%s:%d-%d", name, start, end);
    
    int id;
    id = tbx_name2id(h->idx, name);
    if ( id == -1 ) return 0;
    
    int l;
    struct list_buffer *header = NULL, *tail = NULL;
    int tail_edge;
    l = retrieve_gea_records_from_region(h, id, start, end, &header, &tail, &tail_edge);

    // 2018-08-10, update filter name list    
    
    // if retrieved buffer did NOT cover the chunk, trying to find more upstream records, used to
    // annotate intergenic variants
    if ( header == NULL || header->data->chromStart > start) {
        int l1;
        struct list_buffer *s = NULL, *t = NULL;
        int edge;
        l1 = retrieve_gea_records_from_region(h, id, start - MAX_GAP_GENE_DISTANCE, start, &s, &t, &edge);
        if ( l1 > 0 ) {
            struct list_buffer *p = s, *g = NULL;
            for ( ;; ) {
                if ( p == NULL ) break;
                // temp point 
                struct list_buffer *q = NULL;                
                if ( is_gene(h->hdr, (struct gea_record*)p->data) ) { q = g; g = p; }
                else q = p;
                p = p->next;
                if (q) {
                    gea_destroy((struct gea_record*)q->data);
                    free(q);
                }
            }
            // update header
            if ( g) {
                g->next = header;
                header = g;
            }
        }
    }

   /* // if the chunk was NOT fully covered, find more downstream records */
    if ( tail == NULL || tail_edge < end ) {
        int l2;
        struct list_buffer *s = NULL, *t = NULL;
        int edge;

        // here we retrieve all records downstream, better idea could be read downstream record one by one and check the
        // biotype at the same time. Because we only need to the most nearest gene downstream.
        l2 = retrieve_gea_records_from_region(h, id, end, end + MAX_GAP_GENE_DISTANCE, &s, &t, &edge);
        
        if ( l2 > 0 ) {
            struct list_buffer *p = s, *g = NULL;
            for ( ;; ) {
                if ( p == NULL ) break;
                // temp point
                if ( g == NULL && is_gene(h->hdr, (struct gea_record*)p->data) ) { g = p; p = p->next; g->next = NULL;}
                else {
                    struct list_buffer *q = NULL;
                    q = p;
                    p = p->next;                    
                    gea_destroy((struct gea_record*)q->data);
                    free(q);
                }
            }
            // g could be NULL
            if ( g) {                
                if (tail == NULL ) { tail = g; header = tail;}
                else { tail->next = g; tail = g;}

            }
        }
    }

    if ( header == NULL ) return 0;
    
    // convert list into buffer
    int i, n;
    n = count_list((const void*)header);
    if ( n > h->m_record ) {
        h->m_record = n;
        h->records = realloc(h->records, h->m_record*sizeof(struct gea_record));
    }
    h->n_record = n;

    for ( i = 0; i < n; ++i ) {
        h->records[i] = header->data;
        struct list_buffer *l = header;
        header = header->next;
        free(l);
    }
    
    return h->n_record;
}

/*
 * Variants located in the gene region, protein changes will be checked.
 */

// find the block of variants located, could be overlapped with multiblocks, so both start and end block will be exported
// dist: distance to exon edges, 
static int find_the_block(struct gea_record *v, int *s, int *e, int pos, int32_t *dist)
{
    *s = 0; *e = 0; *dist = 0;
    if ( pos < v->chromStart || pos > v->chromEnd ) return 1;
    int i;
    int left;
    int right;
    int is_exon = 0;
    for ( i = 0; i < v->blockCount; ++i ) {
        int start = v->blockPair[0][i];
        int end   = v->blockPair[1][i];
        left = pos-start;
        right = end - pos;
        // update 2019/06/29: intron count should -1 thereafter
        if ( pos < start) { *e = i; break; }
        else {
            *s = i;
            if ( pos <= end ) { *e = i; is_exon = 1; break; }
        }
    }
    // out of range
    if ( i == v->blockCount && *s > 0 && *e == 0 ) return 1;
    if (is_exon) { // check the distance
        assert(left>=0 && right >=0);
        *dist = left<<16|(uint16_t)right;
    }
    return 0;
}

// find the position on the exon region and offset if variant located in the intron region
// 0 on success, 1 on out of range
static int find_locate(struct gea_hdr *hdr, struct gea_record *v, int *pos, int *offset, int start, int *id, int *dist_exon)
{
    gea_unpack(hdr, v, GEA_UN_TRANS|GEA_UN_CIGAR);
    
    int b1 = -1, b2 = -1;
    *pos = 0;
    *offset = 0;
    // for location out of range set pos and offset to 0, annotate as '?'
    if ( find_the_block(v, &b1, &b2, start, dist_exon ) ) { *pos = 0, *offset = 0; return 1;}

    struct gea_coding_transcript *c = &v->c;
    if ( b1 == b2 ) {
        if ( v->strand == strand_is_plus ) *pos = c->loc[0][b1] + start - v->blockPair[0][b1];
        else *pos = c->loc[1][b1] + (v->blockPair[1][b1] - start);
        *offset = 0;
    }
    else {
        int upstream = start - v->blockPair[1][b1];
        int downstream = v->blockPair[0][b2]-start;
        if ( upstream > downstream ) {
            *pos = c->loc[0][b2];
            *offset = v->strand == strand_is_plus ? -downstream : downstream;
        }
        else {
            *pos = c->loc[1][b1];
            *offset = v->strand == strand_is_plus ? upstream : -upstream;
        }
    }
    // exon ID
    *id = b1;

    // adjust position if located in exon and there are realignments in the transcript sequence
    if ( *offset != 0 ) return 0;
    if ( v->n_cigar == 0 ) return 0;
    if ( v->n_cigar == 1 && CIGAR_EQUAL(v->cigars[0], CIGAR_MATCH_TYPE) ) return 0;

    int adjust = 0;
    int i;
    int match = 0;
    int ins, del;

    for ( i = 0; i < v->n_cigar; ++i ) {            

        if ( CIGAR_EQUAL(v->cigars[i],CIGAR_MATCH_TYPE) ) match += v->cigars[i] >> CIGAR_PACKED_FIELD;            
        else if ( CIGAR_EQUAL(v->cigars[i],CIGAR_DELETE_TYPE) ) {
            del = v->cigars[i] >> CIGAR_PACKED_FIELD;
            // there is no need to check match and *pos, becase for plus strand match will be always smaller than *pos
            // check if this deletion in the target block
            if ( c->loc[v->strand == strand_is_plus? 0 :1][b1] <= match && *pos > match) {
            //if ( c->loc[0][b1] <= match && *pos > match) {
                adjust -= del;
                // if this variant located in the deletion, just put pos to the edge of this gap
                if ( *pos < match +  del ) { *pos += 1; return 0; }
            }
            match += del;
        }
        else if ( CIGAR_EQUAL(v->cigars[i], CIGAR_INSERT_TYPE) ) {
            ins = v->cigars[i] >> CIGAR_PACKED_FIELD;
            if ( c->loc[v->strand == strand_is_plus ? 0 :1][b1] <= match  && *pos > match ) adjust += ins;
        }
        
        //if ( v->strand == strand_is_plus && match >= *pos ) break;            
        //else if ( v->strand == strand_is_minus && match < *pos ) break;

        // if (match >= *pos ) break; in case variant happens at offset location
        if (match > *pos ) break;
        
        // locations in GEA structure has been adjusted, so if matched region in the upstream of the target exon, reset adjust to 0
        if ( v->strand == strand_is_plus && match + adjust < c->loc[0][b1] ) adjust = 0;
        else if ( v->strand == strand_is_minus && match + adjust < c->loc[1][b1] ) adjust = 0;
    }
    *pos += adjust;
    
    return 0;
}
static int trim_capped_sequences(char **_ref, char **_alt, int *lr, int *la, int unit)
{
    assert(unit > 0);
    int i;
    int l = *lr < *la ? *lr : *la;
    char *ref = *_ref;
    char *alt = *_alt;
    //int trimmed;    
    for ( i = 0; i < l; ) {
        if ( l - i < unit ) break;
        if ( strncmp(ref+i, alt+i, unit) != 0 ) break;
        i += unit;
    }
    if ( i > 0 ) {
        *lr -= i;
        *la -= i;
        int j;
        for ( j = 0; j < *lr; ++j ) ref[j] = ref[i+j];
        ref[*lr] = '\0';
        for ( j = 0; j < *la; ++j ) alt[j] = alt[i+j];
        alt[*la] = '\0';
    }
    
    return i;
}
// this function construct the full alternative sequence, do NOT care is frameshift or not
// for minus strand inserted based should be complemented first
static int construct_alternative_sequence(char *ref, char *alt, int lref, int lalt, int pos, char *ori, int lori, char **mut, int *lmut)
{
    *lmut = lori - lref + lalt;
    assert(*lmut >= 0);
    if ( *lmut == 0 ) {
        *mut = NULL;
        return 0;
    }
    *mut = malloc(sizeof(char)*(*lmut+1));

    int i;
    char *p = *mut;
    for ( i = 0; i < pos; ++i, ++p ) *p = ori[i];
    // for insertion, pos point to capped base
    if ( lref == 0) {
        *p = ori[pos];
        p++;
    }
    for ( i = 0; i < lalt; ++i, ++p ) *p = alt[i];
    if ( lref == 0 )         
        for ( i = pos+1; i < lori; ++i, ++p) *p = ori[i];
    else
        for ( i = pos+lref; i < lori; ++i, ++p) *p = ori[i];

    *p = '\0';
    
    return *lmut;
}
// we don't count end location here, because we can easily interupt it with original sequence length
/* static void update_variant_location(struct mc_inf *inf, struct mc_type *type, struct gea_record *v) */
/* { */
/*     struct gea_coding_transcript *c = &v->c; */
/*     // Update the variant location of domain and functional region     */
/*     if ( inf->pos <= c->utr5_length ) { */
/*         type->func1 = func_region_utr5_exon; */
/*         inf->loc = c->utr5_length - inf->pos +1; */
/*         type->con1= mc_utr5_exon; */
/*     } */
/*     else if ( inf->pos <= c->cds_length ) { */
/*         type->func1 = func_region_cds; */
/*         inf->loc = inf->pos - c->utr5_length; */
/*         type->loc_amino = (inf->loc+3)/3; */
/*     } */
/*     else { */
/*         type->func1 = func_region_utr3_exon; */
/*         inf->loc = inf->pos - c->cds_length; */
/*         type->con1 = mc_utr3_exon; */
/*     } */
/* } */

static int predict_molecular_consequence_snp(struct mc *mc, struct mc_type *type, struct mc_inf *inf, struct gea_record *v, int cod, char *ref, char *alt, char *ori, char *mut)
{
    struct gea_coding_transcript *c = &v->c;
    
    // location of amino acid
    inf->loc = inf->pos - c->utr5_length;
    type->loc_amino = (inf->loc+2)/3;
    // if ( lref == lalt && lref == 1 ) {
    type->ori_amino = codon2aminoid(ori,mc->is_mito);

    // no-call
    if ( ori[cod] == *alt ) {
        type->mut_amino = type->ori_amino;
        type->con1 = mc_nocall;
        if ( *mc->ref != ori[cod] ) inf->ref = strndup(ori+cod,1);
        if ( *mc->alt != *alt ) inf->alt = safe_duplicate_string(alt);
        return 0;
    }
    else if ( ori[cod] != *ref) warnings("Inconsistance nucletide between refenence and transcript: %s,%d,%s,%d,%c,%c ", mc->chr, mc->start, v->name, inf->pos, ori[cod], *ref);
    
    // normal case
    type->ori_amino = codon2aminoid(ori,mc->is_mito);
    type->mut_amino = codon2aminoid(mut,mc->is_mito);

    if ( type->ori_amino == type->mut_amino) {
        type->con1 = mc_synonymous;
        // if ( type->ori_amino == C4_Stop ) type->con1 = mc_stop_retained;
        // if ( type->loc_amino == 1) type->con1 = mc_start_retained;        
    }
    else {
        type->con1 = mc_missense;
        // if ( type->ori_amino == C4_Stop ) type->con1 = mc_stop_loss;
        // if ( type->mut_amino == C4_Stop ) type->con1 = mc_stop_gained;
        // else if ( type->loc_amino == 1 ) type->con1 = mc_start_loss;
    }

    if ( *mc->ref != *ref ) inf->ref = safe_duplicate_string(ref);
    if ( *mc->alt != *alt ) inf->alt = safe_duplicate_string(alt);
    
    return 0;
}

static int trim_amino_acid_ends(int *lori_aa, int *ori_aa, int *lmut_aa, int *mut_aa, int *head, int *tail, int e)
{
    int i;
    //for (i = 0; i < *lori_aa && i < *lmut_aa && (*ori_aa)[i] == (*mut_aa)[i]; ++i);
    for (i = 0; i < *lori_aa && i < *lmut_aa && ori_aa[i] == mut_aa[i]; ++i);
    //if ( *lori_aa == i || *lmut_aa == i) break;
    //  if ( (*ori_aa)[i] != (*mut_aa)[i] ) break;
    *head = i;
    *tail = 0;
    if (i) {
        *lori_aa = *lori_aa -i;
        *lmut_aa = *lmut_aa -i;
        //int j;
        //for ( j = 0; j < *lori_aa; ++j) ori_aa[j] = ori_aa[i+j];
        //for ( j = 0; j < *lori_aa; ++j) mut_aa[j] = mut_aa[i+j];
        memmove(ori_aa, ori_aa + i, *lori_aa *sizeof(int));
        memmove(mut_aa, mut_aa + i, *lmut_aa *sizeof(int));
    }
    if ( e == 2 ) {
        for ( i = 0;; ++i ) {
            if ( *lori_aa -i == 0 || *lmut_aa -i == 0 ) break;
            if ( ori_aa[*lori_aa -i-1] != mut_aa[*lmut_aa -i-1] ) break;
        }
        *tail = i;
    }
    return *head+*tail;
}
// todo: consider Selenocysteine!
// max_len is defined length of coding region, if not set -1, UGA codon will be set to Selenocysteine
static int *convert_bases_to_amino_acid(int l, char *s, int *l_aa, int mito, int max_len)
{
    *l_aa = 0;
    if (l<3) return NULL;
    
    if (max_len == -1) { // do not check Sec
        int i;
        int *aa = (int*)calloc(l/3,sizeof(int));
        for (i=0; i<l/3; i++) {
            aa[i] = codon2aminoid(s+i*3,mito);
            if ( aa[i] == C4_Stop ) {
                i++; // put stop codon at the end
                break;
            }
        }
        *l_aa = i;
        return aa;
    }
    else {        
        assert(max_len <= l);
        int i;
        int *aa = (int*)calloc(l/3, sizeof(int));
        for (i = 0; i < max_len/3-1; i++) {
            char *cod = s+i*3;
            aa[i] = codon2aminoid(cod, mito);
            if (aa[i] == C4_Stop) {
                if ((*cod == 'T' || *cod == 'U') && *(cod+1) == 'G' && *(cod+2) == 'A') aa[i] = C4_Sec; // Selenocysteine
                else break;
            }            
        }
        
        // check stop codon
        if (i == max_len/3-1) { // todo: confused transcript, need double check here.
            for (; i < l/3; ++i) {
                aa[i] = codon2aminoid(s+i*3, mito);
                if (aa[i] == C4_Stop) break;
            }
        }
        *l_aa = i+1;
        return aa;
    }
}
//
// the offset only used to realign the alternative sequence and reference sequence, but would NOT change the state of variant types
//
// ex1: insert CGA
//  ATCGACGAT             GAT..      ## at this moment should be interpret as insert GAC instead of CGA consider of read frames
//  ATCGACGACGAT     ==>  GACGAT...
//
// ex2: insert CG
//  ATCGCGCTA             CTA..
//  ATCGCGCGCTA      ==>  CGCTA..
// 
// ex3: delete CGA                   ## at this moment should interpret as delete GAC instead of CGA consider of read frames
//  ATCGACGAT             GACGAT..
//  ATCGAT           ==>  GAT..
//
// ex4: delete C
//  ATCCACGAT             CACGAT.
//  ATCACGAT         ==>  ACGAT..
//
// ex5: delete CGA
//  ATCGATGAT             ATCGATGAT
//  ATTGAT           ==>  ATTGAT..
//  
// ex5: delins G   ## since ref and alt bases have been trimmed during mc_init(), offset should always be 0 for delins
//  ATCGACGAT             ATCGACGAT..
//  ATGCGAT          ==>  ATGCGAT..
// After trimmed capped sequence, location of amino acid will be changed, even the sequences can also be updated, but the state and the length of 
// insert or delete bases will be the same.    
//
// So only need to update the positions and affected bases
//
static int trim_capped_sequences_and_check_mutated_end(struct mc_handler *h, struct mc *mc, struct mc_type *type, struct mc_inf *inf, struct gea_record *v, int lref, int lalt, char *ref, char *alt, int *lori, char **ori, int *lmut, char **mut, int *cod, int cap)
{
    struct gea_coding_transcript *c = &v->c;

    // transcript level
    // trim Longest Common 3mer Prefixes
    int offset = trim_capped_sequences(ori, mut, lori, lmut, 3);
    
    if (lref > 0) inf->ref = strdup(ref);
    if (lalt > 0) inf->alt = strdup(alt);

    //if (offset) {
        // adjust offset by trim capped bases, cod is 0 based, convert offset to 0based by -1 
        // offset = offset - cod -1;

    int i;
    for ( i=0; i<*lmut && i<*lori; ++i) {
        if ( (*mut)[i] != (*ori)[i] ) break;
    }
    assert(i<3);
    offset += i;

    if ( offset > 0 && offset-*cod > 0) {
        //
        offset =cap == 0 ? offset - *cod : offset-*cod -1;
        inf->pos = inf->pos + offset;
        
        if ( inf->pos > c->cds_length ) {
            type->con1 = mc_utr3_exon;
            type->func1 = func_region_utr3_exon;
            inf->loc = inf->pos - c->cds_length;
            inf->end_loc = inf->end_pos + offset - c->cds_length;
            return 0; // out of coding region, no need to check the amino acid changes any more.
        }
        inf->loc += offset;
        inf->end_loc += offset;
        inf->pos += offset;
        inf->end_pos += offset;
        if ( inf->pos > v->c.cds_length) {
            type->func1 = type->func2 = func_region_utr3_exon;
            type->con1 = mc_utr3_exon;
            inf->loc = inf->pos - v->c.cds_length;
            inf->end_loc = inf->end_pos - v->c.cds_length;
            return 0;
        }
        else if ( inf->end_pos > v->c.cds_length ) {
            type->func2 = func_region_utr3_exon;
            type->con1 = mc_exon_loss;
            inf->end_loc = inf->end_pos - v->c.cds_length;
            return 0;
        }
        // update cod
        *cod = (inf->loc-1)%3;

        //type->loc_amino = (inf->loc+2)/3;

        // update ref and alt sequence
        if ( lref > 0 && strncmp(*ori+i, ref, lref) != 0 ) memcpy(inf->ref, *ori+i, lref);        
        if ( lalt > 0 && strncmp(*mut+i, alt, lalt) != 0 ) memcpy(inf->alt, *mut+i, lalt);
    }

    // Sometime after round the codon frame, length of alternative sequence may smaller than 3, in this case transcript has no UTR3 and
    // there is a deletion at the stop codon. At this stituation, retrieve sequence from reference, and concat the retrieved tail to
    // alternative sequence.    
    if ( *lmut < 3 ) {
        // expand 10 base in the downstream in the reference sequence
        faidx_t *fai = fai_load(h->reference_fname);
        if (fai) {
            char *seq;
            int   l;
            if ( v->strand == strand_is_plus ) seq = faidx_fetch_seq(fai, mc->chr, v->cEnd, v->cEnd+10, &l);
            else {
                seq = faidx_fetch_seq(fai, mc->chr, v->cStart - 10, v->cStart, &l);
                compl_seq(seq, l);
            }
            if ( l && seq ) {
                kstring_t str = {0,0,0};
                kputs(*mut, &str);
                kputs(seq, &str);
                type->mut_amino = codon2aminoid(str.s,mc->is_mito);
                free(str.s);
            }
            
            if ( type->mut_amino == C4_Stop) type->con1 = mc_stop_retained;
            else type->con1 = mc_stop_loss;
            
            fai_destroy(fai);
            return 0;
        }
        else {
            // even no reference sequence specified, we could interpret it as a stop-loss
            warnings("%s : %s.", h->reference_fname, strerror(errno));
            type->con1 = mc_stop_loss;
            return 0;
        }
    }
    
    return 1;
}
// For deletion, lmut should always be 0, so no need to check it here
static int predict_molecular_consequence_deletion(struct mc_handler *h, struct mc *mc, struct mc_type *type, struct mc_inf *inf, struct gea_record *v, int lref, char *ref, int lori, int lmut, char *ori, char *mut)
{
    int cod;
    cod = (inf->loc-1)%3;
        
    // re-alignment of reference and alternative alleles
    if ( trim_capped_sequences_and_check_mutated_end(h, mc, type, inf, v, lref, 0, ref, NULL, &lori, &ori, &lmut, &mut, &cod, 0) == 0 ) return 0;
    
    // location of amino acids
    type->loc_amino = (inf->loc+2)/3; // 1 based
    type->loc_end_amino = (inf->end_loc+2)/3;
    
    int *ori_aa, *mut_aa;
    int lori_aa, lmut_aa;

    int max_len = v->c.cds_length - v->c.utr5_length; // real cds length
    max_len = max_len - (type->loc_amino-1)*3;
    // max_len = (max_len+2)/3*3; // weird adjust
    // debug_print("%d\t%d\t%d\t%d\t%s", v->c.cds_length, v->c.utr5_length, inf->loc, max_len, v->name);
    ori_aa = convert_bases_to_amino_acid(lori, ori, &lori_aa, mc->is_mito, max_len);
    mut_aa = convert_bases_to_amino_acid(lmut, mut, &lmut_aa, mc->is_mito, -1);
    if (ori_aa == NULL) error("%s", v->name);
    if ( lori_aa == 0 ) return 0;
    
    // check the amino acids length
    //if ( lori < 10000 && ((v->c.cds_length - v->c.utr5_length - inf->loc)/3 +1 > lori_aa )) {
    // debug_print("%d\t%d\t%d\t%d", v->c.cds_length, v->c.utr5_length, inf->loc, lori_aa);

    if ( lori < MAX_GAP_GENE_DISTANCE && (max_len/3 > lori_aa )) {
        inf->inframe_stop = 1;
        //warnings("In frame stop codon found. %s, %s:%d", inf->transcript, mc->chr, mc->start);
    }
    
    if ( lref%3 ) { // frameshift
        type->con1 = mc_frameshift_truncate;
        
        int head, tail;
        trim_amino_acid_ends(&lori_aa, ori_aa, &lmut_aa, mut_aa, &head, &tail, 1);
        if ( head ) type->loc_amino += head; // if trimmed head, update new location
        type->ori_amino = ori_aa[0];
        if ( lmut_aa == 0 ) {
            type->fs = -1;
            return 0;
        }
        if ( type->ori_amino == C4_Stop ) type->con1 = mc_frameshift_elongation;
        type->mut_amino = mut_aa[0];
        type->fs = mut_aa[lmut_aa-1] == C4_Stop ? lmut_aa : -1;
        //debug_print("%d,%d", type->mut_amino, lmut_aa);
    }
    else { // inframe deletion
        
        if ( cod == 0 ) {
            /*
              For inframe deletions, there is NO need to consider the length of deleted amino acids, since no amino acid codon will be broken during this process. We should only find the first changed
              codon as start codon, and update the start and end positions on the amino acids if necessary.
            */
            int l = lref/3-1; // 0 based
            if ( l >= lori_aa ) {
                warnings("Amibigious interpret of transcript. %s, %s:%d", inf->transcript, mc->chr, mc->start);
                l = lori_aa -1 > 0 ? lori_aa -1 : 0;
            }
            type->con1 = mc_conservation_inframe_deletion;
            int head, tail;
            trim_amino_acid_ends(&lori_aa, ori_aa, &lmut_aa, mut_aa, &head, &tail, 1);
            type->loc_amino += head;
            type->loc_end_amino = type->loc_amino + l;
            type->ori_amino = ori_aa[0];
            type->ori_end_amino =ori_aa[l];
            type->n = 0; // check again?? 2019/07/01
        }
        else {
            /*
              For disruptive deletions, original codons will be broken, and a new codon will be create after deletion. So we should consider the deleted codons and mutated codon, check it is a delins
              or deletion variant.
             */
            type->con1 = mc_disruption_inframe_deletion;
            int head, tail;
            trim_amino_acid_ends(&lori_aa, ori_aa, &lmut_aa, mut_aa, &head, &tail, 1);

            if ( head > 0 ) { // mutated codon has been trimmed, so only treat it as standard deletion, else we should treat it as del-ins
                type->loc_amino += head;
                int l = lref/3 -1; // 0 based
                type->loc_end_amino = type->loc_amino + l;
                type->ori_amino = ori_aa[0];
                type->ori_end_amino = ori_aa[l];
                
                type->n = 0;
            }
            else {
                // now we also should check if mutated codon same with end codon, if it is same, treat it as deletioon else del-ins
                int l = lref/3; // 0 based
                if ( ori_aa[l] == mut_aa[0] ) {
                    l -= 1;
                    type->loc_end_amino = type->loc_amino +l;
                    type->ori_amino = ori_aa[0];
                    type->ori_end_amino = ori_aa[l];
                    type->n = 0;
                }
                else { 
                    type->loc_end_amino = type->loc_amino +l;
                    type->ori_amino = ori_aa[0];
                    type->ori_end_amino = ori_aa[l];
                    type->n = 1;
                    type->aminos = calloc(1,sizeof(int));
                    type->aminos[0] = mut_aa[0];
                }
            }
        }

    }
    if ( lori_aa ) free(ori_aa);
    if ( lmut_aa ) free(mut_aa);

    return 0;
}

static int predict_molecular_consequence_delins(struct mc_handler *h, struct mc *mc, struct mc_type *type, struct mc_inf *inf, struct gea_record *v, int lref, char *ref, int lalt, char *alt, int lori, int lmut, char *ori, char *mut)
{
    int cod = (inf->loc -1)%3;
    // re-alignment of reference and alternative alleles
    if ( trim_capped_sequences_and_check_mutated_end(h, mc, type, inf, v, lref, lalt, ref, alt, &lori, &ori, &lmut, &mut, &cod, 0) == 0 ) return 0;

    inf->end_loc = inf->loc + lref -1;

    type->loc_amino = (inf->loc+2)/3;
    type->loc_end_amino = (inf->end_loc+2)/3;

    //if ( strncmp(ref, ori, lref) != 0 ) inf->ref = strndup(ori, lref);
    //if ( strncmp(alt, mut, lalt) != 0 ) inf->alt = strndup(mut, lalt);
    
    int *ori_aa, *mut_aa;
    int lori_aa, lmut_aa;
    int max_len = v->c.cds_length - v->c.utr5_length - inf->loc;
    ori_aa = convert_bases_to_amino_acid(lori, ori, &lori_aa, mc->is_mito, max_len);
    mut_aa = convert_bases_to_amino_acid(lmut, mut, &lmut_aa, mc->is_mito, -1);

    //assert(lori_aa >0);
    if ( lori_aa == 0 ) return 0;

    if ( lref == lalt ) {
      inframe_indel:
        type->con1 = mc_inframe_indel;
        //int l1 = (cod + lref)/3*3+1; // 1 based
        int l1 = (cod + lref-1)/3+1; // 1 based
        //int l2 = (cod + lalt)/3*3+1;
        int l2 = (cod + lalt-1)/3+1;
        //assert(lori_aa>l1);
        // update:2019/06/24
        // l1 is change position, if outrange of coding region, return -1
        if ( l1 > lori_aa ) return -1;
        
        // stop codon founded the inserted codon when lmut_aa < l2
        
        int head, tail;
        trim_amino_acid_ends(&l1, ori_aa, &l2, mut_aa, &head, &tail, 2);
        type->loc_amino += head;
        type->loc_end_amino = type->loc_amino + l1 -1;
        type->ori_amino = ori_aa[0];
        type->ori_end_amino = ori_aa[l1-1];
        // update: inframe delins 2019/06/29
        if (l1 == 1 && l2 == 1) { // A>B
            type->n = 1;
            type->mut_amino = mut_aa[0];
        }
        else { // AdelinsB
            type->n = l2;
            type->aminos = malloc(sizeof(int)*type->n);
            int i;
            for ( i=0; i<type->n; ++i) {
                type->aminos[i] = mut_aa[i];
                if ( type->aminos[i] == 0 ) { // stop-codon
                    //type->con1 = mc_frameshift_truncate;
                    type->con1 = i == 0 ? mc_stop_gained : mc_frameshift_truncate;
                    type->fs = i+1;
                    type->n = i+1;
                    type->mut_amino = mut_aa[0];
                    break;
                }
            }
        }
    }
    else {
        if ( (lref-lalt)%3 == 0 ) goto inframe_indel;
        // p.Arg267AlnfsXX
        type->con1 = lref > lalt ? mc_frameshift_truncate : mc_frameshift_elongation;
        int head, tail;
        trim_amino_acid_ends(&lori_aa, ori_aa, &lmut_aa, mut_aa, &head, &tail, 1);
        type->loc_amino += head;

        type->ori_amino = ori_aa[0];  
        type->mut_amino = mut_aa[0];
        type->fs = mut_aa[lmut_aa-1] == C4_Stop ? lmut_aa : -1;
    }

    if ( lori_aa>0 || ori_aa ) free(ori_aa);
    if ( lmut_aa>0 || mut_aa ) free(mut_aa);
    return 0;
}

static int predict_molecular_consequence_insertion(struct mc_handler *h, struct mc *mc, struct mc_type *type, struct mc_inf *inf, struct gea_record *v, int lalt, char *alt, int lori, int lmut, char *ori, char *mut)
{

    /* 
       TODO: improve this function!!!!
    */

    
    int cod = (inf->loc-1)%3;
    // re-alignment of reference and alternative alleles
    // inf->loc will be updated in this function
    if ( trim_capped_sequences_and_check_mutated_end(h, mc, type, inf, v, 0, lalt, NULL, alt, &lori, &ori, &lmut, &mut, &cod, 1) == 0 )
        return 0;

    inf->end_loc = inf->loc+1;
    //type->loc_amino = (inf->loc+2)/3;
    type->loc_amino = inf->loc/3+1;
    //type->loc_end_amino = type->loc_amino+1; 
    
    int *ori_aa, *mut_aa;
    int lori_aa, lmut_aa;
    int max_len = v->c.cds_length - v->c.utr5_length - inf->loc;
    ori_aa = convert_bases_to_amino_acid(lori, ori, &lori_aa, mc->is_mito, max_len);
    mut_aa = convert_bases_to_amino_acid(lmut, mut, &lmut_aa, mc->is_mito, -1);

    if ( lori_aa == 0 ) {
        return 0;
    }

    if ( lalt%3 ) { // frameshift
        type->con1 = mc_frameshift_elongation;
        
        int head, tail;
        trim_amino_acid_ends(&lori_aa, ori_aa, &lmut_aa, mut_aa, &head, &tail, 1);
        if ( head ) type->loc_amino += head; // if trimmed head, update new location
        // no end location for frameshift
        type->loc_end_amino = -1;        
        type->ori_amino = ori_aa[0];
        // incase deletion of mutated allele
        if ( lmut_aa <= 1 ) goto insert_a_stop_gained;
        
        type->mut_amino = mut_aa[0];        
        type->fs = mut_aa[lmut_aa-1] == C4_Stop ? lmut_aa : -1;
        
    }
    else { // inframe Insertion
        // the start of insertions is the affected amino acids, no capped base consider
        if ( cod == 0 ) { // p.1_2insXXX
            // related amino acid frames
            //int l = (cod + lalt+2)/3; // 1 based
            type->con1 = mc_conservation_inframe_insertion;

            int head, tail;
            trim_amino_acid_ends(&lori_aa, ori_aa, &lmut_aa, mut_aa, &head, &tail, 1);

            // todo: check duplicate

            if (head) type->loc_amino += head;           
            type->ori_amino = ori_aa[0];
            type->loc_end_amino = type->loc_amino +1;
            type->n = lalt/3;
            if ( lori_aa < type->n ) {
                type->n = lori_aa;
                if ( type->n <= 1 ) goto insert_a_stop_gained;
            }
            type->aminos = malloc(sizeof(int)*type->n);
            
            type->ori_end_amino = -1;            
            int i;
            for ( i=0; i<type->n; ++i) {                
                type->aminos[i] = mut_aa[i];
                if ( type->aminos[i] == 0 ) {
                    //type->con1 = mc_frameshift_truncate; // insert a stop base, transcript truncated
                    if ( i == 0 ) goto insert_a_stop_gained;
                    type->con1 = mc_frameshift_truncate;
                    type->fs = i+1;
                    type->mut_amino = mut_aa[0];
                    type->n = i+1;
                    break;
                }
            }
        }
        else {
            type->con1 = mc_disruption_inframe_insertion;
            int head, tail;
            trim_amino_acid_ends(&lori_aa, ori_aa, &lmut_aa, mut_aa, &head, &tail, 1);
            type->ori_amino = ori_aa[0];  
            type->n = lalt/3; // disrupted amino acid should be emited
            if ( lmut_aa <= type->n ) {
                type->n = lmut_aa;
                if ( type->n <= 1 ) goto insert_a_stop_gained;
            }

            if ( head > 0 ) { // inframe insert, trim start
                type->loc_amino += head-1;
                type->loc_end_amino = type->loc_amino+1;
                type->aminos = malloc(sizeof(int)*type->n);
                int i;
                for ( i=0; i < type->n; ++i ) type->aminos[i] = mut_aa[i];
            }
            else if (mut_aa[lalt/3] == ori_aa[0] ) { // inframe insert, trim end
                type->loc_end_amino = type->loc_amino +1;
                type->aminos = malloc(sizeof(int)*type->n);
                int i;
                for ( i=0; i < type->n; ++i ) type->aminos[i] = mut_aa[i];                
            }
            else { //del-insert                
                // assume delins first
                type->n = lalt/3+1;

                /* if ( lmut_aa < type->n ) { */
                /*     type->n = lmut_aa; */
                /*     if ( type->n <= 1 ) goto insert_a_stop_gained; */
                /* } */
                /* if ( mut_aa[type->n-1] == ori_aa[0]) { // if deleted amino acid is same with inserted tail, convernt to insert */
                /*     type->loc_end_amino = type->loc_amino; */
                /*     type->loc_amino --; */
                /*     type->n--; */
                /* } */
                /* else { */
                type->loc_amino = type->loc_end_amino = (inf->loc+2)/3;
                type->ori_amino = ori_aa[0];                                    
                type->aminos = malloc(sizeof(int)*type->n);
                int i;
                for ( i=0; i < type->n; ++i ) {
                    type->aminos[i] = mut_aa[i];
                    if ( type->aminos[i] == 0 ) {
                        //type->con1 = mc_frameshift_truncate; // insert a stop base, transcript truncated
                        type->con1 = i == 0 ? mc_stop_gained : mc_frameshift_truncate;
                        type->fs = i+1;
                        type->mut_amino = mut_aa[0];
                        break;
                    }
                }
            }
        }
    }

    if ( lori_aa ) free(ori_aa);
    if ( lmut_aa ) free(mut_aa);
    return 0;

  insert_a_stop_gained:
    if ( lori_aa ) free(ori_aa);
    if ( lmut_aa ) free(mut_aa);
    type->con1 = mc_stop_gained;
    type->mut_amino = 0;
    type->n = 0;
    type->fs = -1;
    return 0;
}

/*
 * This function used to parse reference allele and alternative allele sequences. Notice that, all the capped codon will be trimmed.
 * Original sequences generated by bounding the reference allele to disrupted cod   * Mutated sequences reconstruct from original sequences, for frameshift, all remained sequences will be export in this function.
 * For simple case, like SNV, amino acid changes and variant type prediction will be performed in this function; complex case such as
 * frameshift will be consider in compare_reference_and_alternative_amino_acids()
 *
 * return 0 on variant type checked.
 *        1 on non-checked, need do future prediction.
 *       -1 on failure.
 *       
 */
static int compare_reference_and_alternative_allele (struct mc_handler *h,struct mc_inf *inf, struct mc_type *type,struct mc *mc, struct gea_record *v, int lref, char *ref, int lalt, char *alt, int *lori, char **ori, int *lmut, char **mut)
{
    int pos = inf->pos;
    struct gea_coding_transcript *c = &v->c;

//    debug_print("%s\t%d\t%s\t%s\t%s", mc->chr, mc->start, mc->ref, mc->alt, inf->transcript);
    
    // For variants in coding region, check the amino acid changes.
    int cds_pos = pos - c->utr5_length;

    assert(pos > 0 && cds_pos >0);
    // variant start in the amino codon, 0 based, [0,3)
    int cod = (cds_pos-1) % 3;
    
    // codon inition position in the transcript, 0 based.
    // NOTICE: One to 3 base(s) will be CAPPED for the start of insertion; remove caps to check amino acid changes.
    //int start = cod == 0 ? gl->utr5_length + 1 : (cds_pos-1)/3*3 + gl->utr5_length;
    int start = (cds_pos-1)/3*3 + c->utr5_length;
    // retrieve affected sequences and downstream
    char *name = v->name;

    *lori = 0;
    *lmut = 0;
    // original sequnence from RNA sequence, start round from first aa changed position;
    // start aa location may be changed becase realignment, but ori_seq will be "stable" (reset if duplicate) in this function;
    // ori_seq is a temp sequence, will be free before leave this function.
    //
    // update 10th Jun. 2019, retrieve orginal sequence of one exon, not adjust across exons
    *ori = faidx_fetch_seq(h->rna_fai, name, start, start + MAX_TRANSCRIPT_LENGTH, lori);
    
    if ( *ori == NULL || *lori == 0 ) return -1;
    
    // Sometime trancated transcript records will disturb downstream analysis.
    if ( *lori < 3 ) {
        warnings("Transcript record %s probably truncated.", name);
        //if ( *ori ) free(*ori);
        return -1;
    }    

    // number of frames involved in this change, fn_ref for reference frame, fn_alt for alternative frame
    // For example, an insertion of codon between two codons will not change any reference frame, at this time fn_ref == 0 && fn_alt == 1;
    //              a deletion codon between two codons will descrease one frame on the ref, fn_ref == 1 && fn_alt == 0;
    //              a deletion overlapped with two codon will influence two frames, fn_ref == 2 && since new frame will be build fn_alt == 1;
    // int fn_ref, fn_alt;
    construct_alternative_sequence(ref, alt, lref, lalt, cod, *ori, *lori, mut, lmut);

    //static int trim_capped_sequences_and_check_mutated_end(struct mc_handler *h, struct mc *mc, struct mc_type *type, struct mc_inf *inf, struct gea_record *v, int *lori, char **ori, int *lmut, char **mut)
    //if ( mc->type != var_type_snp ) {
        // consider the offset of bases only when insert or delete
        //if ( trim_capped_sequences_and_check_mutated_end(h, mc, type, inf, v, lori, ori, lmut, mut, cod) == 0 )
        //return 0;
    //}

    int ret = -1;
    switch (mc->type) {
        // For SNV, only one codon involved, check the alternative codon and predict the variant type directly.
        case var_type_snp : // fast check SNV, for most cases
            ret = predict_molecular_consequence_snp(mc, type, inf, v, cod, ref, alt, *ori, *mut);
            break;

            // for complex cases, we will build the alternative allele sequence and the amino acids
        case var_type_del:
            ret = predict_molecular_consequence_deletion(h, mc, type, inf, v, lref, ref, *lori, *lmut, *ori, *mut);
            break;
            
        case var_type_ins:
            ret = predict_molecular_consequence_insertion(h, mc, type, inf, v, lalt, alt, *lori, *lmut, *ori, *mut);
            break;
            
        case var_type_delins:
            ret = predict_molecular_consequence_delins(h, mc, type, inf, v, lref, ref, lalt, alt, *lori, *lmut, *ori, *mut);
            break;

        default:
            error("Failed to predict molecular conseqeunce because of: Unknown type %d.", mc->type);
            // break;
    }

    if (ret) return ret;
    // check start and stop codon, update variant type
    if ( type->ori_amino == type->mut_amino) {
        if ( type->ori_amino == C4_Stop ) type->con1 = mc_stop_retained;
        if ( type->loc_amino == 1) type->con1 = mc_start_retained;
    }
    else {
        if ( type->ori_amino == C4_Stop ) type->con1 = mc_stop_loss;
        if ( type->mut_amino == C4_Stop ) type->con1 = mc_stop_gained;
        else if ( type->loc_amino == 1 ) type->con1 = mc_start_loss;
    }
    return 0;
}
static int protein_update_state(struct mc_handler *h, struct mc *mc, struct mc_inf *inf, struct mc_type *type, struct gea_record *v, int lref, int lalt, char *ref, char *alt)
{
    int  lori, lmut;
    char *ori =NULL, *mut = NULL;
    //struct gea_coding_transcript *c = &v->c;
    
    // compare the reference allele and alternative allele sequence, variant position may be update in the function
    if ( compare_reference_and_alternative_allele(h, inf, type, mc, v, lref, ref, lalt, alt, &lori, &ori, &lmut, &mut) == -1 ) {
        warnings("Failed to predict variant type of %s:%d:%s>%s", mc->chr, mc->start, ref, alt);
        if ( lmut > 0 ) free(mut);
        if ( lori > 0 ) free(ori);
        return 1;
    }
    if ( lmut > 0 ) free(mut);
    if ( lori > 0 ) free(ori);
    return 0;
}

// return 1 on whole gene deletion
//        2 on exome deletion
//        0 otherwise
static int whole_gene_deletion_state_update(struct mc *mc, struct gea_record *v)
{
    if ( mc->start <= v->chromStart ) {
        if ( mc->end >= v->chromEnd ) return 1;
        if ( v->blockCount > 1 && mc->end >= v->blockPair[1][0] ) return 2;
    }
    return 0;
}
static int noncoding_transcript_update_molecular_conseqeunce_state(struct mc_handler *h, struct mc *mc, struct mc_inf *inf, struct mc_type *type, struct gea_record *v, char *ref, int lref, char *alt, int lalt)
{
    memset(type, 0, sizeof(struct mc_type));

    if ( type->func1 == func_region_cds ) type->func1 = func_region_noncoding_exon;
    if ( type->func2 == func_region_cds ) type->func2 = func_region_noncoding_exon;

    //int i; // exon / intron ID
    //int pos = inf->pos;
    // set amino acid length to 0
    inf->aa_length = 0;
    
    // exome, just check splice sites
    if ( inf->offset == 0 ) {
        type->con1 = mc_noncoding_exon;
        
        // just check cooridinate is fine
        if ( mc->start < v->blockPair[0][type->count] + splice_site_options.range || mc->end > v->blockPair[1][type->count] - splice_site_options.range)
            type->con1 = mc_noncoding_splice_region;
    }
    // check if variantion located in the splice sites around the edge of utr and cds regions    
    else if ( inf->offset < 0 ) {
        type->con1 = mc_noncoding_intron;
        if ( inf->offset > -splice_site_options.range ) type->con1 = mc_splice_acceptor;
    }
    
    else if ( inf->offset > 0 ) {
        type->con1 = mc_coding_intron;
        if ( inf->offset < splice_site_options.range ) type->con1 = mc_splice_donor;
    }
    // update exon number for minus strand
    
    return 0;
}

static int coding_transcript_update_molecular_consequence_state(struct mc_handler *h, struct mc *mc, struct mc_inf *inf, struct mc_type *type, struct gea_record *v, char *ref, int lref, char *alt, int lalt)
{
    // update the protein changes and variant type state
    if ( protein_update_state(h, mc, inf, type, v, lref, lalt, ref, alt ) ) type->con1 = mc_unknown;

    return 0;
}
static int transcript_function_update(struct mc_handler *h, int is_coding, struct gea_record *v, int *ex, int *pos, int *offset, int position, enum func_region_type *func, enum mol_con *con1, enum mol_con *con_splice, int *loc, int *dist_exon)
{
    if ( find_locate(h->hdr, v, pos, offset, position, ex, dist_exon) ) {
        *func = func_region_outrange;
        return 1; // outside of transcripts ??
    }
    
    // init molecular consequence
    *con1 = mc_unknown;
    *con_splice = mc_unknown; 

    struct gea_coding_transcript *c = &v->c;
    
    if ( is_coding ) {
        if ( *pos <= c->utr5_length ) {
            /*
             * Todo: we only checked the edge of exomes, but for position inside it, insertions and deletions should aslo be consider.
             */
            *loc = c->utr5_length - *pos +1;

            if (*offset) { 
                *con1 = mc_utr5_intron;
                *func = func_region_utr5_intron;
                // check splice sites
                if ( *offset < 0 && *offset > -splice_site_options.range ) *con_splice = mc_splice_acceptor;
                else if ( *offset >0 && *offset < splice_site_options.range ) *con_splice = mc_splice_donor;                
            }
            else {
                
                
                *con1 = mc_utr5_exon;
                *func = func_region_utr5_exon;
            }
        }
        else if ( *pos <= c->cds_length ) {
            *loc = *pos - c->utr5_length;
            if ( *offset ) {
                *con1 = mc_coding_intron;
                *func = func_region_intron;
                if ( *offset < 0 && *offset > -splice_site_options.range ) *con_splice = mc_splice_acceptor;
                else if (*offset >0 && *offset < splice_site_options.range ) *con_splice = mc_splice_donor;                
            }
            else {
                // there is no coding_exon type, to interpret molecular type for coding region thereafter
                *func = func_region_cds;
            }
        }
        else {
            *loc = *pos - c->cds_length;
            if ( *offset ) {
                *con1 = mc_utr3_intron;
                *func = func_region_utr3_intron;
                // check splice sites
                if ( *offset < 0 && *offset > -splice_site_options.range ) *con_splice = mc_splice_acceptor;
                else if (*offset >0 && *offset < splice_site_options.range ) *con_splice = mc_splice_donor;                
            }
            else {
                *con1 = mc_utr3_exon;
                *func = func_region_utr3_exon;
            }
        }
    }
    else {
        *loc = *pos;
        if ( *offset ) {
            *func = func_region_noncoding_intron;
            *con1 = mc_noncoding_intron;
            if ( *offset < 0 && *offset > -splice_site_options.range ) *con_splice = mc_splice_acceptor;
            else if (*offset >0 && *offset < splice_site_options.range ) *con_splice = mc_splice_donor;                
        }
        else {
            *con1 = mc_noncoding_exon;
            *func = func_region_noncoding_exon;
        }
    }

    // check splice sites in exon
    if ( *offset == 0 ) {
        if ( position < v->blockPair[0][*ex]+splice_site_options.range || position>v->blockPair[1][*ex]-splice_site_options.range) *con_splice = mc_exon_splice_sites;
    }

    return 0;
}
//static int mc_anno_trans_core(struct mc_handler *h, struct mc *n, struct mc_core *trans, struct gea_record *v)
static int transcript_molecular_consequence_update(struct mc_handler *h, struct mc *n, struct mc_core *trans, struct gea_record *v)
{
    //assert(n->start >= v->chromStart && n->end <= v->chromEnd);
    // update hgvs locations
    struct mc_inf  *inf = &trans->inf;
    struct mc_type *type = &trans->type;
    const char     *bt = h->hdr->id[GEA_DT_BIOTYPE][v->biotype].key;
    
    assert ( strcmp(bt, "mRNA") == 0 || strcmp(bt, "ncRNA") == 0);

    
    inf->strand = v->strand == strand_is_plus ? '+' : '-';
    inf->transcript = safe_duplicate_string(v->name);
    inf->gene = safe_duplicate_string(v->geneName);


    /*
      Normlise the variant position, some insertion or deletion happened at the tandom repeat regions and cover exon-intron edge, 
      here trying to normlize the inserted or deleted bases to the last position consider the strand.
    */
    // for - strand, already normlized by variantion caller, we just take care of the + strand
        
    // init reference sequences and alternative sequences
    int ref_length = n->ref == NULL ? 0    : strlen(n->ref);
    int alt_length = n->alt == NULL ? 0    : strlen(n->alt);
    char *ref_seq  = n->ref == NULL ? NULL : safe_duplicate_string(n->ref);
    char *alt_seq  = n->alt == NULL ? NULL : safe_duplicate_string(n->alt);
    // ref_seq and alt_seq will be reverse if - strand
    // char *ref, *alt, start, end;
    int start = n->start;
    int end = n->end;
    if (n->start != n->end) { // indel, try to normlise
        if (inf->strand == '+') {
            faidx_t *fai = fai_load(h->reference_fname);
            if (fai) {
                char *seq;
                int l;
                seq = faidx_fetch_seq(fai, n->chr, start, MAX_NORMLIZE_LENGTH, &l);
                if (l < MAX_NORMLIZE_LENGTH) {// too short, escape normlise
                    free(seq);
                }
                else {
                    kstring_t alt = {0,0,0};
                    if (alt_length) kputs(alt_seq, &alt);
                    kputs(seq+ref_length, &alt);
                    int i;
                    for (i = 0; i < l && l < alt.l; i++)
                        if (seq[i] != alt.s[i]) break;
                    if (i) {
                        if (ref_length) memcpy(ref_seq, seq+i, ref_length *sizeof(char));
                        if (alt_length) memcpy(alt_seq, alt.s+i, alt_length*sizeof(char));
                        start+=i;
                        end += i;
                    }
                    free(alt.s);
                }
                fai_destroy(fai);
            }
        }
    }
    // debug_print("%s\t%s\t%d\t%d", ref_seq, alt_seq, start, end);
    //debug_print("%s\t%d\t%s\t%s\t%s", n->chr, n->start, n->ref, n->alt, inf->transcript);
     
    int ex1 = 0, ex2 = 0;
    
    // update Amino acid length
    /* if ( strcmp(bt, "mRNA") == 0 ) inf->aa_length = (c->cds_length - c->utr5_length)/3; */
    /* else inf->aa_length = 0;       */

    // no matter whole gene/exome deletion, interupt as exon_loss
    // Update: try to generate HGVS names for exon loss...
    if ( whole_gene_deletion_state_update(n, v) )  type->con1 = mc_exon_loss;
    //{ type->con1 = mc_exon_loss; return 0; }    
    
    int is_coding = 0;
    if ( strcmp(bt, "mRNA") == 0 ) is_coding = 1;

    
    // set amino acid length
    if (is_coding) {
        inf->aa_length = (v->c.cds_length -v->c.utr5_length)/3-1; // update: 2019/06/29 emit stop codon
    }

    transcript_function_update(h, is_coding, v, &ex1, &inf->pos, &inf->offset, start, &type->func1, &type->con1, &type->con2, &inf->loc, &inf->dist_exon);
    
    type->count = ex1; // For now, count may NOT be the real exome number, considering the backward strand.
    
    if (start != end ) {
        // unannotated location will export as '?'
        enum mol_con con1 = mc_unknown;
        enum mol_con con2 = mc_unknown;
        
        transcript_function_update(h, is_coding, v, &ex2, &inf->end_pos, &inf->end_offset, end, &type->func2, &con1, &con2, &inf->end_loc, &inf->dist_exon);
        
        if ( v->strand == strand_is_minus ) {
            int t = inf->pos;
            inf->pos = inf->end_pos;
            inf->end_pos = t;
            t = inf->offset;
            inf->offset = inf->end_offset;
            inf->end_offset = t;
            t = inf->loc;
            inf->loc = inf->end_loc;
            inf->end_loc = t;
            t = type->func2;
            type->func2 = type->func1;
            type->func1 = t;
            uint16_t x = (uint16_t)inf->dist_exon;
            inf->dist_exon = ((int32_t)x)<<16|(uint16_t)(inf->dist_exon>>16);
        }
        if ( con2 != mc_unknown ) type->con2 = con2;
        if ( type->func1 != type->func2 && n->type != var_type_ins) { type->con1 = mc_exon_loss; return 0;}
    } 
    else {
        inf->end_pos = inf->pos;        
        inf->end_offset = inf->offset;
        inf->end_loc = inf->loc;
        type->func2 = type->func1;
    }

    if ( type->con1 != mc_unknown ) {
        if ( v->strand == strand_is_minus ) {        
            if (ref_length) {
                inf->ref = safe_duplicate_string(ref_seq);
                compl_seq(inf->ref, strlen(inf->ref));
            }
            if (alt_length) {
                inf->alt = safe_duplicate_string(alt_seq);
                compl_seq(inf->alt, strlen(inf->alt));
            } 
        }

        return 0;
    }
    
    // For variant in coding region, check amino acid(s) changes.
    
    // for reverse strand, complent sequence
    if ( inf->strand == '-' ) {        
        if ( ref_seq ) compl_seq(ref_seq, ref_length);
        if ( alt_seq ) compl_seq(alt_seq, alt_length);
    }
    
    if ( strcmp(bt, "mRNA") == 0 )
        coding_transcript_update_molecular_consequence_state (h, n, inf, type, v, ref_seq, ref_length, alt_seq, alt_length);
    else 
        noncoding_transcript_update_molecular_conseqeunce_state (h, n, inf, type, v, ref_seq, ref_length, alt_seq, alt_length);
    
    if ( ref_seq != NULL ) free(ref_seq);
    if ( alt_seq != NULL ) free(alt_seq);

    // update variant functional types
    // check_func_vartype(h, n, n->n_tran, v);
    return 0;
}

static int up_downstream_gene_update(struct mc *n, struct gea_record *last, struct gea_record *next)
{
    struct intergenic_core *inter = &n->inter;
    
    if ( last != NULL && next != NULL ) {
        int mid = (next->chromStart + last->chromEnd)/2;
        if ( n->start < mid ) goto near_last_gene;
        else goto near_next_gene;
    }

    if ( last != NULL && next == NULL ) goto near_last_gene;

    if ( last == NULL && next != NULL ) goto near_next_gene;

    return 1;

  near_last_gene:
    inter->gene = last->geneName == NULL ? NULL : safe_duplicate_string(last->geneName);
    if ( last->strand == strand_is_plus ) {
        // forward strand, variant located in the downstream of gene, gene length should be added for TSS 
        inter->TSS_dist = n->start - last->chromStart;
        inter->con1 = inter->TSS_dist <= 1000 ? mc_downstream_1KB : mc_downstream_10KB;
    }
    else {
        // backward strand, var in upstream and treat chromEnd as TSS
        inter->TSS_dist = n->start - last->chromEnd;
        inter->con1 = inter->TSS_dist <= 1000 ? mc_upstream_1KB : mc_upstream_10KB ;
    }
    return 0;

  near_next_gene:
    inter->gene = next->geneName == NULL ? NULL : safe_duplicate_string(next->geneName);
    if ( next->strand == strand_is_plus ) {
        // forward strand, variant located in the downstream of gene, gene length should be added for TSS 
        inter->TSS_dist = n->start - next->chromStart;
        inter->con1 = inter->TSS_dist <= 1000 ? mc_upstream_1KB : mc_upstream_10KB;
    }
    else {
        // backward strand, var in upstream and treat chromEnd as TSS
        inter->TSS_dist = n->start - next->chromEnd;
        inter->con1 = inter->TSS_dist <= 1000 ? mc_downstream_1KB : mc_downstream_10KB;
    }
    return 0;
}

// return current iterate in the buffer
//
// check the chunk covered this variant
// R1    |===================|
// R2     |====|
// R3           |=====|
// R4                   |======|
// R5                              |====|
//               1
//
//                           2
//                               3
//                                         4
//
// this function will give the interval of chunck covered variant
static int count_all_overlapped_records(struct mc_handler *h, struct mc *n, int *c, int *a, int *intragenic_flag)
{
    int i; // used for iter records only, avoid to use it again
    *c = *a = *intragenic_flag = 0;
    
    // check how mang genes and/or motifs overlaped with this variant, it is the first round to check the records, we will check twice
    for ( i = h->i_record; i < h->n_record; ++i ) {
    //for ( i = 0; i < h->n_record; ++i ) {
        struct gea_record *v = h->records[i];

        // update end_pos_for_skip marker
        if ( h->end_pos_for_skip == 0 ) h->end_pos_for_skip = v->chromEnd;
        //if ( n->start > h->end_pos_for_skip) {// if ( h->end_pos_for_skip < v->chromEnd ) {
        //}
        if ( h->end_pos_for_skip < v->chromEnd && n->start > h->end_pos_for_skip ) {
            int j;
            for ( j = h->i_record; j < i; ++j ) {
                struct gea_record *g = h->records[j];
                if ( g->chromEnd > n->start ) {
                    h->i_record = j;
                    h->end_pos_for_skip = g->chromEnd;
                    break;
                }
            }
            if ( j == i ) {
                h->i_record = i;
                h->end_pos_for_skip = v->chromEnd;
            }
        }

        const char *bt = h->hdr->id[GEA_DT_BIOTYPE][v->biotype].key;
        // only consider coding transcript
        // Since upstream gene records had be filled, and all records in the chunk keep in coordinate,
        // last_gene will alway point to the last gene record or be NULL.
        if ( (strcmp(bt, "mRNA") == 0 || strcmp(bt, "Gene") == 0 ) && v->chromEnd < n->start) h->last_gene = v;

        // check if variant located in this record
        if ( n->start > v->chromEnd ) continue; 
        if ( n->end < v->chromStart ) break;
        
        if ( strcmp(bt, "mRNA") == 0 || strcmp(bt, "ncRNA") == 0 ) (*a)++;

        // if no Gene record in database, skip to check intragenic variant
        if ( strcmp(bt, "Gene") == 0 ) *intragenic_flag = 1;
        // count all records
        (*c)++;

    }
    return i;
}

// 0 on sucess, 1 on failure
static int next_gene_record_update(struct mc_handler *h, int i)
{
    int j;
    if ( i < h->n_record ) {
        for (j = i; j < h->n_record; ++j) {
            struct gea_record *v = (struct gea_record*)h->records[j];
            char *bt = (char*)h->hdr->id[GEA_DT_BIOTYPE][v->biotype].key;
            if ( strcmp(bt, "mRNA") && strcmp(bt, "ncRNA") && strcmp(bt, "Gene") ) continue;
            h->next_gene = v;
            return 0;
        }
    }
    return 1;
}
// 0 on intragenic, 1 on otherwise
static int intragenic_variant_state_update(struct mc *n, int intragenic_flag)
{
    struct intergenic_core *inter = &n->inter;
    // if intragenic flag is set and no transcript overlapped with it, interupt it as intragenic variant
    if ( intragenic_flag ) {
        // c should always greater than 0
        inter->con1 = mc_intragenic;
        return 0;
    }
    return 1;
}
static int tfbs_variant_state_update(struct mc *n, struct mc_handler *h, struct gea_record *v)
{
    struct intergenic_core *inter = &n->inter;
    char *bt = (char*)h->hdr->id[GEA_DT_BIOTYPE][v->biotype].key;
    //inter->con1 =  inter->con2 = mc_unknown;
    if ( strcmp(bt, "TFBS") == 0 ) { 
        if (n->start <= v->chromStart && n->end >= v->chromEnd ) inter->con1 = mc_tfbs_ablation;
        else inter->con1 = mc_tfbs_variant;
        return 0; 
    }
    return 1;
}
// c is the total count of records, include motifs and genes
static int intergenic_variant_state_update(struct mc *n, struct mc_handler *h, int c)
{
    struct intergenic_core *inter = &n->inter;
    //int upstream_gap = 0, downstream_gap = 0;
    // init molecular consequence 
    inter->con1 = inter->con2 = mc_intergenic;
       
    // motifs, or gene record for intragenic variant
    // the rules to convert inter->con1 and inter->con2 into molecular consequence string:
    // if con1 != intergenic { print con1; if con2 != intergenic; print +con2; }
    // else print con2; 
    // if c== 0, the function should never work            
    int j;
    for ( j = h->i_record; j < c; ++j ) { // check motif and gene record
        struct gea_record *v = (struct gea_record*)h->records[j];
        if ( tfbs_variant_state_update(n, h, v) == 0 ) return 0;
    }
    return 1;
}
//
static int transcripts_variant_state_update(struct mc *n, struct mc_handler *h, int a, int c)    
{
    // if overlapped with gene, interpret the amino acid changes and variant types
    // for variant in the intron region, motif should also be checked.
    int j, i;
    int tfbs_region_skip = 0;
    n->trans = malloc(a*sizeof(struct mc_core));    
    n->n_tran = 0;    
    
    for ( i = h->i_record, j = 0; i < h->n_record && j < c; ++i,++j ) {
        struct gea_record *v = (struct gea_record*)h->records[i];
        const char *bt = h->hdr->id[GEA_DT_BIOTYPE][v->biotype].key;        
        if ( n->end < v->chromStart || n->start > v->chromEnd) continue;
        // for this part we only consider transcripts
        // intron region also could be overlapped with motifs
        if ( strcmp(bt, "TFBS") == 0 ) {
            if ( tfbs_region_skip == 0 ) {
                if ( tfbs_variant_state_update(n, h, v) == 0 ) tfbs_region_skip = 1;
            }
            continue;
        }
        if ( strcmp(bt, "mRNA") && strcmp(bt, "ncRNA") ) continue;
        
        // clean core structure for updating transcript record
        struct mc_core *trans = &n->trans[n->n_tran];
        memset(trans, 0, sizeof(*trans));
        //int ret;
        if ( transcript_molecular_consequence_update(h, n, trans, v) ) continue;
        if ( v->strand == strand_is_plus ) trans->type.count++;
        else {
            trans->type.count = v->blockCount-trans->type.count;
            // if intron, --count
            if (trans->inf.offset) trans->type.count--;
        }
        n->n_tran++;
    }
    return n->n_tran;
}
static int intergenic_update_molecular_consequence_state(struct mc *n, struct mc_handler *h, int *a, int *c)
{
    // used to check if variant located in gene region but not overlapped with any transcript,
    // if it set and no transcript overlapped with this variant, we could assume this is a intragenic variant
    int intragenic_flag;
    
    int i = count_all_overlapped_records(h, n, c, a, &intragenic_flag);
    
    // if no gene overlapped with this variant, interupt as intergenic variant, check the up/downstream gene
    if ( *a == 0 ) { // variants located in intergenic regions

        // update general intergenic information
        intergenic_variant_state_update(n, h, *c);
        
        // intragenic variant
        if ( intragenic_variant_state_update(n, intragenic_flag) == 0 ) return 0;

        // find the upstream nearest gene if possibe
        next_gene_record_update(h, i);

        // check the gaps between variant and nearby genes
        up_downstream_gene_update(n, h->last_gene, h->next_gene);
        
        return 0;
    }
    // return records should be count
    return i;
}
// There is no need to consider the difference of motif records overlapped with the variant, because we only need to interpet
// it is in motif region.
int mc_anno_trans_chunk(struct mc *mc, struct mc_handler *h)
{
    int c; // iterater for total of motif and genes
    int a; // iterater for transcript records only
    int n;
    n = intergenic_update_molecular_consequence_state(mc, h, &a, &c);

    if ( n == 0 ) return -1;
    
    // ** transcripts **
    return transcripts_variant_state_update(mc, h, a, n);
}

// todo: export details of a variant
/* static int variants_describe_details(struct mc_handler *h, char *chr, int pos, char *ref, char *alt) */
/* { */
/*     int lref = strlen(ref); */
/*     int lalt = strlen(alt); */
    
/*     struct mc *mc = mc_init(chr, pos, pos+lref, ref, alt); */

/*     if ( mc_anno_trans_chunk(mc, h) == -1) return -1; */
    
/* } */
/*
  Predefined Tags :
  Gene, Transcript, HGVSnom, ExonIntron, VarType, MC, TFBS, 
  ##INFO=<ID=Gene,Number=A,Type=String,Description="Gene name(s).">
  ##INFO=<ID=Transcript,Number=A,Type=String,Description="Transcript name(s).">
  ##INFO=<ID=HGVSnom,Number=A,Type=String,Description="HGVS nomenclature for the description of DNA sequence variants">
  ##INFO=<ID=VarType,Number=A,Type=String,Description="Variant type."> 
  ##INFO=<ID=MC,Number=A,Type=String,Description="Molecular consequence, similar with VarType but give more details."> 
  ##INFO=<ID=ExonIntron,Number=A,Type=String,Description="Exon/CDS or intron id on transcripts.">
  ##INFO=<ID=AAlength,Number=A,Type=String,Description="Amino acid length for each transcript. 0 for noncoding transcript.">"
  ##INFO=<ID=TFBS,Number=1,Type=Flag,Description="Set if the variant cover or overlap transcript factor binding sites.">
  ##INFO=<ID=DistanceGene,Number=1,Type=Integer,Description="Distance to nearby gene, this value could be positive for upstream and negative for downstream.">

Supported for special test:
  ##INFO=<ID=ANNOVARname,Number=A,Type=String,Description="Variant description in ANNOVAR format.">
  ##INFO=<ID=IVSnom,Number=A,Type=String,Description="Old style nomenclature for the description of intron variants. Not recommand to use it.">
  ##INFO=<ID=Oldnom,Number=A,Type=String,Description="Old style nomenclature, compared with HGVSnom use gene position instead of UTR/coding position.">
 */

struct anno_mc {
    htsFile *fp;
    faidx_t *rna_fai;
    tbx_t   *idx;
    void    *name_hash;
    int      n_col;
    struct anno_col *cols;
};

static int anno_mc_update_buffer(struct anno_mc_file *file, bcf_hdr_t *hdr, bcf1_t *line)
{
    bcf_unpack(line, BCF_UN_INFO);
    
    int i;
    // clean buffer
    for ( i = 0; i < file->n_allele; ++i ) mc_destroy(file->files[i]);
    file->n_allele = line->n_allele-1; // emit ref
    file->files = malloc(file->n_allele*sizeof(struct mc));
    for ( i = 0; i < file->n_allele; ++i )
        file->files[i] = mc_init(bcf_seqname(hdr, line), line->pos+1, line->pos+line->rlen, line->d.allele[0], line->d.allele[i+1]);
    return 0;
}
static int empty_tag_string(kstring_t *str)
{
    int i;
    for (i=0; i<str->l; ++i ) 
        if (str->s[i] != '|' && str->s[i] != ',' && str->s[i] != '.') return 0;
    return 1;
}
static void generate_hgvsnom_string_empty(kstring_t *str)
{
    kputc('.', str);    
}
static void generate_hgvsnom_string_NoncodingExon(struct mc const *mc, struct mc_type const *type, struct mc_inf const *inf, kstring_t *str)
{
    ksprintf(str, "%s:n.%d", inf->transcript, inf->loc);
    //if ( inf->offset > 0 ) ksprintf(str, "+%d", inf->offset);
    //else if ( inf->offset < 0 ) ksprintf(str, "%d", inf->offset);
    char *ref = inf->ref != NULL ? inf->ref : mc->ref;
    char *alt = inf->alt != NULL ? inf->alt : mc->alt;
        
    if ( inf->end_loc != inf->loc ) {
        ksprintf(str, "_%d", inf->end_loc);
        //if ( inf->end_offset > 0 ) ksprintf(str, "+%d", inf->end_offset);
        //else if ( inf->end_offset < 0 ) ksprintf(str, "%d", inf->end_offset);
        if ( ref == NULL ) { // insertion
            ksprintf(str, "ins%s", alt);
        }        
        else {
            if ( alt == NULL ) { // deletion
                ksprintf(str, "del%s", ref);
            }
            else { // del-ins
                ksprintf(str, "delins%s", alt);
            }
        }
    }
    else {
        //if ( ref == NULL ) error("For developer: bad range for noncoding insertion. %s:%d", mc->chr, mc->start);
        if ( ref == NULL ) {
            warnings("Bad range for insertion. %s:%d, %s.", mc->chr, mc->start, inf->transcript);
            kputc('?', str);
        }
        else if ( alt == NULL ) ksprintf(str, "del%s", ref);
        else {
            if ( strlen(alt) == 1 ) ksprintf(str, "%s>%s", ref, alt);
            else ksprintf(str, "del%sins%s", ref, alt);
        }
    }
}
static void generate_hgvsnom_string_intron2(struct mc const *mc, struct mc_type const *type, struct mc_inf const *inf, kstring_t *str)
{
    ksprintf(str, "%s:", inf->transcript);
    /*
    if ( type->func1 == func_region_utr5_intron ) ksprintf(str, "c.-%d", inf->loc);
    else if ( type->func1 == func_region_utr3_intron ) ksprintf(str, "c.*%d", inf->loc);
    else if ( type->func1 == func_region_intron ) ksprintf(str, "c.%d", inf->loc);
    else if ( type->func1 == func_region_noncoding_intron) ksprintf(str, "n.%d", inf->loc);
    else error("For developer: only intron region should come here!");
    */
    ksprintf(str, "IVS%d", type->count);
    if ( inf->offset > 0 ) ksprintf(str, "+%d", inf->offset);
    else if ( inf->offset < 0 ) ksprintf(str, "%d", inf->offset);
    else error("For developer: offset should NOT be 0 for intron.");
 
    char *ref = inf->ref != NULL ? inf->ref : mc->ref;
    char *alt = inf->alt != NULL ? inf->alt : mc->alt;
    
    if ( inf->end_offset != inf->offset ) { // since this function handler all intron variants, we only need to check the offsets
        kputc('_', str);
        /*
        if ( type->func1 == func_region_utr5_intron ) ksprintf(str, "-%d", inf->end_loc);
        else if ( type->func1 == func_region_utr3_intron ) ksprintf(str, "*%d", inf->end_loc);
        else if ( type->func1 == func_region_intron ) ksprintf(str, "%d", inf->end_loc);
        else if ( type->func1 == func_region_noncoding_intron) ksprintf(str, "%d", inf->end_loc);
        else if ( type->func1 == func_region_cds ) ksprintf(str, "%d", inf->end_loc);
        */
        if ( inf->end_offset > 0 ) ksprintf(str, "+%d", inf->end_offset);
        else if ( inf->end_offset < 0 ) ksprintf(str, "%d", inf->end_offset);
        
        if ( ref == NULL ) { // insertion
            ksprintf(str, "ins%s", alt);
        }        
        else {
            if ( alt == NULL ) { // deletion
                ksprintf(str, "del%s", ref);
            }
            else { // del-ins
                ksprintf(str, "delins%s", alt);
            }
        }
    }
    else {
        if (ref == NULL ) error("For developer: bad range for noncoding insertion. %s:%d", mc->chr, mc->start);
        if (alt == NULL ) ksprintf(str,"del%s",ref);
        else {
            if ( strlen(alt) == 1 ) ksprintf(str, "%s>%s", ref, alt);
            else ksprintf(str, "del%sins%s", ref, alt);
        }
    }
}
static void generate_hgvsnom_string_intron(struct mc const *mc, struct mc_type const *type, struct mc_inf const *inf, kstring_t *str)
{
    ksprintf(str, "%s:", inf->transcript);
    if ( type->func1 == func_region_utr5_intron ) ksprintf(str, "c.-%d", inf->loc);
    else if ( type->func1 == func_region_utr3_intron ) ksprintf(str, "c.*%d", inf->loc);
    else if ( type->func1 == func_region_intron ) ksprintf(str, "c.%d", inf->loc);
    else if ( type->func1 == func_region_noncoding_intron) ksprintf(str, "n.%d", inf->loc);
    else if ( type->func1 == func_region_cds) ksprintf(str, "c.%d", inf->loc);
    else error("For developer: only intron region should come here! %s:%d", mc->chr, mc->start);
    
    if ( inf->offset > 0 ) ksprintf(str, "+%d", inf->offset);
    else if ( inf->offset < 0 ) ksprintf(str, "%d", inf->offset);
    // else error("For developer: offset should NOT be 0 for intron."); // delete because insert at splice position
 
    char *ref = inf->ref != NULL ? inf->ref : mc->ref;
    char *alt = inf->alt != NULL ? inf->alt : mc->alt;
    
    if ( inf->end_offset != inf->offset ) { // since this function handler all intron variants, we only need to check the offsets
        kputc('_', str);
        if ( type->func1 == func_region_utr5_intron ) ksprintf(str, "-%d", inf->end_loc);
        else if ( type->func1 == func_region_utr3_intron ) ksprintf(str, "*%d", inf->end_loc);
        else if ( type->func1 == func_region_intron ) ksprintf(str, "%d", inf->end_loc);
        else if ( type->func1 == func_region_noncoding_intron) ksprintf(str, "%d", inf->end_loc);
        else if ( type->func1 == func_region_cds ) ksprintf(str, "%d", inf->end_loc);

        if ( inf->end_offset > 0 ) ksprintf(str, "+%d", inf->end_offset);
        else if ( inf->end_offset < 0 ) ksprintf(str, "%d", inf->end_offset);
        
        if ( ref == NULL ) { // insertion
            ksprintf(str, "ins%s", alt);
        }        
        else {
            if ( alt == NULL ) { // deletion
                ksprintf(str, "del%s", ref);
            }
            else { // del-ins
                ksprintf(str, "delins%s", alt);
            }
        }
    }
    else {
        if (ref == NULL ) error("For developer: bad range for noncoding insertion. %s:%d", mc->chr, mc->start);
        if (alt == NULL ) ksprintf(str,"del%s",ref);
        else {
            if ( strlen(alt) == 1 ) ksprintf(str, "%s>%s", ref, alt);
            else ksprintf(str, "del%sins%s", ref, alt);
        }
    }
}

static void generate_hgvsnom_string_UtrExon(struct mc const *mc, struct mc_type const *type, struct mc_inf const *inf, kstring_t *str)
{
    ksprintf(str, "%s:", inf->transcript);
    if ( type->func1 == func_region_utr5_exon ) ksprintf(str, "c.-%d", inf->loc);
    else if ( type->func1 == func_region_utr3_exon ) ksprintf(str, "c.*%d", inf->loc);
    //else if ( type->func1 == func_region_intron ) ksprintf(str, "c.%d", inf->loc);
    else kputs("c.?", str);
        //error("For developer: only UTR region should come here!");

    char *ref = inf->ref == NULL ? mc->ref : inf->ref;
    char *alt = inf->alt == NULL ? mc->alt : inf->alt;
    
    if ( inf->end_loc != inf->loc ) {
        kputc('_', str);

        if ( type->func1 == func_region_utr5_exon ) ksprintf(str, "-%d", inf->end_loc);
        else if ( type->func1 == func_region_utr3_exon ) ksprintf(str, "*%d", inf->end_loc);
        else if ( type->func1 == func_region_intron ) ksprintf(str, "%d", inf->end_loc);
        else if ( type->func1 == func_region_cds ) ksprintf(str, "%d", inf->end_loc);
        else kputc('?', str);
        // end position could cover intron region
        if ( inf->end_offset > 0 ) ksprintf(str, "+%d", inf->end_offset);
        else if ( inf->end_offset < 0 ) ksprintf(str, "%d", inf->end_offset);

        if ( mc->ref == NULL ) { // insertion
            ksprintf(str, "ins%s", alt);
        }        
        else {
            if ( mc->alt == NULL ) { // deletion
                ksprintf(str, "del%s", ref);
            }
            else { // del-ins
                ksprintf(str, "delins%s", alt);
            }
        }
    }
    else {
        if ( ref == NULL ) {
            warnings("Bad range for insertion. %s:%d", mc->chr, mc->start);
            kputc('?', str);
        }
        else if ( alt == NULL ) ksprintf(str, "del%s", ref);
        else {
            if ( strlen(alt) == 1 ) ksprintf(str, "%s>%s", ref, alt);
            else ksprintf(str, "del%sins%s", ref, alt);
        }
    }
}
static void generate_hgvsnom_string_CodingDelins(struct mc const *mc, struct mc_type const *type, struct mc_inf const *inf, kstring_t *str)
{

    char *ref = inf->ref != NULL ? inf->ref : mc->ref;
    char *alt = inf->alt != NULL ? inf->alt : mc->alt;

    ksprintf(str, "%s:c.%d", inf->transcript, inf->loc);
    
    if ( inf->end_loc != inf->loc ) {
        kputc('_', str);
        if ( type->func1 == func_region_cds ) ksprintf(str, "%d", inf->end_loc);
        else if ( type->func1 == func_region_utr3_exon ) ksprintf(str, "*%d", inf->end_loc);
        else if ( type->func1 == func_region_intron ) {
            ksprintf(str, "%d", inf->end_loc);
            // end position could cover intron region
            if ( inf->end_offset > 0 ) ksprintf(str, "+%d", inf->end_offset);
            else if ( inf->end_offset < 0 ) ksprintf(str, "%d", inf->end_offset);
        }
        
        if ( inf->ref == NULL ) { // insertion
            ksprintf(str, "ins%s", alt);
        }        
        else {
            if ( inf->alt == NULL ) { // deletion
                // ksprintf(str, "del%s",ref);
                ksprintf(str, "del");
            }
            else { // del-ins
                ksprintf(str, "delins%s", alt);
            }
        }
    }
    else {
        if ( ref == NULL ) {
            error("For developer: bad range for insertion. %s:%d", mc->chr, mc->start);
        }
        else if ( alt == NULL ) kputs("del", str); //ksprintf(str, "del%s",ref);
        else {
            if ( strlen(alt) == 1 && strlen(ref) == 1 ) ksprintf(str, "%s>%s", ref, alt);
            else ksprintf(str, "del%sins%s", ref, alt);
        }
    }

    if ( inf->inframe_stop == 1 ) {
        ksprintf(str, "(?)");
    }
    // amino acid changes
    else if ( type->con1 == mc_frameshift_truncate || type->con1 == mc_frameshift_elongation )
        
        if ( type->fs > 0 ) {
            if ( type->ori_amino == C4_Stop ) { // ext
                ksprintf(str, "(p.%s%d%sext*%d/p.%s%d%sext*%d)",
                         codon_names[type->ori_amino],
                         type->loc_amino,
                         codon_names[type->mut_amino],
                         type->fs-1,
                         codon_short_names[type->ori_amino], type->loc_amino, codon_short_names[type->mut_amino], type->fs-1);                
            }
            else {
                ksprintf(str, "(p.%s%d%sfs*%d/p.%s%d%sfs*%d)",
                         codon_names[type->ori_amino],
                         type->loc_amino,
                         codon_names[type->mut_amino],
                         type->fs,
                         codon_short_names[type->ori_amino], type->loc_amino, codon_short_names[type->mut_amino], type->fs);
            }
        }
        else
            ksprintf(str, "(p.%s%d%s/p.%s%d%s)",
                     codon_names[type->ori_amino],
                     type->loc_amino,
                     codon_names[type->mut_amino],
                     codon_short_names[type->ori_amino], type->loc_amino, codon_short_names[type->mut_amino]);
       
            
    else { // inframe delins
        if ( type->loc_amino == type->loc_end_amino ) {
            if ( type->n == 1 ) {
                ksprintf(str, "(p.%s%d%s/p.%s%d%s)", codon_names[type->ori_amino], type->loc_amino, codon_names[type->mut_amino], codon_short_names[type->ori_amino], type->loc_amino, codon_short_names[type->mut_amino]);
            }
            else if ( type->n ) {
                int i;
                ksprintf(str, "(p.%s%ddelins", codon_names[type->ori_amino], type->loc_amino);
                for ( i = 0; i < type->n; ++i ) kputs(codon_names[type->aminos[i]], str);
                kputc(')', str);
                //ksprintf(str, "(/p.%s%ddelins)", codon_short_names[type->ori_amino], type->loc_amino);
                //for ( i = 0; i < type->n; ++i ) kputs(codon_short_names[type->aminos[i]], str);
            }
            else {
                ksprintf(str, "(p.%s%ddel/p.%s%ddel)", codon_names[type->ori_amino], type->loc_amino, codon_short_names[type->ori_amino], type->loc_amino);
            }
        }
        else {
            if ( type->n ) {
                if (type->con1 == mc_conservation_inframe_insertion ) {
                  inframe_insertion:
                    ksprintf(str, "(p.%d_%dins", type->loc_amino, type->loc_end_amino);
                    int i;
                    for ( i = 0; i < type->n; ++i ) kputs(codon_names[type->aminos[i]], str);
                    kputc('/', str);
                    ksprintf(str, "p.%d_%dins", type->loc_amino, type->loc_end_amino);
                    for ( i = 0; i < type->n; ++i ) kputs(codon_short_names[type->aminos[i]], str);
                    kputc(')', str);
                }
                else if ( type->con1 == mc_disruption_inframe_insertion ) {
                    if ( type->loc_amino == type->loc_end_amino ) {
                        int i;
                        ksprintf(str, "(p.%s%ddelins", codon_names[type->loc_amino], type->loc_amino);
                        for ( i = 0; i < type->n; ++i ) kputs(codon_names[type->aminos[i]], str);
                        kputc('/', str);
                        ksprintf(str, "p.%s%ddelins", codon_short_names[type->loc_amino], type->loc_amino);
                        for ( i = 0; i < type->n; ++i ) kputs(codon_short_names[type->aminos[i]], str);
                        kputc(')', str);
                    }
                    else goto inframe_insertion;                    
                }
                else {
                    // disruption_inframe_deletion
                    int i;
                    ksprintf(str, "(p.%d_%ddelins", type->loc_amino, type->loc_end_amino);
                    for ( i = 0; i < type->n; ++i ) kputs(codon_names[type->aminos[i]], str);
                    ksprintf(str, "/p.%d_%ddelins", type->loc_amino, type->loc_end_amino);
                    for ( i = 0; i < type->n; ++i ) kputs(codon_short_names[type->aminos[i]], str);
                    kputc(')', str);
                }
            }
            else { 
                if (type->con1 == mc_conservation_inframe_deletion ) {
                    ksprintf(str, "(p.%s%d_%s%ddel/p.%s%d_%s%ddel)", codon_names[type->ori_amino], type->loc_amino, codon_names[type->ori_end_amino], type->loc_end_amino,
                             codon_short_names[type->ori_amino], type->loc_amino, codon_short_names[type->ori_end_amino], type->loc_end_amino
                    );
                }
                else if ( type->con1 == mc_disruption_inframe_deletion) {
                    ksprintf(str, "(p.%s%d_%s%ddel/p.%s%d_%s%ddel)", codon_names[type->ori_amino], type->loc_amino, codon_names[type->ori_end_amino], type->loc_end_amino,
                                 codon_short_names[type->ori_amino], type->loc_amino, codon_short_names[type->ori_end_amino], type->loc_end_amino
                            );
                }
                else {
                    ksprintf(str, "(p.%s%d%s/p.%s%d%s)", codon_names[type->ori_amino], type->loc_amino, codon_names[type->mut_amino],codon_short_names[type->ori_amino], type->loc_amino, codon_short_names[type->mut_amino]);
                }
            }
        }
    }
    
    
}
static void generate_hgvsnom_string_CodingSNV(struct mc *mc, struct mc_type *type, struct mc_inf *inf, kstring_t *str)
{
    char *ref = inf->ref != NULL ? inf->ref : mc->ref;
    char *alt = inf->alt != NULL ? inf->alt : mc->alt;

    ksprintf(str, "%s:c.%d%s>%s(p.%s%d%s/p.%s%d%s)", inf->transcript, inf->loc, ref, alt, codon_names[type->ori_amino], type->loc_amino, codon_names[type->mut_amino], codon_short_names[type->ori_amino], type->loc_amino, codon_short_names[type->mut_amino]);
}
static void generate_hgvsnom_string_NoCall(struct mc *mc, struct mc_type *type, struct mc_inf *inf, kstring_t *str)
{
    // debug_print("%s\t%d\t%s\t%s\t%s", mc->chr, mc->start, mc->ref, mc->alt, inf->transcript);
    char *ref = inf->ref != NULL ? inf->ref : mc->ref;
    char *alt = inf->alt != NULL ? inf->alt : mc->alt;
    int lref = ref == NULL ? 0 : strlen(ref);
    int lalt = alt == NULL ? 0 : strlen(alt);
    
    if (type->func1 == func_region_cds ) {
        if (lref == 1 && lalt == 1 ) generate_hgvsnom_string_CodingSNV(mc, type, inf, str);
        else generate_hgvsnom_string_CodingDelins(mc, type, inf, str);        
    }
    else {
        ksprintf(str, "%s:", inf->transcript);
        if ( type->func1 == func_region_utr5_exon ) ksprintf(str, "c.-%d", inf->loc);
        else if ( type->func1 == func_region_utr3_exon ) ksprintf(str, "c.*%d", inf->loc);
        else if ( type->func1 == func_region_noncoding_exon ) ksprintf(str, "n.%d", inf->loc);
        else error("For developer: only UTR region should come here! %d", type->func1);

        if ( inf->end_loc != inf->loc ) {
            kputc('_', str);

            if ( type->func1 == func_region_utr5_exon ) ksprintf(str, "-%d", inf->end_loc);
            else if ( type->func1 == func_region_utr3_exon ) ksprintf(str, "*%d", inf->end_loc);
            else  ksprintf(str, "%d", inf->end_loc);

            // end position could cover intron region
            if ( inf->end_offset > 0 ) ksprintf(str, "+%d", inf->end_offset);
            else if ( inf->end_offset < 0 ) ksprintf(str, "%d", inf->end_offset);
        }
        kputc('=', str);
    }
}
static void generate_hgvsnom_string_exonLoss(struct mc *mc, struct mc_type *type, struct mc_inf *inf, kstring_t *str)
{
    if ( type->func1 == func_region_noncoding_exon || type->func1 == func_region_noncoding_intron) {
        ksprintf(str, "%s:n.%d", inf->transcript, inf->loc);
        if ( inf->offset > 0 ) ksprintf(str, "+%d", inf->offset);
        else if (inf->offset < 0 ) ksprintf(str, "%d", inf->offset);
    }    
    else if ( type->func1 == func_region_cds || type->func1 == func_region_intron) {
        ksprintf(str, "%s:c.%d", inf->transcript, inf->loc);
        if ( inf->offset > 0 ) ksprintf(str, "+%d", inf->offset);
        else if (inf->offset < 0 ) ksprintf(str, "%d", inf->offset);
    }
    else if ( type->func1 == func_region_utr5_exon || type->func1 == func_region_utr5_intron ) {
        ksprintf(str, "%s:c.-%d", inf->transcript, inf->loc);
        if ( inf->offset > 0 ) ksprintf(str, "+%d", inf->offset);
        else if (inf->offset < 0 ) ksprintf(str, "%d", inf->offset);
    }
    else if ( type->func1 == func_region_utr3_exon || type->func1 == func_region_utr3_intron ) {
        ksprintf(str, "%s:c.*%d", inf->transcript, inf->loc);
        if ( inf->offset > 0 ) ksprintf(str, "+%d", inf->offset);
        else if (inf->offset < 0 ) ksprintf(str, "%d", inf->offset);
    }
    else {
        generate_hgvsnom_string_empty(str);
        return; // do not generate any nomenclature
    }

    kputc('_',str);
    
    if ( type->func2 == func_region_noncoding_exon || type->func2 == func_region_noncoding_intron ||
         type->func2 == func_region_cds || type->func2 == func_region_intron ){ 
        ksprintf(str, "%d", inf->end_loc);
        if ( inf->end_offset > 0 ) ksprintf(str, "+%d", inf->end_offset);
        else if (inf->end_offset < 0 ) ksprintf(str, "%d", inf->end_offset);        
    }
    else if ( type->func2 == func_region_utr5_exon || type->func2 == func_region_utr5_intron) {
        ksprintf(str, "-%d", inf->end_loc);
        if ( inf->end_offset > 0 ) ksprintf(str, "+%d", inf->end_offset);
        else if (inf->end_offset < 0 ) ksprintf(str, "%d", inf->end_offset);        
    }
    else if ( type->func2 == func_region_utr3_exon || type->func2 == func_region_utr3_intron ) {
        ksprintf(str, "*%d", inf->end_loc);
        if ( inf->end_offset > 0 ) ksprintf(str, "+%d", inf->end_offset);
        else if (inf->end_offset < 0 ) ksprintf(str, "%d", inf->end_offset);
    }
    else {
        kputc('?', str);
    }
    kputs("del", str);
}
static char *generate_hgvsnom_string(struct mc *h)
{
    if ( h->n_tran == 0 ) return NULL;
    
    kstring_t str = {0,0,0};    

    int i;
    for ( i = 0; i < h->n_tran; ++i ) {

        if ( i ) kputc('|', &str);
        
        struct mc_type *type = &h->trans[i].type;
        struct mc_inf *inf = &h->trans[i].inf;
        

        switch (type->con1) {
            
            case mc_unknown:
            case mc_intergenic_1KB: //intergenic_1KB_variant
            case mc_upstream_1KB: // upstream_1K_of_gene
            case mc_upstream_10KB: //upstream_10K_of_gene
            case mc_downstream_1KB: //downstream_1K_of_gene
            case mc_downstream_10KB: //downstream_10K_of_gene
            case mc_intergenic: // intergenic_
            case mc_tfbs_ablation: //
            case mc_tfbs_variant:
            case mc_intragenic:
            case mc_whole_gene:
                generate_hgvsnom_string_empty(&str);
                break;

            
            case mc_noncoding_splice_region:
            case mc_noncoding_exon:
                generate_hgvsnom_string_NoncodingExon(h, type, inf, &str);
                break;

            case mc_noncoding_intron:
            case mc_coding_intron:
            case mc_utr5_intron:
            case mc_utr3_intron:
            case mc_donor_region:
            case mc_splice_donor:
            case mc_splice_acceptor:
                generate_hgvsnom_string_intron(h, type, inf, &str);
                break;

            case mc_utr5_exon:
            case mc_utr3_exon:
                generate_hgvsnom_string_UtrExon(h, type, inf, &str);
                break;

            case mc_exon_loss:
                generate_hgvsnom_string_exonLoss(h, type, inf, &str);
                break;
                
            case mc_stop_gained:
            case mc_frameshift_truncate:
            case mc_frameshift_elongation:
            case mc_inframe_indel:
            case mc_conservation_inframe_deletion:
            case mc_conservation_inframe_insertion:
            case mc_disruption_inframe_deletion:
            case mc_disruption_inframe_insertion:
            case mc_coding_variant:
                generate_hgvsnom_string_CodingDelins(h, type, inf, &str);
                break;
                
            case mc_start_loss:
            case mc_start_retained:
            case mc_stop_loss:
            case mc_stop_retained:
            case mc_missense:
            case mc_synonymous:
                // generate_hgvsnom_string_CodingSNV(h, type, inf, &str);
                generate_hgvsnom_string_CodingDelins(h, type, inf, &str);
                break;

            case mc_nocall:
                generate_hgvsnom_string_NoCall(h, type, inf, &str);
                break;
                
            default:
                error("Failed to recognise variant type. %s.", MCT[type->con1].lname);
        }

    }
    
    if ( empty_tag_string(&str) ) {
        free(str.s);
        return NULL;
    }
    
    return str.s;
}
/*
static int generate_annovar_string(struct mc *h, struct mc_type *type, struct mc_inf *inf, kstring_t *str)
{
    char *ref = inf->ref != NULL ? inf->ref : h->ref;
    char *alt = inf->alt != NULL ? inf->alt : h->alt;

    if ( inf->offset != 0 ) return 0;
    if ( str->l ) kputc(',', str);

    char *name = strdup(inf->transcript);
    int l, i;
    l = strlen(name);    
    for ( i = 0; i < l; ++i ) { if (name[i] == '.') name[i] = '\0'; }
    ksprintf(str, "%s:%s:", inf->gene, name);
    ksprintf(str, "exon%d:c.", type->count);
    free(name);
    if ( inf->end_loc != inf->loc ) {
        ksprintf(str, "%d_%d", inf->loc, inf->end_loc);        
        if ( inf->ref == NULL ) { // insertion
            ksprintf(str, "ins%s", alt);
        }        
        else {
            if ( inf->alt == NULL ) { // deletion
                ksprintf(str, "del%s",ref);
            }
            else { // del-ins
                ksprintf(str, "delins%s", alt);
            }
        }
    }
    else {
        if ( ref == NULL ) error("For developer: bad range for insertion. %s:%d", h->chr, h->start);
        if ( alt == NULL ) ksprintf(str, "%ddel%s", inf->loc, ref); //`kputs("del", str);
        else {
            if ( strlen(alt) == 1 && strlen(ref) == 1 ) ksprintf(str, "%s%d%s", ref, inf->loc, alt);
            else ksprintf(str, "%ddelins%s", inf->loc, alt);
        }
    }

    kputc(':', str);
    if ( inf->inframe_stop == 1 ) return 0;
    // amino acid changes
    else if ( type->con1 == mc_frameshift_truncate || type->con1 == mc_frameshift_elongation ) {
        if ( type->fs > 0) 
            ksprintf(str, "p.%s%dfs*", codon_short_names[type->ori_amino], type->loc_amino);
        else
            ksprintf(str, "p.%s%d%s", codon_short_names[type->ori_amino], type->loc_amino, codon_short_names[type->mut_amino]);
    }
    else { // inframe delins
        if ( type->loc_amino == type->loc_end_amino ) {
            if ( type->n ) {
                int i;
                ksprintf(str, "p.%s%ddelins", codon_short_names[type->ori_amino], type->loc_amino);
                for ( i = 0; i < type->n; ++i ) kputs(codon_short_names[type->aminos[i]], str);
            }
            else {
                ksprintf(str, "p.%s%ddel", codon_short_names[type->ori_amino], type->loc_amino);
            }
        }
        else {
            if ( type->n ) {
                if (type->con1 == mc_conservation_inframe_insertion ) {
                  inframe_insertion:
                    ksprintf(str, "p.%d_%dins", type->loc_amino, type->loc_end_amino);
                    int i;
                    for ( i = 0; i < type->n; ++i ) kputs(codon_short_names[type->aminos[i]], str);
                }
                else if ( type->con1 == mc_disruption_inframe_insertion ) {
                    if ( type->loc_amino == type->loc_end_amino ) {
                        int i;
                        ksprintf(str, "(p.%s%ddelins", codon_short_names[type->loc_amino], type->loc_amino);
                        for ( i = 0; i < type->n; ++i ) kputs(codon_short_names[type->aminos[i]], str);
                    }
                    else goto inframe_insertion;                    
                }
                else {
                    // disruption_inframe_deletion
                    int i;
                    ksprintf(str, "p.%d_%ddelins", type->loc_amino, type->loc_end_amino);
                    for ( i = 0; i < type->n; ++i ) kputs(codon_short_names[type->aminos[i]], str);
                }
            }
            else { 
                if (type->con1 == mc_conservation_inframe_deletion ) {
                    ksprintf(str, "p.%d_%ddel", type->loc_amino, type->loc_end_amino);
                }
                else if ( type->con1 == mc_disruption_inframe_deletion) {
                    ksprintf(str, "p.%d_%ddel",  type->loc_amino, type->loc_end_amino);
                }
                else {
                    ksprintf(str, "p.%s%d%s",codon_short_names[type->ori_amino], type->loc_amino, codon_short_names[type->mut_amino]);
                }
            }
        }
    }
    
    return 0;
}
// generate variant string in annovar format
// OR4F5:NM_001005484:exon1:c.T809A:p.V270E

static char *generate_annovar_name(struct mc *h)
{
    if ( h->n_tran == 0 ) return NULL;
    
    kstring_t str = {0,0,0};    

    int i;
    for ( i = 0; i < h->n_tran; ++i ) {        
        struct mc_type *type = &h->trans[i].type;
        struct mc_inf *inf = &h->trans[i].inf;        

        switch (type->con1) {
            
            case mc_unknown:
            case mc_intergenic_1KB: //intergenic_1KB_variant
            case mc_upstream_1KB: // upstream_1K_of_gene
            case mc_upstream_10KB: //upstream_10K_of_gene
            case mc_downstream_1KB: //downstream_1K_of_gene
            case mc_downstream_10KB: //downstream_10K_of_gene
            case mc_intergenic: // intergenic_
            case mc_tfbs_ablation: //
            case mc_tfbs_variant:
            case mc_intragenic:
            case mc_exon_loss:
            case mc_whole_gene:
            case mc_noncoding_splice_region:
            case mc_noncoding_exon:
            case mc_noncoding_intron:
            case mc_coding_intron:
            case mc_utr5_intron:
            case mc_utr3_intron:
            case mc_donor_region:
            case mc_splice_donor:
            case mc_splice_acceptor:
            case mc_utr5_exon:
            case mc_utr3_exon:
                break;

            case mc_stop_gained:
            case mc_frameshift_truncate:
            case mc_frameshift_elongation:
            case mc_inframe_indel:
            case mc_conservation_inframe_deletion:
            case mc_conservation_inframe_insertion:
            case mc_disruption_inframe_deletion:
            case mc_disruption_inframe_insertion:
            case mc_coding_variant:
            case mc_start_loss:
            case mc_start_retained:
            case mc_stop_loss:
            case mc_stop_retained:
            case mc_missense:
            case mc_synonymous:
            case mc_nocall:
                generate_annovar_string(h, type, inf, &str);
                break;
                
            default:
                error("Failed to recognise variant type. %s.", MCT[type->con1].lname);
        }
    }
    
    if ( empty_tag_string(&str) ) {
        free(str.s);
        return NULL;
    }
    
    return str.s;
}
*/
static char *generate_gene_string(struct mc *h)
{
    kstring_t str = {0,0,0};
    if ( h->n_tran == 0 ) {
        if ( h->inter.gene ) {
            kputs(h->inter.gene, &str);
            return str.s;
        }
        else return NULL;
    }
    
    int i;
    for ( i = 0; i < h->n_tran; ++i ) {
        if ( i ) kputc('|', &str);
        struct mc_inf *inf = &h->trans[i].inf;
        if ( inf->gene ) kputs(inf->gene, &str);
        else kputc('.', &str);
    }
    
    if (empty_tag_string(&str)) {
        free(str.s);
        return NULL;
    }

    return str.s;
}
static char *generate_transcript_string(struct mc *h)
{
    kstring_t str = {0,0,0};
    int i, k, l;
    for ( i = 0; i < h->n_tran; ++i ) {
        if ( i ) kputc('|', &str);
        struct mc_inf *inf = &h->trans[i].inf;
        kputs(inf->transcript, &str);
        l = strlen(inf->transcript);
        for ( k = 0; k < l; ++k )
            if ( inf->transcript[k] == '.') break;
    }
    if (empty_tag_string(&str)) {
        free(str.s);
        return NULL;
    }

    return str.s;
}
static char *generate_vartype_string(struct mc *h)
{
    kstring_t str = {0,0,0};
    if ( h->n_tran == 0 ) {
        kputs(MCT[h->inter.con1].sname, &str);
        return str.s;
    }

    int i;
    for ( i = 0; i < h->n_tran; ++i ) {
        if ( i ) kputc('|', &str);
        struct mc_type *type = &h->trans[i].type;
        kputs(MCT[type->con1].sname, &str);
        if ( type->con2 != mc_unknown ) {
            kputc('+', &str);
            kputs(MCT[type->con2].sname, &str);
        }
    }
    if (empty_tag_string(&str)) {
        free(str.s);
        return NULL;
    }

    return str.s;    
}
static char *generate_molecular_consequence_string(struct mc *h)
{
    kstring_t str = {0,0,0};
    int i;
    struct intergenic_core *inter = &h->inter;
    for ( i = 0; i < h->n_tran; ++i ) {
        if (i) kputc('|', &str);
        struct mc_type *type = &h->trans[i].type;
        if ( type->con1 == mc_unknown ) {
            if ( type->con2 != mc_unknown ) ksprintf(&str, "+%s", MCT[type->con2].lname);
            else kputs("unknown", &str);
        }
        else {
            kputs(MCT[type->con1].lname, &str);
            if ( type->con2 != mc_unknown ) ksprintf(&str, "+%s", MCT[type->con2].lname);
            if ( (type->con1 == mc_noncoding_intron || type->con1 == mc_coding_intron) && inter->con1 != mc_unknown ) ksprintf(&str, "+%s", MCT[inter->con1].lname);
        }
    }

    if ( i == 0 ) {
        assert(inter->con1 != mc_unknown);
        kputs(MCT[inter->con1].lname, &str);
    }

    if (empty_tag_string(&str)) {
        free(str.s);
        return NULL;
    }

    return str.s;
}
static char *generate_molecular_consequence_string_uniq(struct mc *h)
{
    kstring_t str = {0,0,0};
    int i;
    struct intergenic_core *inter = &h->inter;
    struct anno_stack *s = anno_stack_init();
    
    for ( i = 0; i < h->n_tran; ++i ) {
        struct mc_type *type = &h->trans[i].type;
        if ( type->con1 != mc_unknown ) anno_stack_push(s, (char*)MCT[type->con1].lname);
        if ( type->con2 != mc_unknown ) anno_stack_push(s, (char*)MCT[type->con2].lname);
    }

    if ( i == 0 ) {
        assert(inter->con1 != mc_unknown);
        anno_stack_push(s, (char*)MCT[inter->con1].lname);
    }
    for ( i = 0; i < s->l; ++i) {
        if ( i ) kputc('+', &str);
        kputs(s->a[i], &str);
    }
    anno_stack_destroy(s);
    return str.s;
}
static char *generate_exonintron_string(struct mc *h)
{
    kstring_t str = {0,0,0};
    int i;
    for ( i = 0; i < h->n_tran; ++i ) {
        if ( i ) kputc('|', &str);
        struct mc_type *type = &h->trans[i].type;
        struct mc_inf *inf = &h->trans[i].inf;
        if ( inf->offset != 0 )
            ksprintf(&str, "I%d", type->count);
        else {
            ksprintf(&str, "E%d", type->count);
            if ( type->count2 != 0 )
                ksprintf(&str, "/C%d", type->count2);
        }
    }
    return str.s;
}
/*
static int generate_upstream_downstream_gap_value(struct mc *h ) {
    return h->inter.TSS_dist;
}
*/
static char *generate_ivsnom_string(struct mc *h)
{
    if ( h->n_tran == 0) return NULL;
    
    kstring_t str = {0,0,0};
    int i;
    // int empty = 1;
    for ( i = 0; i < h->n_tran; ++i ) {        
        if ( i ) kputc('|', &str);
        struct mc_type *type = &h->trans[i].type;
        struct mc_inf *inf = &h->trans[i].inf;

        switch (type->con1) {
            case mc_noncoding_intron:
            case mc_coding_intron:
            case mc_utr5_intron:
            case mc_utr3_intron:
            case mc_donor_region:
            case mc_splice_donor:
            case mc_splice_acceptor:
                generate_hgvsnom_string_intron2(h, type, inf, &str);
                break;

            default:
                generate_hgvsnom_string_empty(&str);
                break;
        }
    }

    if (empty_tag_string(&str)) {
        free(str.s);
        return NULL;
    }
    return str.s;
}
/*
static char *generate_oldnom_string(struct hgvs *h)    
{
    kstring_t str = {0,0,0};
    int i;
    for ( i = 0; i < h->n_tran; ++i ) {
        if ( i ) kputc('|', &str);
        struct hgvs_inf *inf = &h->trans[i].inf;
        ksprintf(&str, "%s", inf->transcript);
        int l = strlen (inf->transcript);
        int k;
        for ( k = 0; k < l; ++k )
            if ( inf->transcript[k] == '.') break;
        if ( k == l && inf->version > 0 ) ksprintf(&str, ".%d", inf->version);
        kputc(':', &str);
        if ( inf->pos != 0 ) ksprintf(&str, "n.%d", inf->pos);
        else kputs("n.?", &str);

        if ( inf->offset > 0 ) ksprintf(&str, "+%d", inf->offset);
        else if ( inf->offset < 0 ) ksprintf(&str, "%d", inf->offset);

        if ( h->start != h->end) {
            kputc('_', &str);
            ksprintf(&str, "%d", inf->end_pos);
            if ( inf->end_offset > 0 ) ksprintf(&str,"+%d", inf->end_offset);
            else if ( inf->end_offset < 0 ) ksprintf(&str,"%d", inf->end_offset);
        }
        char *ref, *alt;
        if ( inf->strand == '+' ) {
            ref = h->ref ? safe_duplicate_string(h->ref) : NULL;
            alt = h->alt ? safe_duplicate_string(h->alt) : NULL;
        }
        else {
            ref = h->ref ? rev_seqs(h->ref, strlen(h->ref)) : NULL;
            alt = h->alt ? rev_seqs(h->alt, strlen(h->alt)) : NULL;
        }
        if ( h->type == var_type_snp ) ksprintf(&str, "%s>%s", ref, alt);
        else if ( h->type == var_type_del ) ksprintf(&str, "del%s", ref);
        else if ( h->type == var_type_ins ) ksprintf(&str, "ins%s", alt);        
        else if ( h->type == var_type_delins ) ksprintf(&str, "%s>%s", ref, alt);
        else {
            error("Failed to parse HGVS nom.");
        }
        if ( ref ) free(ref);
        if ( alt ) free(alt);
    }
    return str.s;
}
*/
static char *generate_aalength_string(struct mc *h)
{
    kstring_t str = {0,0,0};
    int i;
    for ( i = 0; i < h->n_tran; ++i ) {
        if ( i ) kputc('|', &str);
        kputw(h->trans[i].inf.aa_length, &str);
    }
    return str.s;
}


#define BRANCH(_f, _h, _l, _c, _func) do {      \
        int i;                                  \
        kstring_t str = {0,0,0};                \
        int empty = 1;                          \
        for ( i = 0; i < _f->n_allele; ++i ) {  \
            struct mc *f = _f->files[i];        \
            if ( i > 0 ) kputc(',', &str);      \
            if ( f == NULL ) {                  \
                kputc('.', &str);               \
                continue;                       \
            }                                                           \
            if ( f->type == var_type_unknow || f->type == var_type_nonref ) { \
                kputc('.', &str);                                       \
                continue;                                               \
            }                                                           \
            char *name = _func(f);                                      \
            if (name) {                                                 \
                kputs(name, &str);                                      \
                free(name);                                             \
            }                                                           \
            empty = 0;                                                  \
        }                                                               \
        if ( empty ) {                                                  \
            free(str.s);                                                \
            return 1;                                                   \
        }                                                               \
        int ret;                                                        \
        if ( _c->replace == REPLACE_MISSING ) {                          \
            ret = bcf_get_info_string(_h, _l, _c->hdr_key, &_f->tmps, &_f->mtmps); \
            if ( ret > 0 && (_f->tmps[0]!= '.' || _f->tmps[1] != 0 ) ) return 1; \
        }                                                               \
        ret = bcf_update_info_string_fixed(_h, _l, _c->hdr_key, str.s); \
        free(str.s);                                                    \
        return ret;                                                     \
    } while(0)


static int anno_hgvs_nom(struct anno_mc_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col)
{
    BRANCH(file, hdr, line, col, generate_hgvsnom_string);
}
static int anno_hgvs_ivsnom(struct anno_mc_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col)
{
    BRANCH(file, hdr, line, col, generate_ivsnom_string);
}
static int anno_hgvs_gene(struct anno_mc_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col)
{
    BRANCH(file, hdr, line, col, generate_gene_string);
}
static int anno_hgvs_trans(struct anno_mc_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col)
{
    BRANCH(file, hdr, line, col, generate_transcript_string);
}
static int anno_hgvs_vartype(struct anno_mc_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col)
{
    BRANCH(file, hdr, line, col, generate_vartype_string);
}
static int anno_hgvs_mc(struct anno_mc_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col)
{
    BRANCH(file, hdr, line, col, generate_molecular_consequence_string);
}
static int anno_hgvs_mc1(struct anno_mc_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col)
{
    BRANCH(file, hdr, line, col, generate_molecular_consequence_string_uniq);
}
static int anno_hgvs_exon(struct anno_mc_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col)
{
    BRANCH(file, hdr, line, col, generate_exonintron_string);
}
static int anno_hgvs_aalength(struct anno_mc_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col)
{
    BRANCH(file, hdr, line, col, generate_aalength_string);
}

#undef BRANCH

static int anno_hgvs_TSSdist(struct anno_mc_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col)
{
    int i;
    int *d = calloc(file->n_allele, sizeof(int));
    int empty = 1;
    for ( i = 0; i < file->n_allele; ++i ) {
        struct mc *f = file->files[i];
        struct intergenic_core *inter = &f->inter;
        d[i] = 0;
        if (inter->TSS_dist != 0) {
            d[i] = inter->TSS_dist;
            empty = 0;
        }
    }
    if ( empty == 1 ) {
        free(d);
        return 1;
    }

    int ret;
    ret = bcf_update_info_int32_fixed(hdr, line, col->hdr_key, d, file->n_allele);
    free(d);
    return ret;
}

static int anno_mc_setter2(struct anno_mc_file *file, bcf_hdr_t *hdr, bcf1_t *line)
{
    int i;

    // init transcript
    for (i = 0; i < file->n_allele; ++i ) {
        struct mc *f = file->files[i];
        
        if ( f == NULL ) continue;

        if ( f->type == var_type_unknow || f->type == var_type_nonref ) continue;

        mc_anno_trans_chunk(f, file->h);
    }

    // annotate attributes
    for ( i = 0; i < file->n_col; ++i ) {
        struct anno_col *col = &file->cols[i];

        // HGVS nomenclature, please refer to http://www.hgvs.org/ and http://varnomen.hgvs.org/
        if ( strcmp(col->hdr_key, "HGVSnom") == 0 )
            anno_hgvs_nom(file, hdr, line, col);
        // IVS nomenclature, old style
        else if ( strcmp(col->hdr_key, "IVSnom") == 0)
            anno_hgvs_ivsnom(file, hdr, line, col);
        
        // transcript names, seperated by '|'
        else if ( strcmp(col->hdr_key, "Transcript") == 0 )
            anno_hgvs_trans(file, hdr, line, col);

        // variant type, not standard
        else if ( strcmp(col->hdr_key, "VarType") == 0 )
            anno_hgvs_vartype(file, hdr, line, col);

        // Gene names, respond to transcripts, seperated by '|'
        else if ( strcmp(col->hdr_key, "Gene") == 0 )
            anno_hgvs_gene(file, hdr, line, col);

        // molecular consequence, refer to http://www.sequenceontology.org/
        // report name for each transcript, seperated by '|'
        else if ( strcmp(col->hdr_key, "MolecularConsequence") == 0 )
            anno_hgvs_mc(file, hdr, line, col);

        // molecular consequence, but report one value only
        else if ( strcmp(col->hdr_key, "MC1") == 0 )
            anno_hgvs_mc1(file, hdr, line, col);

        // exon or intron number for each transcript
        else if ( strcmp(col->hdr_key, "ExonIntron") == 0 )
            anno_hgvs_exon(file, hdr, line, col);

        // amino acid length
        else if ( strcmp(col->hdr_key, "AAlength") == 0 )
            anno_hgvs_aalength(file, hdr, line, col);

        // distance between intergenic/intragenic variant and nearby TSS
        else if ( strcmp(col->hdr_key, "TSSdistance") == 0 )
            anno_hgvs_TSSdist(file, hdr, line, col);
        
    }

    return 0;
}
//static int anno_hgvs_setter_info(struct anno_hgvs_file *file, bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col)
/*
static int anno_mc_setter(struct anno_mc_file *file, bcf_hdr_t *hdr, bcf1_t *line)
{
    int i, j;
    struct mc_handler *h = file->h;
    int empty = 1;
    kstring_t *str = malloc(file->n_col*sizeof(kstring_t));
    for ( i = 0; i < file->n_col; ++i ) memset(&str[i], 0, sizeof(kstring_t));

    for ( i = 0; i < file->n_allele; ++i ) {
        struct mc *f = file->files[i];
        if ( i > 0 ) 
            for ( j = 0; j < file->n_col; ++j ) kputc(',', &str[j]);
        if ( f == NULL ) { // for N
            for ( j = 0; j < file->n_col; ++j ) kputc('.', &str[j]);
            continue;
        }
        if ( f->type == var_type_unknow || f->type == var_type_nonref ) {
            for ( j = 0; j < file->n_col; ++j ) kputc('.', &str[j]);
            continue;
        }
        //int ret;
        int ret = mc_anno_trans_chunk(f, h);

        //continue;
        
        if ( ret == 0 ) {
            for ( j = 0; j < file->n_col; ++j ) kputc('.', &str[j]);                
            continue;
        }
        
        for ( j = 0; j < file->n_col; ++j ) {
            struct anno_col *col = &file->cols[j];
            if ( strcmp(col->hdr_key, "HGVSnom") == 0 ) {
                char *name = generate_hgvsnom_string(f);
                if (name) {
                    kputs(name, &str[j]);
                    free(name);
                }
            }
            else if ( strcmp(col->hdr_key, "Gene") == 0 ) {
                char *name = generate_gene_string(f);
                if (name) {
                    kputs(name, &str[j]);
                    free(name);
                }
            }
            else if ( strcmp(col->hdr_key, "Transcript") == 0 ) {
                char *name = generate_transcript_string(f);
                if (name) {
                    kputs(name, &str[j]);
                    free(name);                    
                }
            }
            else if ( strcmp(col->hdr_key, "VarType") == 0 ) {
                //char *name = generate_vartype_string(f);
                char *name = generate_vartype_string(f);
                if (name) {
                    kputs(name, &str[j]);
                    free(name);
                }
            }
            else if ( strcmp(col->hdr_key, "MolecularConsequence") == 0 ) {
                char *name = generate_molecular_consequence_string(f);
                if (name) {
                    kputs(name, &str[j]);
                    free(name);
                }
            }
            else if ( strcmp(col->hdr_key, "ANNOVARname") == 0 ) {
                char *name = generate_annovar_name(f);
                if (name) {
                    kputs(name, &str[j]);
                    free(name);
                }
            }
            else if ( strcmp(col->hdr_key, "MC1") == 0 ) {
                char *name = generate_molecular_consequence_string_uniq(f);
                if (name) {
                    kputs(name, &str[j]);
                    free(name);
                }
            }
            
            else if ( strcmp(col->hdr_key, "ExonIntron") == 0 ) {
                char *name = generate_exonintron_string(f);
                if (name) {
                    kputs(name, &str[j]);
                    free(name);
                }
            }
            else if ( strcmp(col->hdr_key, "AAlength") == 0 ) {
                char *name = generate_aalength_string(f);
                if (name) {
                    kputs(name, &str[j]);
                    free(name);
                }
            }            
        }
        empty = 0;
    }
    if ( empty == 1 ) {
        for ( j = 0; j < file->n_col; ++j ) free(str[j].s);
        free(str);
        return 1;
    }
    // update INFO
    for ( i = 0; i < file->n_col; ++i ) {
        if ( empty_tag_string(&str[i]) ) continue;
        struct anno_col *col = &file->cols[i];
        if ( col->replace == REPLACE_MISSING ) {
            int ret = bcf_get_info_string(hdr, line, col->hdr_key, &file->tmps, &file->mtmps);
            if ( ret > 0 && (file->tmps[0]!= '.' || file->tmps[1] != 0 ) ) continue;
        }
        bcf_update_info_string_fixed(hdr, line, col->hdr_key, str[i].s);
    }
    
    for ( i = 0; i < file->n_col; ++i ) free(str[i].s);
    free(str);
    return 0;
}
*/
struct anno_mc_file *anno_mc_file_duplicate(struct anno_mc_file *f)
{
    struct anno_mc_file *d = malloc(sizeof(*d));
    memset(d, 0, sizeof(*d));
    d->h = mc_handler_duplicate(f->h);
    d->n_allele = 0;
    d->files = NULL;
    d->n_col = f->n_col;
    d->cols = malloc(d->n_col*sizeof(struct anno_col));
    int i;
    for ( i = 0; i < d->n_col; ++i) anno_col_copy(&f->cols[i], &d->cols[i]);
    return d;
}

void anno_mc_file_destroy(struct anno_mc_file *f, int l)
{
    int i;
    // clear buffer
    for ( i = 0; i < f->n_allele; ++i ) mc_destroy(f->files[i]);
    if ( f->n_allele ) free(f->files);
    mc_handler_destroy(f->h, l);
    //if ( f->tmps) free(f->tmps);
    for ( i = 0; i < f->n_col; ++i ) free(f->cols[i].hdr_key);
    free(f->cols);
    free(f);
}

// do not update this function for now, to avoid potential bugs
struct anno_mc_file *anno_mc_file_init(bcf_hdr_t *hdr, const char *column, const char *data, const char *rna, const char *reference, const char *name_list)
{
    if ( data == NULL || rna == NULL ) return NULL;
        
    struct anno_mc_file *f = malloc(sizeof(*f));
    memset(f, 0, sizeof(*f));
    
    if ( column == NULL ) {
      full_column:
        f->n_col = 6;
        f->cols = malloc(8*sizeof(struct anno_col));
        memset(f->cols, 0, 8*sizeof(struct anno_col));
        f->cols[0].hdr_key = safe_duplicate_string("HGVSnom");
        f->cols[1].hdr_key = safe_duplicate_string("Gene");
        f->cols[2].hdr_key = safe_duplicate_string("Transcript");
        f->cols[3].hdr_key = safe_duplicate_string("VarType");
        f->cols[4].hdr_key = safe_duplicate_string("ExonIntron");
        f->cols[5].hdr_key = safe_duplicate_string("AAlength");        
        //f->cols[5].hdr_key = safe_duplicate_string("IVSnom");
        //f->cols[6].hdr_key = safe_duplicate_string("Oldnom");

    }
    else {
        kstring_t str = {0,0,0};
        kputs(column, &str);
        int i, n;
        int *s = ksplit(&str, ',', &n);
        if ( n == 0 ) {
            warnings("Failed to parse columns, try to annotate all tags. %s.", column);
            goto full_column;            
        }
        f->cols = malloc(n*sizeof(struct anno_col));
        for ( i = 0; i < n; ++i ) {
            char *ss = str.s + s[i];
            struct anno_col *col = &f->cols[f->n_col];
            memset(col, 0, sizeof(struct anno_col));
            //col->icol = -1;
            col->replace = REPLACE_MISSING;
            if ( *ss == '+' ) ss++;
            else if (*ss == '-') { col->replace = REPLACE_EXISTING; ss++; }
            if ( ss[0] == '\0') continue;
            if ( strncmp(ss, "INFO/", 5) == 0 ) ss += 5;
            if ( strcmp(ss, "HGVSnom")
                 && strcmp(ss, "Gene")
                 && strcmp(ss, "Transcript")
                 && strcmp(ss, "VarType")
                 && strcmp(ss, "ExonIntron")
                 && strcmp(ss, "IVSnom")
                 && strcmp(ss, "Oldnom")
                 && strcmp(ss, "AAlength")
                 && strcmp(ss, "ANNOVARname")
                 && strcmp(ss, "MolecularConsequence")
                 && strcmp(ss, "MC1")
                 && strcmp(ss, "TSSdistance")
                ){
                warnings("Do NOT support tag %s.", ss);
                continue;
            }
            col->hdr_key = safe_duplicate_string(ss);            
            f->n_col++;
        }
        free(str.s);
        free(s);
    }    

    // update header
#define BRANCH(_key, _description) do {                                 \
        int id;                                                         \
        id = bcf_hdr_id2int(hdr, BCF_DT_ID, _key);                      \
        if (id == -1) {                                                 \
            bcf_hdr_append(hdr, _description);                          \
            bcf_hdr_sync(hdr);                                          \
            id = bcf_hdr_id2int(hdr, BCF_DT_ID, _key);                  \
            assert(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id));        \
        }                                                               \
    } while(0)
    int i;
    for ( i = 0; i < f->n_col; ++i ) {
        struct anno_col *col = &f->cols[i];
        if ( strcmp(col->hdr_key, "HGVSnom") == 0 ) 
            BRANCH("HGVSnom", "##INFO=<ID=HGVSnom,Number=A,Type=String,Description=\"HGVS nomenclature for the description of DNA sequence variants\">");
        else if ( strcmp(col->hdr_key, "ANNOVARname") == 0 )
            BRANCH("ANNOVARname", "##INFO=<ID=ANNOVARname,Number=A,Type=String,Description=\"Variant description in ANNOVAR format.\">");
        else if ( strcmp(col->hdr_key, "Gene") == 0 ) 
            BRANCH("Gene","##INFO=<ID=Gene,Number=A,Type=String,Description=\"Gene names\">");
        else if ( strcmp(col->hdr_key, "Transcript") == 0 ) 
            BRANCH("Transcript","##INFO=<ID=Transcript,Number=A,Type=String,Description=\"Transcript names\">");
        else if ( strcmp(col->hdr_key, "VarType") == 0 ) 
            BRANCH("VarType", "##INFO=<ID=VarType,Number=A,Type=String,Description=\"Variant type.\">");
        else if ( strcmp(col->hdr_key, "MolecularConsequence") == 0 ) 
            BRANCH("MolecularConsequence", "##INFO=<ID=MolecularConsequence,Number=A,Type=String,Description=\"Predicted molecular consequence of variant.\">");
        else if ( strcmp(col->hdr_key, "MC1") == 0 ) 
            BRANCH("MC1", "##INFO=<ID=MC1,Number=A,Type=String,Description=\"Predicted molecular consequence of variant.\">");
        else if ( strcmp(col->hdr_key, "ExonIntron") == 0 ) 
            BRANCH("ExonIntron","##INFO=<ID=ExonIntron,Number=A,Type=String,Description=\"Exon/CDS or intron id on transcripts.\">");
        else if ( strcmp(col->hdr_key, "IVSnom") == 0 ) 
            BRANCH("IVSnom", "##INFO=<ID=IVSnom,Number=A,Type=String,Description=\"Old style nomenclature for the description of intron variants. Not recommand to use it.\">");
        else if ( strcmp(col->hdr_key, "Oldnom") == 0 ) 
            BRANCH("Oldnom", "##INFO=<ID=Oldnom,Number=A,Type=String,Description=\"Old style nomenclature, compared with HGVSnom use gene position instead of UTR/coding position.\">");
        else if ( strcmp(col->hdr_key, "AAlength") == 0 ) 
            BRANCH("AAlength", "##INFO=<ID=AAlength,Number=A,Type=String,Description=\"Amino acid length for each transcript. 0 for noncoding transcript.\">");
        else if ( strcmp(col->hdr_key, "TSSdistance") == 0 )
            BRANCH("TSSdistance", "##INFO=<ID=TSSdistance,Number=A,Type=Integer,Description=\"The distance between noncoding variant to nearby transcription start site. + for downstream, - for upstream.\">");
    }
#undef BRANCH

    f->h = mc_handler_init(rna, data, reference, name_list);
    return f;
}

// IMPROVE: find the gap length of nearby genes
int anno_mc_chunk(struct anno_mc_file *f, bcf_hdr_t *hdr, struct anno_pool *pool)
{
    // Update buffer for this chunk, the retrieve strategy is :
    // 1) retrieve records overlapped this chunk, if the chunk is fully covered by buffer exit this function
    // 2) if partly cover the chunk which could be interpret as start or end record located in intergenic region, update buffer by retrieving the most nearby gene (limited to 10K)

    // Notice : Gene, MOTIFs or other regulatory region will be filled in buffer
    mc_handler_fill_buffer_chunk(f->h, (char*)bcf_seqname(hdr, pool->readers[pool->i_chunk]), pool->curr_start, pool->curr_end+1);

    // annotate each record in the chunk
    int i;
    for ( i = pool->i_chunk; i < pool->n_chunk; ++i ) {
        bcf1_t *line = pool->readers[i];

        // Do we need to annotate REF allele ?
        if ( bcf_get_variant_types(line) == VCF_REF ) continue;
        
        // update buffer for one variant per time, gene and regulatory element will be treat seperately
        // regulatory element will be checked only if variant located outside of exome (intron and intergenic region will be checked)
        anno_mc_update_buffer(f, hdr, line);
        
        // generate tags
        anno_mc_setter2(f, hdr, line);
    }
    
    return 0;
}


#ifdef ANNO_MC_MAIN
#include "anno_thread_pool.h"
#include "anno_pool.h"
#include "number.h"
#include <unistd.h>

int usage()
{
    fprintf(stderr, "anno_mc [options] in.vcf\n");
    fprintf(stderr, " -data <genepred_plus.gz>   GenePredext database.\n");
    fprintf(stderr, " -rna  <rna.fa>             RNA sequence in fasta format.\n");
    fprintf(stderr, " -ref  <reference.fa>       Genome reference sequence in fasta format.\n");
    fprintf(stderr, " -tag <tag,tag>             Specify tags.\n");    
    fprintf(stderr, " -t [1]                     Threads.\n");
    fprintf(stderr, " -O <u|v|b|z>               Output format.\n");
    fprintf(stderr, " -o <output.vcf>            Output file.\n");
    fprintf(stderr, " -r [1000]                  Record per thread per time.\n");
    //fprintf(stderr, " -gene                      Selected gene list for annotation.\n");
    fprintf(stderr, " -name                      Selected transcript list for annotation.\n");
    return 1;
}
static const char *hts_bcf_wmode(int file_type)
{
    if ( file_type == FT_BCF ) return "wbu";    // uncompressed BCF
    if ( file_type & FT_BCF ) return "wb";      // compressed BCF
    if ( file_type & FT_GZ ) return "wz";       // compressed VCF
    return "w";                                 // uncompressed VCF
}
struct args {
    const char *input_fname;
    const char *output_fname;
    const char *data_fname;
    const char *rna_fname;
    const char *reference_fname;
    int n_thread;
    htsFile *fp_input;
    htsFile *fp_out;
    bcf_hdr_t *hdr_out;
    int n_record;
    //void    *gene_hash;
    void    *name_hash;
    struct anno_mc_file **files;
} args = {
    .input_fname  = NULL,
    .output_fname = NULL,
    .data_fname   = NULL,
    .rna_fname    = NULL,
    .reference_fname = NULL,
    .n_thread     = 1,
    .files        = NULL,
    .fp_input     = NULL,
    .fp_out       = NULL,
    .hdr_out      = NULL,
    .n_record     = 100,
    //.gene_hash    = NULL,
    .name_hash   = NULL,
};

void *anno_mc(void *arg, int idx) {
    struct anno_pool *pool = (struct anno_pool*)arg;
    struct args *args = (struct args*)pool->arg;
    struct anno_mc_file *f = args->files[idx];
    int i;
    //for ( i = 0; i < pool->n_reader; ++i) {
    //bcf_unpack(pool->readers[i], BCF_UN_INFO);
    while ( pool->n_chunk != pool->n_reader ) {
        update_chunk_region(pool);
        anno_mc_chunk(f, args->hdr_out, pool);
    }
//}
    return pool;
}

int parse_args(int argc, char **argv)
{
    int i;
    if ( argc == 1 )
        return usage();

    const char *tags        = 0;
    const char *thread      = 0;
    const char *record      = 0;
    const char *output_type = 0;
    //const char *gene_list_fname = 0;
    const char *trans_list_fname = 0;
    
    for ( i = 1; i < argc; ) {
        
        const char *a = argv[i++];
        const char **var = 0;
        
        if ( strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0 ) return usage();
        else if ( strcmp(a, "-data") == 0 ) var = &args.data_fname;
        else if ( strcmp(a, "-rna") == 0 ) var = &args.rna_fname;
        else if ( strcmp(a, "-tag") == 0 || strcmp(a, "-column") == 0) var = &tags;
        else if ( strcmp(a, "-t") == 0 ) var = &thread;
        else if ( strcmp(a, "-o") == 0 ) var = &args.output_fname;
        else if ( strcmp(a, "-O") == 0 ) var = &output_type;
        else if ( strcmp(a, "-r") == 0 ) var = &record;
        //else if ( strcmp(a, "-gene") == 0 || strcmp(a, "-data") == 0 ) var = &gene_list_fname;
        else if ( strcmp(a, "-name") == 0 ) var = &trans_list_fname;
        else if ( strcmp(a, "-ref") == 0 ) var = &args.reference_fname;
        if ( var != 0 ) {
            if ( i == argc) error("Missing an argument after %s", a);
            *var = argv[i++];
            continue;
        }        
        if ( args.input_fname != 0 ) error("Unknown argument, %s", a);
        if ( a[0] == '-' && a[1] ) error("Unknown parameter, %s", a);
        args.input_fname = a;
    }
    
    if ( args.data_fname == 0 ) error("Please specify bed database with -data.");
    if ( args.reference_fname == 0 ) error("Please specify genome reference database with -ref.");
    if ( args.rna_fname == 0 ) error("Please specify rna reference database with -rna.");
    if ( args.input_fname == 0 && (!isatty(fileno(stdin))) ) args.input_fname = "-";
    if ( args.input_fname == 0 ) error("No input file.");
    args.fp_input = hts_open(args.input_fname, "r");
    if ( args.fp_input == NULL ) error("%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(args.fp_input);
    if ( type.format != vcf && type.format != bcf ) error("Unsupport input format, only accept VCF/BCF.");

    int out_type = FT_VCF;
    if ( output_type != 0 ) {
        switch (output_type[0]) {
            case 'b': out_type = FT_BCF_GZ; break;
            case 'u': out_type = FT_BCF; break;
            case 'z': out_type = FT_VCF_GZ; break;
            case 'v': out_type = FT_VCF; break;
            default: error("The output type \"%d\" unrecognised", out_type);
        }
    }
    args.fp_out = args.output_fname == 0 ? hts_open("-", hts_bcf_wmode(out_type)) : hts_open(args.output_fname, hts_bcf_wmode(out_type));
    
    if ( thread ) args.n_thread = str2int((char*)thread);
    if ( record ) args.n_record = str2int((char*)record);
    if ( args.n_thread < 1) args.n_thread = 1;
    if ( args.n_record < 1 ) args.n_record = 100;

    args.hdr_out = bcf_hdr_read(args.fp_input);

    if ( args.hdr_out == NULL ) error("Failed to parse header of input.");

//    if ( gene_list_fname ) {
    //args.name_hash = name_hash_init(gene_list_fname);
    //  if ( args.gene_hash == NULL ) warnings("Failed to load gene list. %s", gene_list_fname);
//}
    if ( trans_list_fname ) {
        args.name_hash = name_hash_init(trans_list_fname);
        if ( args.name_hash == NULL ) warnings("Failed to load transcript list. %s", trans_list_fname);
    }
    
    args.files = malloc(args.n_thread*sizeof(void*));
    args.files[0] = anno_mc_file_init(args.hdr_out, tags, args.data_fname, args.rna_fname, args.reference_fname, trans_list_fname);
    //args.files[0]->gene_hash = args.gene_hash;
    //args.files[0]->name_hash = args.name_hash;
    
    for ( i = 1; i < args.n_thread; ++i ) args.files[i] = anno_mc_file_duplicate(args.files[0]);
    bcf_hdr_write(args.fp_out, args.hdr_out);
    return 0;
}

int bcfanno_mc()
{
    struct thread_pool *p = thread_pool_init(args.n_thread);
    struct thread_pool_process *q = thread_pool_process_init(p, args.n_thread*2, 0);
    struct thread_pool_result  *r;
    
    for (;;) {
        struct anno_pool *arg = anno_reader(args.fp_input, args.hdr_out, args.n_record);
        if ( arg->n_reader == 0 )
            break;
        arg->arg = &args;
        int block;
        do {
            block = thread_pool_dispatch2(p, q, anno_mc, arg, 1);
            if ( ( r = thread_pool_next_result(q) ) ) {
                // generate output
                struct anno_pool *data = (struct anno_pool*)r->data;
                int i;
                for ( i = 0; i < data->n_reader; ++i ) {
                    bcf_write1(args.fp_out, args.hdr_out, data->readers[i]);
                    bcf_destroy(data->readers[i]);
                }
                free(data->readers);
                thread_pool_delete_result(r, 1);
            }
            // flush output
        } while ( block == -1);
    }

    thread_pool_process_flush(q);
    while ( (r = thread_pool_next_result(q)) ) {
        // generate output
        struct anno_pool *data = (struct anno_pool*)r->data;
        int i;
        for ( i = 0; i < data->n_reader; ++i ) {
            bcf_write1(args.fp_out, args.hdr_out, data->readers[i]);
            bcf_destroy(data->readers[i]);
        }
        free(data->readers);
        thread_pool_delete_result(r, 1);
    }
    thread_pool_process_destroy(q);
    thread_pool_destroy(p);
    return 0;
}

void release_memory()
{
    hts_close(args.fp_input);
    hts_close(args.fp_out);

    if ( args.name_hash ) name_hash_destroy(args.name_hash);
    bcf_hdr_destroy(args.hdr_out);
    int i;
    for ( i = 0; i < args.n_thread; ++i )
        anno_mc_file_destroy(args.files[i], i);
    free(args.files);
}

int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    bcfanno_mc();

    release_memory();

    return 0;
}

#endif

