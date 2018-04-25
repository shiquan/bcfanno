#include "utils.h"
#include "gea.h"
#include "number.h"
#include "htslib/khash.h"
#include "htslib/bgzf.h"
#include "htslib/tbx.h"
#include "htslib/kstring.h"
#include "htslib/khash_str2int.h"
#include "htslib/vcf.h"
#include "hts_internal.h"
#include "htslib/kseq.h"
#include "htslib/hfile.h"

KHASH_MAP_INIT_STR(vdict, struct gea_id_info)
typedef khash_t(vdict) vdict_t;

static struct gea_id_info gea_idinfo_def = {-1, 0, -1, -1, -1, -1, NULL};

extern int gea_hdr_format(const struct gea_hdr *hdr, kstring_t *str);
extern void gea_clear(struct gea_record *v);

// parse GEA header
static const char *dump_char(char *buffer, char c)
{
    switch (c) {
    case '\n': strcpy(buffer, "\\n"); break;
    case '\r': strcpy(buffer, "\\r"); break;
    case '\t': strcpy(buffer, "\\t"); break;
    case '\'':
    case '\"':
    case '\\':
        sprintf(buffer, "\\%c", c);
        break;
    default:
        if (isprint_c(c)) sprintf(buffer, "%c", c);
        else sprintf(buffer, "\\x%02X", (unsigned char) c);
        break;
    }
    return buffer;
}

static char *find_chrom_header_line(char *s)
{
    char *nl;
    if (strncmp(s, "#chrom\t", 7) == 0) return s;
    else if ((nl = strstr(s, "\n#chrom\t")) != NULL) return nl+1;
    else return NULL;
}
int gea_hdr_add_sample(struct gea_hdr *h, const char *s)
{
    if ( !s ) return 0;

    const char *ss = s;
    while ( !*ss && isspace_c(*ss) ) ss++;
    if ( !*ss ) error("Empty sample name: trailing spaces/tabs in the header line?");
    vdict_t *d = (vdict_t*)h->dict[GEA_DT_SAMPLE];
    int ret;
    char *sdup = strdup(s);
    int k = kh_put(vdict, d, sdup, &ret);
    if (ret) { // absent
        kh_val(d, k) = gea_idinfo_def;
        kh_val(d, k).id = kh_size(d) - 1;
    }
    else {
        error_print("Duplicated sample name '%s'", s);
        free(sdup);
        return -1;
    }
    int n = kh_size(d);
    h->samples = (char**) realloc(h->samples,sizeof(char*)*n);
    h->samples[n-1] = sdup;
    h->dirty = 1;
    return 0;
}

int gea_hdr_parse_sample_line(struct gea_hdr *h, const char *str)
{
    int ret = 0;
    int i = 0;
    const char *p, *q;
    // add samples
    for (p = q = str;; ++q) {
        if (*q != '\t' && *q != 0 && *q != '\n') continue;
        if (++i > GEA_SAMPLE_COL) {
            char *s = (char*)malloc(q - p + 1);
            strncpy(s, p, q - p);
            s[q - p] = 0;
            if ( gea_hdr_add_sample(h,s) < 0 ) ret = -1;
            free(s);
        }
        if (*q == 0 || *q == '\n') break;
        p = q + 1;
    }
    gea_hdr_add_sample(h,NULL);
    return ret;
}
int gea_hdr_sync(struct gea_hdr *h)
{
    int i;
    for (i = 0; i < GEA_DICT_ALL; i++) {
        vdict_t *d = (vdict_t*)h->dict[i];
        khint_t k;
        if ( h->n[i] < kh_size(d) ) {        
            // this should be true only for GEA_DT_SAMPLE
            h->n[i] = kh_size(d);
            h->id[i] = (struct gea_id_pair*) realloc(h->id[i], kh_size(d)*sizeof(struct gea_id_pair));
        }
        for (k=kh_begin(d); k<kh_end(d); k++) {
            if (!kh_exist(d,k)) continue;
            h->id[i][kh_val(d,k).id].key = kh_key(d,k);
            h->id[i][kh_val(d,k).id].val = &kh_val(d,k);
        }
    }
    h->dirty = 0;
    return 0;
}

void gea_hrec_destroy(struct gea_hrec *hrec)
{
    free(hrec->key);
    if ( hrec->value ) free(hrec->value);
    int i;
    for (i=0; i<hrec->n_key; i++) {
        free(hrec->keys[i]);
        free(hrec->vals[i]);
    }
    free(hrec->keys);
    free(hrec->vals);
    free(hrec);
}

// Copies all fields except IDX.
struct gea_hrec *gea_hrec_dups(struct gea_hrec *hrec)
{
    struct gea_hrec *out = (struct gea_hrec*) calloc(1,sizeof(*out));
    out->type = hrec->type;
    if ( hrec->key ) out->key = strdup(hrec->key);
    if ( hrec->value ) out->value = strdup(hrec->value);
    out->n_key = hrec->n_key;
    out->keys = (char**) malloc(sizeof(char*)*hrec->n_key);
    out->vals = (char**) malloc(sizeof(char*)*hrec->n_key);
    int i, j = 0;
    for (i=0; i<hrec->n_key; i++) {
        if ( hrec->keys[i] && !strcmp("IDX",hrec->keys[i]) ) continue;
        if ( hrec->keys[i] ) out->keys[j] = strdup(hrec->keys[i]);
        if ( hrec->vals[i] ) out->vals[j] = strdup(hrec->vals[i]);
        j++;
    }
    if ( i!=j ) out->n_key -= i-j;   // IDX was omitted
    return out;
}

void gea_hrec_add_key(struct gea_hrec *hrec, const char *str, int len)
{
    int n = ++hrec->n_key;
    hrec->keys = (char**)realloc(hrec->keys, sizeof(char*)*n);
    hrec->vals = (char**)realloc(hrec->vals, sizeof(char*)*n);
    assert(len);
    hrec->keys[n-1] = (char*)malloc((len+1)*sizeof(char));
    memcpy(hrec->keys[n-1], str,len);
    hrec->keys[n-1][len] = 0;
    hrec->vals[n-1] = NULL;
}

void gea_hrec_set_val(struct gea_hrec *hrec, int i, const char *str, int len, int is_quoted)
{
    if ( !str ) { hrec->vals[i] = NULL; return; }
    if ( hrec->vals[i] ) free(hrec->vals[i]);
    if ( is_quoted ) {
        hrec->vals[i] = (char*)malloc((len+3)*sizeof(char));
        hrec->vals[i][0] = '"';
        memcpy(&hrec->vals[i][1],str,len);
        hrec->vals[i][len+1] = '"';
        hrec->vals[i][len+2] = 0;
    }
    else {
        hrec->vals[i] = (char*) malloc((len+1)*sizeof(char));
        memcpy(hrec->vals[i],str,len);
        hrec->vals[i][len] = 0;
    }
}
void gea_hrec_add_idx(struct gea_hrec *hrec, int idx)
{
    int n = ++hrec->n_key;
    hrec->keys = (char**) realloc(hrec->keys, sizeof(char*)*n);
    hrec->vals = (char**) realloc(hrec->vals, sizeof(char*)*n);
    hrec->keys[n-1] = strdup("IDX");
    kstring_t str = {0,0,0};
    kputw(idx, &str);
    hrec->vals[n-1] = str.s;
}

int gea_hrec_find_key(struct gea_hrec *hrec, const char *key)
{
    int i;
    for ( i=0; i<hrec->n_key; i++ ) 
        if ( !strcasecmp(key, hrec->keys[i]) ) return i;
    return -1;
}
static inline int is_escaped(const char *min, const char *str)
{
    int n = 0;
    while ( --str>=min && *str=='\\' ) n++;
    return n%2;
}

struct gea_hrec *gea_hdr_parse_line(const struct gea_hdr *h, const char *line, int *len)
{
    const char *p = line;
    if (p[0] != '#' || p[1] != '#') { *len = 0; return NULL; }
    p += 2;

    const char *q = p;
    while ( *q && *q != '=' && *q != '\n') q++;
    int n = q-p;
    if ( *q != '=' || !n ) { *len = q-line+1; return NULL; } // wrong format

    struct gea_hrec *hrec = (struct gea_hrec*)calloc(1,sizeof(*hrec));
    hrec->key = (char*)malloc(sizeof(char)*(n+1));
    memcpy(hrec->key, p, n);
    hrec->key[n] = 0;

    p = ++q;
    // generic field, e.g ##geatoolsVersion=0.1.0
    if ( *p != '<' ) {
        while ( *q && *q != '\n') q++;
        hrec->value = (char*)malloc((q-p+1)*sizeof(char));
        memcpy(hrec->value, p, q-p);
        hrec->value[q-p] = 0;
        *len = q - line + (*q ? 1 : 0 ); // skip \n but not \0
        return hrec;
    }

    // structure line, e.g.
    // ##sample=<ID="A1_fibroblast_H3K23me3",individual="A1",cellType="fibroblast",featureType="H3K23me3",experiment="ChIPseq",Description="Fibroblast cell.">
    // ##FORMAT=<ID="expLevel",Number=1,Type=Float,Range=0:100,Description="Expression level.">
    int nopen = 1;
    while ( *q && *q != '\n' && nopen > 0 ) {
        p = ++q;
        while (*q && *q==' ' ) { p++; q++; }
        // ^[A-Za-z_][0-9A-Za-z_.]*$
        if ( p==q && *q && (isalpha_c(*q) || *q=='_')) {
            q++;
            while (*q && (isalnum_c(*q) || *q=='_' || *q=='.')) q++;
        }
        n = q-p;
        int m = 0;
        while ( *q && *q==' ' ) { q++; m++; }
        if ( *q != '=' || !n ) {
            // wrong format
            while ( *q && *q!='\n') q++;
            error("Could not parse the header line: \"%.*s\"",(int)(q-line), line);
            *len = q-line + (*q?1:0);
            gea_hrec_destroy(hrec);
            return NULL;
        }
        gea_hrec_add_key(hrec, p, q-p-m);
        p = ++q;
        while ( *q && *q==' ' ) { p++; q++; }
        int quoted = *p=='"' ? 1 : 0;
        if ( quoted ) p++, q++;
        while ( *q && *q != '\n' ) {
            if ( quoted ) {
                if ( *q=='"' && !is_escaped(p,q) ) break;
            }
            else {
                if ( *q=='<' ) nopen++;
                if ( *q=='>' ) nopen--;
                if ( !nopen ) break;
                if ( *q==',' && nopen == 1 ) break;
            }
            q++;        
        }
        const char *r = q;
        while ( r > p && r[-1] == ' ' ) r--; // emit blank ends
        gea_hrec_set_val(hrec, hrec->n_key-1, p, r-p, quoted);
        if ( quoted && *q == '"' ) q++;
        if ( *q=='>' ) { nopen--; q++; }        
    }
    // skip to end of line
    int nonspace = 0;
    p = q;
    while ( *q && *q != '\n') { nonspace |= !isspace(*q); q++; }
    if ( nonspace ) 
        warnings("Dropped trailing junk from header line '%.*s'", (int)(q-line), line);

    *len = q - line + (*q?1:0);
    return hrec;
}
static int gea_hdr_set_idx(struct gea_hdr *hdr, int dict_type, const char *tag, struct gea_id_info *idinfo)
{
    // If available, preserve existing IDX
    if ( idinfo->id==-1 )
        idinfo->id = hdr->n[dict_type]++;
    else if ( idinfo->id < hdr->n[dict_type] && hdr->id[dict_type][idinfo->id].key ) 
        error("Conflicting IDX=%d lines in the header dictionary, the new tag is %s", idinfo->id, tag);

    if ( idinfo->id >= hdr->n[dict_type] ) hdr->n[dict_type] = idinfo->id+1;
    hts_expand0(struct gea_id_info,hdr->n[dict_type],hdr->m[dict_type],hdr->id[dict_type]);

    // NB: the next kh_put call can invalidate the idinfo pointer, therefore
    // we leave it unassigned here. It must be set explicitly in gea_hdr_sync.
    hdr->id[dict_type][idinfo->id].key = tag;
    return 0;
}
// returns: 1 when hdr needs to be synced, 0 otherwise
int gea_hdr_register_hrec(struct gea_hdr *hdr, struct gea_hrec *hrec)
{
    // contig
    int i, j = 0, ret;
    khint_t k;
    char *str;

#define BRANCH(tag) do {                        \
        i = gea_hrec_find_key(hrec, "ID");\
        if ( i<0 ) return 0;\
        str = strdup(hrec->vals[i]);\
        vdict_t *d = (vdict_t*)hdr->dict[tag];\
        khint_t k = kh_get(vdict, d, str);\
        if ( k != kh_end(d) ) { free(str); return 0; }\
        k = kh_put(vdict, d, str, &ret);\
        int idx = gea_hrec_find_key(hrec, "IDX");\
        if ( idx != -1 ) {\
            char *tmp = hrec->vals[idx];\
            idx = strtol(hrec->vals[idx], &tmp, 10);\
            if ( *tmp || idx < 0 || idx >= INT_MAX -1) {\
                warnings("Error parsing the IDX tag, skipping.");\
                return 0;\
            }\
        }\
        kh_val(d, k) = gea_idinfo_def;\
        kh_val(d, k).id = idx;\
        kh_val(d, k).hrec = hrec;\
        gea_hdr_set_idx(hdr, tag, kh_key(d,k), &kh_val(d,k));    \
        if (idx == -1)\
            gea_hrec_add_idx(hrec, kh_val(d,k).id);\
    } while(0)
        
    if ( !strcmp(hrec->key, "contig") ) {
        hrec->type = GEA_HL_CTG;

        BRANCH(GEA_DT_CTG);

        // Get the contig ID ($str) and length ($j)        
        i = gea_hrec_find_key(hrec, "length");
        if ( i<0 ) j = 0;
        else if ( sscanf(hrec->vals[i], "%d", &j)!= 1 ) return 0;
        vdict_t *d = (vdict_t*)hdr->dict[GEA_DT_CTG];
        khint_t k = kh_get(vdict, d, str);
        kh_val(d, k).info = j;        
        return 1; 
    }
    else if ( !strcmp(hrec->key, "INFO") ) {
        hrec->type = GEA_HL_INFO;
    }
    else if ( !strcmp(hrec->key, "bioType") ) {
        hrec->type = GEA_HL_BIOTYPE;
        BRANCH(GEA_DT_BIOTYPE);
        return 1;
    }
    else if ( !strcmp(hrec->key, "FORMAT") ) {
        hrec->type = GEA_HL_FMT;
    }
    else if ( !strcmp(hrec->key, "sample") ) {
        hrec->type = GEA_HL_SAMPLE;

        BRANCH(GEA_DT_SAMPLE);

        int biotype = -1;
        int individual = -1;
        int experiment = -1;
        int celltype = -1;

#define LEAF(_key, _id, _type) do {                                     \
            i = gea_hrec_find_key(hrec, _key);                          \
            if ( i >= 0 ) {                                             \
                vdict_t *d = (vdict_t*)hdr->dict[_type];                \
                khint_t k = kh_get(vdict, d, hrec->vals[i]);            \
                if ( k == kh_end(d) ) {                                 \
                    warnings("%s '%s' is not defined in the header.", _key, hrec->vals[i]);\
                    kstring_t tmp = {0,0,0};                            \
                    ksprintf(&tmp, "##%s=<ID=%s>", _key, hrec->vals[i]); \
                    int l;                                              \
                    struct gea_hrec *h = gea_hdr_parse_line(hdr, tmp.s, &l); \
                    if ( gea_hdr_register_hrec(hdr, h) == 0 )           \
                        warnings("Failed to update bioType '%s' in the header.", hrec->vals[i]); \
                    free(tmp.s);                                        \
                    k = kh_get(vdict, d, hrec->vals[i]);                \
                    _id = kh_val(d, k).id;                              \
                }                                                       \
                else _id = kh_val(d, k).id;                             \
            }                                                           \
        } while(0)

        LEAF("bioType",    biotype,    GEA_DT_BIOTYPE);
        LEAF("individual", individual, GEA_DT_INDVI);
        LEAF("experiment", experiment, GEA_DT_EXPR);
        LEAF("cellType",   celltype,   GEA_DT_CELLTYPE);
        
        vdict_t *d = (vdict_t*)hdr->dict[GEA_DT_SAMPLE];
        i = gea_hrec_find_key(hrec, "ID");
        khint_t k = kh_get(vdict, d, hrec->vals[i]);
        
        kh_val(d, k).biotype  = biotype;
        kh_val(d, k).indiv    = individual;
        kh_val(d, k).exper    = experiment;
        kh_val(d, k).celltype = celltype;
        
#undef LEAF
        return 1;
    }
    else if ( !strcmp(hrec->key, "individual") ) {
        hrec->type = GEA_HL_INDIVI;
        BRANCH(GEA_DT_INDVI);
        return 1;
    }
    else if ( !strcmp(hrec->key, "experiment") ) {
        hrec->type = GEA_HL_EXPR;
        BRANCH(GEA_DT_EXPR);
        return 1;
    }
    else if ( !strcmp(hrec->key, "cellType") ) {
        hrec->type = GEA_HL_CELL;
        BRANCH(GEA_DT_CELLTYPE);
        return 1;
    }
    else if ( hrec->n_key > 0 ) {
        hrec->type = GEA_HL_STR;
        return 1;
    }
    else return 0; // comments
#undef BRANCH
    
    char *val = NULL;
    uint32_t type = -1, var = -1;
    int num = -1, idx = -1;

    for ( i = 0; i < hrec->n_key; i++) {
        // for all lines, ID should be uniq
        if ( !strcmp(hrec->keys[i], "ID")) val = strdup(hrec->vals[i]);
        else if ( !strcmp(hrec->keys[i], "IDX")) {
            char *tmp = hrec->vals[i];
            idx = strtol(hrec->vals[i], &tmp, 10);
            if ( *tmp || idx < 0 || idx >= INT_MAX -1 ) {
                warnings("Error parsing the IDX tag, skipping.");
                return 0;
            }
        }
        // for INFO/FORMAT
        else if ( !strcmp(hrec->keys[i], "Type") ) {
            if ( !strcmp(hrec->vals[i], "Integer") ) type = GEA_TYPE_INT;
            else if ( !strcmp(hrec->vals[i], "Float") ) type = GEA_TYPE_REAL;
            else if ( !strcmp(hrec->vals[i], "String") ) type = GEA_TYPE_STR;
            else if ( !strcmp(hrec->vals[i], "Character") ) type = GEA_TYPE_STR;
            else if ( !strcmp(hrec->vals[i], "Flag") ) type = GEA_TYPE_FLAG;
            else {
                warnings("The type \"%s\" is not supported, assuming \"String\"", hrec->vals[i]);
                type = GEA_TYPE_STR;
            }
        }
        else if ( !strcmp(hrec->keys[i], "Number") ) {
            if ( !strcmp(hrec->vals[i], "block") ) var = GEA_VL_BLOCK;
            else if ( !strcmp(hrec->vals[i], ".") ) var = GEA_VL_VAR;
            else {
                sscanf(hrec->vals[i], "%d", &num);
                var = GEA_VL_FIXED;
            }
            if ( var != GEA_VL_FIXED ) num = 0xfffff;
        }
    }
    uint32_t info = ((((uint32_t)num) & 0xfffff)<<12 |
                     (var & 0xf) << 8 |
                     (type & 0xf) << 4 |
                     (((uint32_t) hrec->type) & 0xf));

    if ( !val ) return 0;
    vdict_t *d = (vdict_t*)hdr->dict[BCF_DT_ID];
    // str = strdup(hrec->key);
    k = kh_get(vdict, d, val);
    if ( k != kh_end(d) ) {
        // already present
        free(val);
        if ( kh_val(d, k).hrec ) return 0;
        kh_val(d, k).info = info;
        kh_val(d, k).hrec = hrec;
        if ( idx==-1 )
            gea_hrec_add_idx(hrec, kh_val(d, k).id);
        return 1;
    }

    k = kh_put(vdict, d, val, &ret);
    kh_val(d, k) = gea_idinfo_def;
    kh_val(d, k).info = info;
    kh_val(d, k).hrec = hrec;
    kh_val(d, k).id = idx;
    gea_hdr_set_idx(hdr, BCF_DT_ID, kh_key(d,k), &kh_val(d,k));
    if ( idx==-1 )
        gea_hrec_add_idx(hrec, kh_val(d,k).id);
    return 1;
}

int gea_hdr_add_hrec(struct gea_hdr *hdr, struct gea_hrec *hrec)
{
    if ( !hrec ) return 0;

    hrec->type = GEA_HL_GEN;
    if ( !gea_hdr_register_hrec(hdr,hrec) ) {
        // If one of the hashed field, then it is already present
        if ( hrec->type != GEA_HL_GEN ) {
            gea_hrec_destroy(hrec);
            return 0;
        }

        // Is one of the generic fields and already present?
        int i;
        for (i=0; i<hdr->n_hrec; i++) {
            if ( hdr->hrec[i]->type!=GEA_HL_GEN ) continue;
            if ( !strcmp(hdr->hrec[i]->key,hrec->key) && !strcmp(hrec->key,"fileformat") ) break;
            if ( !strcmp(hdr->hrec[i]->key,hrec->key) && !strcmp(hdr->hrec[i]->value,hrec->value) ) break;
        }
        if ( i<hdr->n_hrec ) {
            gea_hrec_destroy(hrec);
            return 0;
        }
    }
    // New record, needs to be added
    int n = ++hdr->n_hrec;
    hdr->hrec = (struct gea_hrec**) realloc(hdr->hrec, n*sizeof(void*));
    hdr->hrec[n-1] = hrec;
    hdr->dirty = 1;

    return hrec->type==GEA_HL_GEN ? 0 : 1;
}

/*
 *  Note that while querying of INFO,FMT,CTG lines is fast (the keys are hashed),
 *  the STR,GEN lines are searched for linearly in a linked list of all header lines.
 *  This may become a problem for GEAs with huge headers, we might need to build a
 *  dictionary for these lines as well.
 */
struct gea_hrec *gea_hdr_get_hrec(const struct gea_hdr *hdr, int type, const char *key, const char *value, const char *str_class)
{
    int i;
    if ( type==GEA_HL_GEN ) {    
        for (i=0; i<hdr->n_hrec; i++) {
            if ( hdr->hrec[i]->type!=type ) continue;
            if ( strcmp(hdr->hrec[i]->key,key) ) continue;
            if ( !value || !strcmp(hdr->hrec[i]->value,value) ) return hdr->hrec[i];
        }
        return NULL;
    }
    else if ( type==GEA_HL_STR ) {
        for ( i=0; i<hdr->n_hrec; i++) {
            if ( hdr->hrec[i]->type!=type ) continue;
            if ( strcmp(hdr->hrec[i]->key,str_class) ) continue;
            int j = gea_hrec_find_key(hdr->hrec[i],key);
            if ( j>=0 && !strcmp(hdr->hrec[i]->vals[j],value) ) return hdr->hrec[i];
        }
        return NULL;
    }
    vdict_t *d = type==GEA_HL_CTG ? (vdict_t*)hdr->dict[GEA_DT_CTG] : (vdict_t*)hdr->dict[GEA_DT_ID];
    khint_t k = kh_get(vdict, d, value);
    if ( k == kh_end(d) ) return NULL;
    return kh_val(d, k).hrec;
}

int gea_hdr_parse(struct gea_hdr *hdr, char *htxt)
{
    int len, needs_sync = 0, done = 0;
    char *p = htxt;

    // Check sanity: "fileformat" string must come as first
    struct gea_hrec *hrec = gea_hdr_parse_line(hdr,p,&len);
    if ( !hrec || !hrec->key || strcasecmp(hrec->key,"fileformat") )
        warnings("The first line should be ##fileformat; is the GEA header broken?");
    needs_sync += gea_hdr_add_hrec(hdr, hrec);
    
    // Parse the whole header
    do {
        while (NULL != (hrec = gea_hdr_parse_line(hdr, p, &len))) {
            needs_sync += gea_hdr_add_hrec(hdr, hrec);
            p += len;
        }

        // Next should be the sample line.  If not, it was a malformed
        // header, in which case print a warning and skip (many GEA
        // operations do not really care about a few malformed lines).
        // In the future we may want to add a strict mode that errors in
        // this case.
        if ( strncasecmp("#chrom\tchromStart",p,10) != 0 ) {
            char *eol = strchr(p, '\n');
            if (*p != '\0') warnings("Could not parse header line: %.*s", eol ? (int)(eol - p) : INT_MAX, p);
            if (eol) p = eol + 1; // Try from the next line.
            else done = -1; // No more lines left, give up.
        }
        else done = 1; // Sample line found

    } while (!done);

    if (done < 0) error("Could not parse the header, sample line not found"); // No sample line is fatal.

    int ret = gea_hdr_parse_sample_line(hdr,p);
    gea_hdr_sync(hdr);
    return ret;
}

int gea_hdr_append(struct gea_hdr *hdr, const char *line)
{
    int len;
    struct gea_hrec *hrec = gea_hdr_parse_line(hdr, (char*) line, &len);
    if ( !hrec ) return -1;
    gea_hdr_add_hrec(hdr, hrec);
    return 0;
}

void gea_hdr_remove(struct gea_hdr *hdr, int type, const char *key)
{
    int i = 0;
    struct gea_hrec *hrec;
    if ( !key ) {
        while ( i<hdr->n_hrec ) {
            if ( hdr->hrec[i]->type!=type ) { i++; continue; }
            hrec = hdr->hrec[i];

            if ( type==GEA_HL_INFO || type==GEA_HL_FMT || type== GEA_HL_CTG ) {            
                int j = gea_hrec_find_key(hdr->hrec[i], "ID");
                if ( j>=0 ) {
                    vdict_t *d = (vdict_t*)hdr->dict[type];
                    khint_t k = kh_get(vdict, d, hdr->hrec[i]->vals[j]);
                    kh_val(d, k).hrec = NULL;
                }
            }

            hdr->dirty = 1;
            hdr->n_hrec--;
            if ( i < hdr->n_hrec )
                memmove(&hdr->hrec[i],&hdr->hrec[i+1],(hdr->n_hrec-i)*sizeof(struct gea_hrec_t*));
            gea_hrec_destroy(hrec);
        }
        return;
    }
    while (1) { 
        if ( type==GEA_HL_INFO || type==GEA_HL_FMT || type==GEA_HL_CTG ) {
            
            hrec = gea_hdr_get_hrec(hdr, type, "ID", key, NULL);
            if ( !hrec ) return;

            for (i=0; i<hdr->n_hrec; i++)
                if ( hdr->hrec[i]==hrec ) break;
            assert( i<hdr->n_hrec );

            vdict_t *d = type==GEA_HL_CTG ? (vdict_t*)hdr->dict[GEA_DT_CTG] : (vdict_t*)hdr->dict[GEA_DT_ID];
            khint_t k = kh_get(vdict, d, key);
            kh_val(d, k).hrec = NULL;
        }
        else {
            for (i=0; i<hdr->n_hrec; i++) {

                if ( hdr->hrec[i]->type!=type ) continue;
                if ( type==GEA_HL_GEN ) {                
                    if ( !strcmp(hdr->hrec[i]->key,key) ) break;
                }
                else {
                    // not all structured lines have ID, we could be more sophisticated as in bcf_hdr_get_hrec()
                    int j = gea_hrec_find_key(hdr->hrec[i], "ID");
                    if ( j>=0 && !strcmp(hdr->hrec[i]->vals[j],key) ) break;
                }
            }
            if ( i==hdr->n_hrec ) return;
            hrec = hdr->hrec[i];
        }

        hdr->n_hrec--;
        if ( i < hdr->n_hrec )
            memmove(&hdr->hrec[i],&hdr->hrec[i+1],(hdr->n_hrec-i)*sizeof(struct gea_hrec*));
        gea_hrec_destroy(hrec);
        hdr->dirty = 1;
    }
}

int gea_hdr_printf(struct gea_hdr *hdr, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(NULL, 0, fmt, ap) + 2;
    va_end(ap);

    char *line = (char*)malloc(n);
    va_start(ap, fmt);
    vsnprintf(line, n, fmt, ap);
    va_end(ap);

    int ret = gea_hdr_append(hdr, line);

    free(line);
    return ret;
}


/**********************
 *** GEA header I/O ***
 **********************/

const char *gea_hdr_get_version(const struct gea_hdr *hdr)
{
    struct gea_hrec *hrec = gea_hdr_get_hrec(hdr, GEA_HL_GEN, "fileformat", NULL, NULL);
    if ( !hrec ) {
        warnings("No version string found, assuming GEA v1.0");
        return "GenomeElementAnnotation v1.0";
    }
    return hrec->value;
}

void gea_hdr_set_version(struct gea_hdr *hdr, const char *version)
{
    struct gea_hrec *hrec = gea_hdr_get_hrec(hdr, GEA_HL_GEN, "fileformat", NULL, NULL);
    if ( !hrec ) {
        int len;
        kstring_t str = {0,0,0};
        ksprintf(&str,"##fileformat=%s", version);
        hrec = gea_hdr_parse_line(hdr, str.s, &len);
        free(str.s);
    }
    else {
        free(hrec->value);
        hrec->value = strdup(version);
    }
    hdr->dirty = 1;
}

struct gea_hdr *gea_hdr_init(const char *mode)
{
    int i;
    struct gea_hdr *h;
    h = (struct gea_hdr*)calloc(1, sizeof(*h));
    if (!h) return NULL;
    for (i = 0; i < GEA_DICT_ALL; ++i)
        if ((h->dict[i] = kh_init(vdict)) == NULL) goto fail;
    if ( strchr(mode,'w') ) 
        gea_hdr_append(h, "##fileformat=GenomeElementAnnotation V1.0");

    return h;

 fail:
    for (i = 0; i < GEA_DICT_ALL; ++i)
        kh_destroy(vdict, h->dict[i]);
    free(h);
    return NULL;
}

void gea_hdr_destroy(struct gea_hdr *h)
{
    int i;
    khint_t k;
    if (!h) return;
    for (i = 0; i < GEA_DICT_ALL; ++i) {
        vdict_t *d = (vdict_t*)h->dict[i];
        if (d == 0) continue;
        for (k = kh_begin(d); k != kh_end(d); ++k)
            if (kh_exist(d, k)) free((char*)kh_key(d, k));
        kh_destroy(vdict, d);
        free(h->id[i]);
    }
    for (i=0; i<h->n_hrec; i++)
        gea_hrec_destroy(h->hrec[i]);
    if (h->n_hrec) free(h->hrec);
    if (h->samples) free(h->samples);
    //free(h->keep_samples);
    //free(h->transl[0]); free(h->transl[1]);
    free(h->mem.s);
    free(h);
}

struct gea_hdr *bea_hdr_read(htsFile *hfp)
{
    if (hfp->is_bgzf == 0) return gea_hdr_read(hfp);

    BGZF *fp = hfp->fp.bgzf;
    uint8_t magic[5];
    struct gea_hdr *h;
    h = gea_hdr_init("r");
    if (!h) error("Failed to allocate bcf header");
        

    if (bgzf_read(fp, magic, 5) != 5) {
        error_print("Failed to read the header (reading BCF in text mode?)");
        gea_hdr_destroy(h);
        return NULL;
    }
    if (strncmp((char*)magic, "BEA\2\2", 5) != 0) {
        error_print("Unrecoginsed magic key.");
        gea_hdr_destroy(h);
        return NULL;
    }
    uint8_t buf[4];
    size_t hlen;
    char *htxt = NULL;
    if (bgzf_read(fp, buf, 4) != 4) goto fail;
    hlen = buf[0] | (buf[1] << 8) | (buf[2] << 16) | (buf[3] << 24);
    if (hlen >= SIZE_MAX) { errno = ENOMEM; goto fail; }
    htxt = (char*)malloc(hlen + 1);
    if (!htxt) goto fail;
    if (bgzf_read(fp, htxt, hlen) != hlen) goto fail;
    htxt[hlen] = '\0'; // Ensure htxt is terminated
    if ( gea_hdr_parse(h, htxt) < 0 ) goto fail;
    free(htxt);
    return h;
 fail:
    error_print("Failed to read BEA header");
    free(htxt);
    gea_hdr_destroy(h);
    return NULL;
}

int bea_hdr_write(htsFile *hfp, struct gea_hdr *h)
{
    if (!h) {
        errno = EINVAL;
        return -1;
    }
    if ( h->dirty ) gea_hdr_sync(h);
    if ( hfp->is_bgzf == 0 )
        return gea_hdr_write(hfp, h);

    kstring_t htxt = {0,0,0};
    gea_hdr_format(h, &htxt);
    kputc('\0', &htxt); // include the \0 byte

    BGZF *fp = hfp->fp.bgzf;
    if ( bgzf_write(fp, "BEA\2\2", 5) !=5 ) return -1;
    uint8_t hlen[4];
    u32_to_le(htxt.l, hlen);
    if ( bgzf_write(fp, hlen, 4) !=4 ) return -1;
    if ( bgzf_write(fp, htxt.s, htxt.l) != htxt.l ) return -1;

    free(htxt.s);
    return 0;
}

int gea_check_format(const char *fn)
{
    htsFile *fp = hts_open(fn, "r");
    if ( fp == NULL ) {
        error_print("Failed to check file type of %s : %s.", fn, strerror(errno));
        return -1;
    }
    // check the file line
    int ret;
    kstring_t *s = &fp->line;
    
    for ( ;; ) {
        ret = hts_getline(fp, KS_SEP_LINE, s);
        if (ret < 0 ) break;
        if ( s->l < 36 ) break;
        if ( !strncmp("##fileformat=GenomeElementAnnotation", s->s, 36) ) {
            hts_close(fp);
            return 0;
        }
        break;
    }

    hts_close(fp);
    return -1;
}
// GEA IO
struct gea_hdr *gea_hdr_read(htsFile *fp)
{
    kstring_t txt, *s = &fp->line;
    int ret;
    struct gea_hdr *h;
    h = gea_hdr_init("r");
    if (!h) error("Failed to allocate gea header.");
    txt.l = txt.m = 0; txt.s = 0;
    while ((ret = hts_getline(fp, KS_SEP_LINE, s)) >= 0) {
        if (s->l == 0) continue;
        if (s->s[0] != '#') {
            error_print("No sample line.");
            goto error;
        }
        if (s->s[1] != '#' && fp->fn_aux) { // insert contigs here
            kstring_t tmp = { 0, 0, NULL };
            hFILE *f = hopen(fp->fn_aux, "r");
            if (f == NULL) {
                error_print("Couldn't open \"%s\"", fp->fn_aux);
                goto error;
            }
            while (tmp.l = 0, kgetline(&tmp, (kgets_func *) hgets, f) >= 0) {
                char *tab = strchr(tmp.s, '\t');
                if (tab == NULL) continue;
                kputs("##contig=<ID=", &txt); kputsn(tmp.s, tab - tmp.s, &txt);
                kputs(",length=", &txt); kputl(atol(tab), &txt);
                kputsn(">\n", 2, &txt);
            }
            free(tmp.s);
            if (hclose(f) != 0) warnings("Failed to close %s", fp->fn_aux);
        }
        kputsn(s->s, s->l, &txt);
        kputc('\n', &txt);
        if (s->s[1] != '#') break;
    }
    if ( ret < -1 ) goto error;
    if ( !txt.s )
        error("Could not read the header");

    if ( gea_hdr_parse(h, txt.s) < 0 ) goto error;

    // check tabix index, are all contigs listed in the header? add the missing ones
    tbx_t *idx = tbx_index_load(fp->fn);
    if ( idx ) {
        int i, n, need_sync = 0;
        const char **names = tbx_seqnames(idx, &n);
        for (i=0; i<n; i++) {
            struct gea_hrec *hrec = gea_hdr_get_hrec(h, GEA_HL_CTG, "ID", (char*) names[i], NULL);
            if ( hrec ) continue;
            hrec = (struct gea_hrec*) calloc(1,sizeof(*hrec));
            hrec->key = strdup("contig");
            gea_hrec_add_key(hrec, "ID", strlen("ID"));
            gea_hrec_set_val(hrec, hrec->n_key-1, (char*) names[i], strlen(names[i]), 0);
            gea_hdr_add_hrec(h, hrec);
            need_sync = 1;
        }
        free(names);
        tbx_destroy(idx);
        if ( need_sync )
            gea_hdr_sync(h);
    }
    free(txt.s);
    return h;

 error:
    free(txt.s);
    if (h) gea_hdr_destroy(h);
    return NULL;
}

int gea_hdr_set(struct gea_hdr *hdr, const char *fname)
{
    int i, n;
    char **lines = hts_readlines(fname, &n);
    if ( !lines ) return 1;
    for (i=0; i<n-1; i++) {
        int k;
        struct gea_hrec *hrec = gea_hdr_parse_line(hdr,lines[i],&k);
        if ( hrec ) gea_hdr_add_hrec(hdr, hrec);
        free(lines[i]);
    }
    gea_hdr_parse_sample_line(hdr,lines[n-1]);
    free(lines[n-1]);
    free(lines);
    gea_hdr_sync(hdr);
    return 0;
}

static void _gea_hrec_format(const struct gea_hrec *hrec, kstring_t *str)
{
    if ( !hrec->value ) {
        int j, nout = 0;
        ksprintf(str, "##%s=<", hrec->key);
        for (j=0; j<hrec->n_key; j++) {
            // do not output IDX if output is VCF
            if ( !strcmp("IDX",hrec->keys[j]) ) continue;
            if ( nout ) kputc(',',str);
            ksprintf(str,"%s=%s", hrec->keys[j], hrec->vals[j]);
            nout++;
        }
        ksprintf(str,">\n");
    }
    else
        ksprintf(str,"##%s=%s\n", hrec->key, hrec->value);
}

void gea_hrec_format(const struct gea_hrec *hrec, kstring_t *str)
{
    _gea_hrec_format(hrec,str);
}

int gea_hdr_format(const struct gea_hdr *hdr, kstring_t *str)
{
    int i;
    for (i=0; i<hdr->n_hrec; i++)
        _gea_hrec_format(hdr->hrec[i], str);

    ksprintf(str, "#chrom\tchromStart\tchromEnd\tname\tbioType\tgeneName\tstrand\tcdsStart\tcdsEnd\tblockCount\tblockStarts\tblockEnds\tINFO");
    if ( gea_hdr_nsamples(hdr) ) {
        ksprintf(str, "\tFORMAT");
        for (i=0; i<gea_hdr_nsamples(hdr); i++)
            ksprintf(str, "\t%s", hdr->samples[i]);
    }
    ksprintf(str, "\n");

    return 0;
}

char *gea_hdr_fmt_text(const struct gea_hdr *hdr, int *len)
{
    kstring_t txt = {0,0,0};
    gea_hdr_format(hdr, &txt);
    if ( len ) *len = txt.l;
    return txt.s;
}

const char **gea_hdr_seqnames(const struct gea_hdr *h, int *n)
{
    vdict_t *d = (vdict_t*)h->dict[GEA_DT_CTG];
    int tid, m = kh_size(d);
    const char **names = (const char**) calloc(m,sizeof(const char*));
    khint_t k;
    for (k=kh_begin(d); k<kh_end(d); k++) {
        if ( !kh_exist(d,k) ) continue;
        tid = kh_val(d,k).id;
        assert( tid<m );
        names[tid] = kh_key(d,k);
    }
    // sanity check: there should be no gaps
    for (tid=0; tid<m; tid++)
        assert(names[tid]);
    *n = m;
    return names;
}

int gea_hdr_write(htsFile *fp, const struct gea_hdr *h)
{
    kstring_t htxt = {0,0,0};
    gea_hdr_format(h, &htxt);
    while (htxt.l && htxt.s[htxt.l-1] == '\0') --htxt.l; // kill trailing zeros
    int ret;
    if ( fp->format.compression!=no_compression )
        ret = bgzf_write(fp->fp.bgzf, htxt.s, htxt.l);
    else
        ret = hwrite(fp->fp.hfile, htxt.s, htxt.l);
    free(htxt.s);
    return ret<0 ? -1 : 0;
}

// Data access routines
int gea_hdr_id2int(const struct gea_hdr *hdr, int which, const char *id)
{
    khint_t k;
    vdict_t *d = (vdict_t*)hdr->dict[which];
    k = kh_get(vdict, d, id);
    return k == kh_end(d) ? -1 : kh_val(d, k).id;
}

// GEA record I/O
int gea_parse(kstring_t *s, const struct gea_hdr *h, struct gea_record *v)
{
    int i = 0;
    char *p, *q, *r, *t;
    kstring_t *str;
    khint_t k;
    ks_tokaux_t aux;
    int max_n_val = 0;
    int32_t *val_a = NULL;
    int ret = -1;

    if (!s || !h || !v || !(s->s))
        return ret;

    // Assumed in lots of places, but we may as well spot this early
    assert(sizeof(float) == sizeof(int32_t));

    gea_clear(v);
    
    str = &v->shared;
    memset(&aux, 0, sizeof(ks_tokaux_t));

    for ( p = kstrtok(s->s, "\t", &aux), i = 0; p; p = kstrtok(0, 0, &aux), ++i) {
        q = (char*)aux.p;
        *q = 0;

        // chromosome
        if (i == 0) { 
            vdict_t *d = (vdict_t*)h->dict[GEA_DT_CTG];
            k = kh_get(vdict, d, p);
            if (k == kh_end(d)) {
                // warnings("Contig '%s' is not defined in the header. (Quick workaround: index the file with tabix.)", p);
                kstring_t tmp = {0,0,0};
                int l;
                ksprintf(&tmp, "##contig=<ID=%s>", p);
                struct gea_hrec *hrec = gea_hdr_parse_line(h, tmp.s, &l);
                free(tmp.s);
                if ( gea_hdr_add_hrec((struct gea_hdr*)h, hrec) )
                    gea_hdr_sync((struct gea_hdr*)h);
                k = kh_get(vdict, d, p);
                v->errcode = GEA_ERR_CTG_UNDEF;
                if (k == kh_end(d)) {
                    warnings("Could not add dummy header for contig '%s'", p);
                    v->errcode |= GEA_ERR_CTG_INVALID;
                    goto err;
                }
            }
            v->rid = kh_val(d, k).id;
        }
        // chromosome start
        else if (i == 1)  v->chromStart = str2int(p); // 0 based
        // chromosome end
        else if (i == 2 ) v->chromEnd = str2int(p);
        // name
        else if ( i == 3 ) {
            if ( !strcmp(p, ".") ) v->name = NULL;
            else v->name = strdup(p);
        }
        // biotype
        else if ( i == 4 ) {
            vdict_t *d = (vdict_t*)h->dict[GEA_DT_BIOTYPE];
            k = kh_get(vdict,d,p);
            if ( k == kh_end(d)) {
                warnings("bioType '%s' is not defined in the header.", p);
                kstring_t tmp = {0,0,0};
                int l;
                ksprintf(&tmp, "##bioType=<ID=%s>", p);
                struct gea_hrec *hrec = gea_hdr_parse_line(h, tmp.s, &l);
                free(tmp.s);

                if ( gea_hdr_add_hrec((struct gea_hdr*)h, hrec) )
                    gea_hdr_sync((struct gea_hdr*)h);
                k = kh_get(vdict, d, p);
                v->errcode = GEA_ERR_TAG_UNDEF;
                if ( k == kh_end(d)) {
                    warnings("Could not add dummy header for biological type '%s'", p);
                    goto err;
                }                
            }
            v->biotype = kh_val(d, k).id;
        }
        // gene name
        else if ( i == 5 ) {
            if ( !strcmp(p, ".") ) v->geneName = NULL;
            else v->geneName = strdup(p);
        }
        // strand
        else if (i == 6 ) {
            if ( !strcmp(p, "-") ) v->strand = strand_is_minus;
            else if ( !strcmp(p, "+") ) v->strand = strand_is_plus;
            else if ( !strcmp(p, "+-") ) v->strand = strand_is_both;
            else if ( !strcmp(p, ".") )  v->strand = strand_is_unknown;
            else {
                warnings("Unknown strand symbol, %s", p);
                v->strand = strand_is_unknown;
            }
        }
        // cdsStart
        else if ( i == 7 ) {
            if ( !strcmp(p, ".") ) v->cStart = -1;
            else v->cStart = str2int(p);  // 0 based
        }
        // cdsEnd
        else if ( i == 8 ) {
            if ( !strcmp(p, ".") ) v->cEnd = -1;
            else v->cEnd = str2int(p);
        }
        // blockCount
        else if ( i == 9 ) {
            if ( !strcmp(p, ".") ) v->blockCount = 0;
            else {
                v->blockCount = str2int(p);
                int j;
                for ( j = 0; j < 2; ++j ) 
                    v->blockPair[j] = (int*)calloc(v->blockCount, sizeof(int));
            }
        }
        // blockStarts
        else if ( i == 10 ) {
            char *ss = p;
            char *se = p;
            int j;
            for ( j = 0; j < v->blockCount; ++j )  {
                while ( se && *se != ',') se++;
                se[0] = '\0';
                v->blockPair[0][j] = str2int(ss) +1; // convert 0 based start to 1 based
                ss = ++se;
            } 
        }
        // blockEnds
        else if ( i == 11 ) {
            char *ss = p;
            char *se = p;
            int j;
            for ( j = 0; j < v->blockCount; ++j )  {
                while ( se && *se != ',') se++;
                se[0] = '\0';
                v->blockPair[1][j] = str2int(ss);
                if ( v->blockPair[1][j] < v->blockPair[0][j] ) 
                    error("%s : block of %d and %d is smaller than 0. Please notice that start is 0 based, end is 1 based.",
                          v->name, v->blockPair[0][j], v->blockPair[1][j]);
                ss = ++se;
            }
        }
        // INFO
        else if ( i == 12 ) {
            char *key;
            vdict_t *d = (vdict_t*)h->dict[GEA_DT_ID];
            v->n_info = 0;
            if ( strcmp(p, ".") ) { // not empty
                if ( *(q-1) == ';' ) *(q-1) = 0;
                for ( r = key = p;; ++r ) {
                    int c;
                    char *val, *end;
                    if (*r != ';' && *r != '=' && *r != 0) continue;
                    val = end = 0;
                    c = *r; *r = 0;
                    if (c == '=') {
                        val = r + 1;
                        for (end = val; *end != ';' && *end != 0; ++end);
                        c = *end; *end = 0;
                    }
                    else end = r;
                    if ( !*key ) {
                        if (c==0) break;
                        r = end;
                        key = r + 1;
                        continue;
                    }  // faulty GEA, ";;" in the INFO
                    k = kh_get(vdict, d, key);
                    if (k == kh_end(d) || kh_val(d, k).info == 0) {
                        warnings("INFO '%s' is not defined in the header, assuming Type=String", key);
                        kstring_t tmp = {0,0,0};
                        int l;
                        ksprintf(&tmp, "##INFO=<ID=%s,Number=1,Type=String,Description=\"Dummy\">", key);
                        struct gea_hrec *hrec = gea_hdr_parse_line(h,tmp.s,&l);
                        free(tmp.s);
                        if ( gea_hdr_add_hrec((struct gea_hdr*)h, hrec) )
                            gea_hdr_sync((struct gea_hdr*)h);
                        //d = (vdict_t*)h->dict[GEA_DT_ID];
                        k = kh_get(vdict, d, key);
                        v->errcode = GEA_ERR_TAG_UNDEF;
                        if (k == kh_end(d)) {
                            error_print("Could not add dummy header for INFO '%s'", key);
                            v->errcode |= BCF_ERR_TAG_INVALID;
                            goto err;
                        }
                    }
                    uint32_t y = kh_val(d, k).info;
                    ++v->n_info;
                    
                    bcf_enc_int1(str, kh_val(d, k).id);
                    if (val == 0) 
                        bcf_enc_size(str, 0, BCF_BT_NULL);
                    // if Flag has a value, treat it as a string
                    else if ((y>>4&0xf) == BCF_HT_FLAG || (y>>4&0xf) == BCF_HT_STR)
                        bcf_enc_vchar(str, end - val, val);
                    else { // int/float value/array
                        int i, n_val;
                        char *t, *te;
                        for (t = val, n_val = 1; *t; ++t) // count the number of values
                            if (*t == ',') ++n_val;
                        // Check both int and float size in one step for simplicity
                        if (n_val > max_n_val) {
                            int32_t *z;
                            z = (int32_t *)realloc((void *)val_a, n_val * sizeof(*z));
                            if (!z) {
                                error_print("Could not allocate memory");
                                v->errcode |= BCF_ERR_LIMITS; // No appropriate code?
                                goto err;
                            }
                            max_n_val = n_val;
                            val_a = z;
                        }
                        if ((y>>4&0xf) == BCF_HT_INT) {
                            for (i = 0, t = val; i < n_val; ++i, ++t) {
                            
                                val_a[i] = strtol(t, &te, 10);
                                // conversion failed
                                if ( te==t ) {                                
                                    val_a[i] = bcf_int32_missing;
                                    while ( *te && *te!=',' ) te++;
                                }
                                t = te;
                            }
                            bcf_enc_vint(str, n_val, val_a, -1);
                        }
                        else if ((y>>4&0xf) == BCF_HT_REAL) {
                            float *val_f = (float *)val_a;
                            for (i = 0, t = val; i < n_val; ++i, ++t) {
                            
                                val_f[i] = strtod(t, &te);
                                if ( te==t ) { // conversion failed
                                    bcf_float_set_missing(val_f[i]);
                                    while ( *te && *te!=',' ) te++;
                                }
                                t = te;
                            }
                            bcf_enc_vfloat(str, n_val, val_f);
                        }
                    }
                    if (c == 0) break;
                    r = end;
                    key = r + 1;
                }
            }
        }
        // FMT
        else if ( i == 13 ) {
            // not support yet
        }
        else if ( i > 13 ) {
            // not support yet
        }
    }
  end:
    ret = 0;

  err:
    if (val_a) free(val_a);
    
    return ret;
}

int gea_read(htsFile *fp, const struct gea_hdr *h,  struct gea_record *v)
{
    int ret;
    ret = hts_getline(fp, KS_SEP_LINE, &fp->line);
    if (ret < 0) return ret;
    return gea_parse(&fp->line, h, v);
}

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


bcf_info_t *gea_get_info_id(const struct gea_hdr *hdr, struct gea_record *b, int id)
{
    int i;
    if ( !(b->unpacked & GEA_UN_INFO) ) gea_unpack(hdr, b, GEA_UN_INFO);
    for ( i = 0; i < b->n_info; ++i ) {
        if ( b->d.info[i].key == id ) return &b->d.info[i];
    }
    return NULL;
}

bcf_info_t *gea_get_info(const struct gea_hdr *hdr, struct gea_record *b, const char *key)
{
    int id = gea_hdr_id2int(hdr, GEA_DT_ID, key);
    if ( !gea_hdr_idinfo_exists(hdr, id) ) return NULL;
    return gea_get_info_id(hdr, b, id);
}

int gea_unpack(const struct gea_hdr *hdr, struct gea_record *b, int which)
{

    uint8_t *ptr = (uint8_t*)b->shared.s, *ptr_ori;
    struct gea_dec *d = &b->d;
    struct gea_coding_transcript *c = &b->c;

    if ( which & GEA_UN_TRANS ) which |= GEA_UN_CIGAR;

    // cigar field stored in the INFO field, so parse INFO before cigar
    if ( which & GEA_UN_CIGAR ) which |= GEA_UN_INFO;
    
    if ( which & GEA_UN_INFO && !(b->unpacked&GEA_UN_INFO) ) {
        if ( !b->shared.l ) return 0;
        hts_expand(bcf_info_t, b->n_info, d->m_info, d->info);
        int i;
        for ( i = 0; i < d->m_info; ++i ) d->info[i].vptr_free = 0;
        for ( i = 0; i < b->n_info; ++i )
            ptr = bcf_unpack_info_core1(ptr, &d->info[i]);
        b->unpacked |= GEA_UN_INFO;
    }

    if ( which & GEA_UN_CIGAR && !(b->unpacked & GEA_UN_CIGAR) ) {
        b->unpacked |= GEA_UN_CIGAR;
        
        // only consider alignment state of transcript
        if ( strcmp(hdr->id[GEA_DT_BIOTYPE][b->biotype].key, "mRNA") != 0 &&
             strcmp(hdr->id[GEA_DT_BIOTYPE][b->biotype].key, "ncRNA") != 0 )
            return 0; 

        // alignmentState
        bcf_info_t *info = gea_get_info(hdr, b, "alignment_state");
        if ( info == NULL ) goto unset_alignment;
            
        char *p = (char *)info->vptr;
        int j, k, l, m = 0;
        l = info->len;
        if ( *p == bcf_str_missing || (*p == '.' && l == 1) ) goto unset_alignment;

        for ( j = 0; j < l; ) {
            k = j;
            for ( ; isdigit(p[j]) && j < l; ++j );
            if ( k == j )
                error("Failed to parse alignment_state, %s %d %d %s",
                      hdr->id[GEA_DT_CTG][b->rid].key, b->chromStart, b->chromEnd, p);
            if ( b->n_cigar == m ) {
                m = m == 0 ? 2 : m+2;
                b->cigars = (int*)realloc(b->cigars, m*sizeof(int));
            }
                
#define BRANCH(x, type) do {                                            \
                int c = str2int_l(p+k, j-k);                            \
                x = (c<<CIGAR_PACKED_FIELD) | (type & CIGAR_MASK_TYPE); \
            } while(0)
              
            switch ( p[j] ) {
                case CIGAR_UNKNOWN_BASE: {
                    if ( m > 0 ) free(b->cigars);
                    b->n_cigar = 0;
                    b->cigars = NULL;              
                    m = 0;
                    break;
                }
                    
                case CIGAR_MATCH_BASE:
                    // Symbol M denote match and mismatch, no gaps between transcript and genome.
                    BRANCH(b->cigars[b->n_cigar], CIGAR_MATCH_TYPE);
                    break;
                    
                case CIGAR_INSERT_BASE:
                    // Insertion in the transcript sequence.
                    BRANCH(b->cigars[b->n_cigar], CIGAR_INSERT_TYPE);
                    break;
                      
                case CIGAR_DELETE_BASE:
                    BRANCH(b->cigars[b->n_cigar], CIGAR_DELETE_TYPE);
                    break;
                    
                default:
                    error_print("Unknown alignment_state type %c", p[j]);
                    return 1;
                }
#undef BRANCH
            // not allocated memory for cigars, get out of realign...
            if ( m == 0 ) break;
            b->n_cigar++;
            j++;
        }
        if ( 0 ) {
          unset_alignment:
            b->n_cigar = 0;
            b->cigars = NULL;
        }
    }
    
    if ( which & GEA_UN_TRANS && !(b->unpacked&GEA_UN_TRANS) ) {
        // check if gene
        if ( strcmp(hdr->id[GEA_DT_BIOTYPE][b->biotype].key, "mRNA") != 0 &&
             strcmp(hdr->id[GEA_DT_BIOTYPE][b->biotype].key, "ncRNA") != 0 )
            return 0; 

        b->unpacked |= GEA_UN_TRANS;
        
        if ( b->strand != strand_is_plus && b->strand != strand_is_minus )
            error("Strand is not specified. Cannot parse location of this gene '%s'.", b->name);
            
        if ( b->blockCount <= 0 || b->cStart <= 0 || b->cEnd <= 0) return 0; // no exome information

        // reset structure
        memset(&b->c, 0, sizeof(struct gea_coding_transcript));
        // Calculate the length of function regions, utr5_length is the length of UTR5, and utr3_length
        // is the length of UTR3, without consider of the strand of transcript.
        // int read_length = 0;
        int utr5_length = 0;
        int utr3_length = 0;
        // int cds_length = 0;
        int loc = 0;
        int exon_start, exon_end, exon_length;
        int is_coding = strcmp(hdr->id[GEA_DT_BIOTYPE][b->biotype].key, "mRNA") == 0 ? 1 : 0;
        int i;
        for ( i = 0; i < 2; i++ ) {
            c->loc[i] = malloc(b->blockCount*sizeof(int));
            if ( c->loc[i] == NULL ) error("Failed to allocate memory.");
        }
                
        // First loop. Purpose of this loop is trying to calculate the forward and backward length.
        // Meanwhile, the related location of the transcripts block edges will also be calculated.
        for ( i = 0; i < b->blockCount; ++i ) {
            exon_start = b->blockPair[0][i]; // 0 based s
            exon_end   = b->blockPair[1][i];

            // Because we have convert 0 based start to 1 based in the first place, so need to plus 1 here.
            exon_length = exon_end - exon_start +1;
                
            // Set related location of the transcript block edges. Assume all the transcripts are plus strand in
            // the first plus, then translocate the blocks if transcript is minus.
            c->loc[0][i] = ++loc;
            loc = loc + exon_length - 1;
            c->loc[1][i] = loc;

            // Add the exon length to dna reference length. The reference length should be the sum of all exons.
            c->reference_length += exon_length;
                
            // If noncoding transcript skip next steps.
            if ( is_coding == 0 ) continue;

            // Count forward length.
            if ( b->blockPair[1][i] <= b->cStart )  utr5_length += exon_length;
            // First cds, exon consist of UTR and cds.
            else if ( b->cStart > b->blockPair[0][i] )  utr5_length += b->cStart - b->blockPair[0][i]+1;
        
            if ( b->cEnd <= b->blockPair[0][i]) utr3_length += exon_length;
            else if ( b->cEnd < b->blockPair[1][i]) utr3_length += b->blockPair[1][i] - b->cEnd;
        }
        
        // Init forward and backward length. For minus strand, reverse the backward and forward length.
        if ( b->strand == strand_is_plus ) {
            c->utr5_length = utr5_length;
            c->cds_length  = c->reference_length - utr3_length;
        }
        else {
            c->utr5_length = utr3_length;
            c->cds_length  = c->reference_length - utr5_length;
            for ( i = 0; i < b->blockCount; ++i ) {
                c->loc[0][i] = c->reference_length - c->loc[0][i] + 1;
                c->loc[1][i] = c->reference_length - c->loc[1][i] + 1;
            }
        }
        
        // realign genome locations
        if ( b->n_cigar > 0 ) {
            int match = 0, del = 0, ins = 0, offset = 0;
            int j = 0;
            int lock_utr5_realign = is_coding ? 0 : 1;
            int lock_cds_realign  = is_coding ? 0 : 1;
                
            for ( ;; ) {
                // Offset need to be added to utr5 length only if first match block smaller than it.
                if ( lock_utr5_realign == 0 && c->utr5_length <= match) {
                    c->utr5_length += offset;
                    lock_utr5_realign = 1;
                }
            
                // Adjust the cds length.
                if ( lock_cds_realign == 0 && c->cds_length >= match) {
                    c->cds_length += offset;
                    lock_cds_realign = 1;
                }

                if ( i == b->blockCount*2 || j == b->n_cigar ) break;
                if ( b->cigars[j] & CIGAR_MATCH_TYPE ) match += b->cigars[j] >> CIGAR_PACKED_FIELD;
                // deletions should be consider as match when count locs
                else if ( b->cigars[j] & CIGAR_DELETE_TYPE ) {
                    del = b->cigars[j] >> CIGAR_PACKED_FIELD;
                    match += del;
                    offset -= del;
                }
                // If insertion is tail-As, skip count offset.
                else if ( b->cigars[j] & CIGAR_INSERT_TYPE ) {
                    ins = b->cigars[j] >> CIGAR_PACKED_FIELD;
                    if ( j != b->n_cigar -1) offset += ins;
                }
                for (; i < b->blockCount*2;) {
                    int *loc = b->strand == strand_is_plus ? &c->loc[i%2][i/2] : &c->loc[i%2 ? 0 : 1][b->blockCount-i/2-1];
                    if (*loc > match) break;
                    *loc += offset;
                    i++;
                }
                j++;
            }            
        }
        else {
            error("Try to unpack locations of %s. Only works for mRNA and ncRNA.", hdr->id[GEA_DT_BIOTYPE][b->biotype].key);
        }
    }
    // todo : unpack format
    return 0;
}

int gea_format(const struct gea_hdr *h, const struct gea_record *v, kstring_t *s)
{
    int i = 0;
    gea_unpack((struct gea_hdr*)h, (struct gea_record*)v, GEA_UN_ALL);
    // Chrom
    kputs(h->id[GEA_DT_CTG][v->rid].key, s);
    // chromStart
    kputc('\t', s); kputw(v->chromStart, s);
    // chromEnd
    kputc('\t', s); kputw(v->chromEnd, s);
    // name
    kputc('\t', s);
    if ( v->name == NULL ) kputc('.', s);
    else kputs(v->name, s);
    // biotype
    kputc('\t', s);
    kputs(h->id[GEA_DT_BIOTYPE][v->biotype].key, s);
    // geneName
    kputc('\t', s);
    if ( v->geneName == NULL ) kputc('.', s);
    else kputs(v->geneName, s);
    // strand
    kputc('\t', s);
    if ( v->strand == strand_is_minus ) kputc('-', s);
    else if ( v->strand == strand_is_plus ) kputc('+', s);
    else if ( v->strand == strand_is_both ) kputs("+-", s);
    else kputc('.', s);
    // cStart
    kputc('\t', s);
    if ( v->cStart == -1 ) kputc('.', s);
    else kputw(v->cStart, s);    
    // cEnd
    kputc('\t', s);
    if ( v->cEnd == -1 ) kputc('.', s);
    else kputw(v->cEnd, s);    
    // blockCount
    kputc('\t', s); kputw(v->blockCount, s);
    // blockStarts
    kputc('\t', s);
    for ( i = 0; i < v->blockCount; ++i ) {
        kputw(v->blockPair[0][i],s);
        kputc(',',s);
    }
    if ( i == 0) kputc('.', s);
    
    // blockEnds
    kputc('\t', s);
    for ( i = 0; i < v->blockCount; ++i ) {
        kputw(v->blockPair[1][i],s);
        kputc(',',s);
    }
    if ( i == 0) kputc('.', s);
    
    // INFO
    kputc('\t', s);
    if (v->n_info) {
        int first = 1;
        for (i = 0; i < v->n_info; ++i) {
            bcf_info_t *z = &v->d.info[i];
            if ( !z->vptr ) continue;
            if ( !first ) kputc(';', s);
            first = 0;
            if (z->key >= h->n[GEA_DT_ID]) {
                error_print("Invalid GEA, the INFO index is too large");
                errno = EINVAL;
                return -1;
            }
            kputs(h->id[GEA_DT_ID][z->key].key, s);
            if (z->len <= 0) continue;
            kputc('=', s);
            if (z->len == 1) {
                switch (z->type) {
                    case BCF_BT_INT8:  if ( z->v1.i==bcf_int8_missing ) kputc('.', s); else kputw(z->v1.i, s); break;
                    case BCF_BT_INT16: if ( z->v1.i==bcf_int16_missing ) kputc('.', s); else kputw(z->v1.i, s); break;
                    case BCF_BT_INT32: if ( z->v1.i==bcf_int32_missing ) kputc('.', s); else kputw(z->v1.i, s); break;
                    case BCF_BT_FLOAT: if ( bcf_float_is_missing(z->v1.f) ) kputc('.', s); else kputd(z->v1.f, s); break;
                    case BCF_BT_CHAR:  kputc(z->v1.i, s); break;
                    default: error("Unexpected type %d", z->type); 
                }
            }
            else bcf_fmt_array(s, z->len, z->type, z->vptr);
        }
        if ( first ) kputc('.', s);
    } else kputc('.', s);

    // todo: FORMAT and individual information
    kputc('\n', s);
    return 0;
}

int gea_write_line(htsFile *fp, kstring_t *line)
{
    int ret;
    if ( line->s[line->l-1]!='\n' ) kputc('\n',line);
    if ( fp->format.compression!=no_compression )
        ret = bgzf_write(fp->fp.bgzf, line->s, line->l);
    else
        ret = hwrite(fp->fp.hfile, line->s, line->l);
    return ret==line->l ? 0 : -1;
}

int gea_write(htsFile *fp, const struct gea_hdr *h, struct gea_record *v)
{
    int ret;
    fp->line.l = 0;
    if (gea_format(h, v, &fp->line) != 0)
        return -1;
    if ( fp->format.compression!=no_compression )
        ret = bgzf_write(fp->fp.bgzf, fp->line.s, fp->line.l);
    else
        ret = hwrite(fp->fp.hfile, fp->line.s, fp->line.l);
    return ret==fp->line.l ? 0 : -1;
}


// GEA site I/O
struct gea_record *gea_init()
{
    struct gea_record *r = (struct gea_record *)calloc(1, sizeof(*r));
    return r;
}
void gea_clear(struct gea_record *v)
{
    int i;
    for (i=0; i<v->d.m_info; i++) {
        if ( v->d.info[i].vptr_free ) {
            free(v->d.info[i].vptr - v->d.info[i].vptr_off);
            v->d.info[i].vptr_free = 0;
        }
    }
    for (i=0; i<v->d.m_fmt; i++) {
        if ( v->d.fmt[i].p_free ) {
            free(v->d.fmt[i].p - v->d.fmt[i].p_off);
            v->d.fmt[i].p_free = 0;
        }
    }
    if ( (v->unpacked | GEA_UN_TRANS) && ( v->blockCount > 0) ) {
        free(v->c.loc[0]); free(v->c.loc[1]);
    }
    if ( v->blockCount > 0) {
        free(v->blockPair[0]); free(v->blockPair[1]);
        v->blockCount = 0;
    }
    if ( v->n_cigar ) { free(v->cigars); v->n_cigar = 0; v->cigars = NULL; }
    if ( v->name ) { free(v->name); v->name = NULL; }
    if ( v->geneName ) { free(v->geneName); v->geneName = NULL; }

    v->rid = v->chromStart = v->chromEnd = v->unpacked = 0;
    v->n_info = v->n_fmt = v->n_sample = 0;
    v->shared.l = v->indiv.l = 0;
    v->d.shared_dirty = 0;
    v->d.indiv_dirty  = 0;
    v->errcode = 0;
}

void gea_destroy(struct gea_record *v)
{
    gea_clear(v);
    free(v);
}

#ifdef GEA_MAIN_TEST

int main(int argc, char **argv)
{
    if ( argc != 2 ) error("Usage: geaview in.gea");

    htsFile *fp = hts_open(argv[1], "r");
    if ( fp == NULL ) error("%s : %s.", argv[1], strerror(errno));
    struct gea_hdr *hdr = gea_hdr_read(fp);
    htsFile *out = hts_open("-", "w");
    struct gea_record *v = gea_init();
    gea_hdr_write(out, hdr);
    for ( ;; ) {
        if ( gea_read(fp,hdr,v) ) break;        
        gea_write(out, hdr, v);
    }
    gea_destroy(v);
    gea_hdr_destroy(hdr);
    hts_close(fp);
    hts_close(out);
    return 0;
}
#elif defined GENE_MAIN

int main(int argc, char **argv)
{
    if ( argc != 2 )
        error("gene_main_test data.gea.gz");

    htsFile *fp = hts_open(argv[1], "r");
    if ( fp == NULL ) error("%s : %s.", argv[1], strerror(errno));
    struct gea_hdr *hdr = gea_hdr_read(fp);
    struct gea_record *v = gea_init();
    for (;;) {
        if ( gea_read(fp, hdr, v) ) break;
        if ( strcmp(hdr->id[GEA_DT_BIOTYPE][v->biotype].key, "mRNA") && strcmp(hdr->id[GEA_DT_BIOTYPE][v->biotype].key, "ncRNA") ) continue;
        gea_unpack(hdr, v, GEA_UN_CIGAR|GEA_UN_TRANS);
        int i;
        for ( i = 0; i < v->blockCount; ++i ) {
            fprintf(stdout, "%s\t%d\t%d\t%c\t%s\t%s\tEX%d\t%d\t%d\t",
                    hdr->id[GEA_DT_CTG][v->rid].key,
                    v->blockPair[0][i]-1,
                    v->blockPair[1][i],
                    v->strand == strand_is_plus ? '+' : '-',
                    v->geneName,
                    v->name,
                    v->strand == strand_is_plus ? i+1 : v->blockCount-i,
                    v->c.loc[0][i],
                    v->c.loc[1][i]
                );
            int k;
            for ( k = 0; k < 2; ++k) {
                if ( v->c.utr5_length > v->c.loc[k][i]) fprintf(stdout, "-%d\t",v->c.utr5_length-v->c.loc[k][i]+1);
                else {
                    if ( v->c.cds_length >= v->c.loc[k][i]) fprintf(stdout, "%d\t",v->c.loc[k][i]-v->c.utr5_length);
                    else fprintf(stdout, "*%d\t",v->c.loc[k][i]-v->c.cds_length+1);
                }
            }
            fputc('\n', stdout);
        }
    }
    gea_hdr_destroy(hdr);
    gea_destroy(v);
    hts_close(fp);
}



#endif
