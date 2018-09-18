#include "utils.h"
#include "motif.h"
#include "motif_encode.h"
#include "anno_pool.h"
#include "number.h"

/* encode bitcodes */

void encode_destory(struct encode *x)
{
    free(x->x);
    free(x);
}

struct encode16 *str_encode16(const char *s, int l)
{
    int i;
    struct encode16 *x = malloc(sizeof(*x));
    x->x = 0;
    x->l = 0;
    for ( i =  0; i < 16 && i < l; ++i ) x->x = x->x<<4|(_enc[s[i]]&0xf);
    return x;
}
struct encode32 *str_encode32(const char *s, int l)
{
    int i;
    struct encode32 *x = malloc(sizeof(*x));
    x->x[0] = x->x[1] = 0;
    x->l = 0;
    for ( i =  0; i < 16 && i < l; ++i ) x->x[0] = x->x[0]<<4|(_enc[s[i]]&0xf);
    for ( i = 16; i < 32 && i < l; ++i ) x->x[1] = x->x[1]<<4|(_enc[s[i]]&0xf);
    return x;
}
struct encode64 *str_encode64(const char *s, int l)
{
    int i;
    struct encode64 *x = malloc(sizeof(*x));
    x->x[0] = x->x[1] = x->x[2] = x->x[3] = 0;
    x->l = 0;
    for ( i =  0; i < 16 && i < l; ++i ) x->x[0] = x->x[0]<<4|(_enc[s[i]]&0xf);
    for ( i = 16; i < 32 && i < l; ++i ) x->x[1] = x->x[1]<<4|(_enc[s[i]]&0xf);
    for ( i = 32; i < 48 && i < l; ++i ) x->x[2] = x->x[2]<<4|(_enc[s[i]]&0xf);
    for ( i = 48; i < 64 && i < l; ++i ) x->x[3] = x->x[3]<<4|(_enc[s[i]]&0xf);
    return x;
}
struct encode16 *str_encode16_rev(const char *s, int l)
{
    int i;
    struct encode16 *x = malloc(sizeof(*x));
    x->x = 0;
    x->l = 0;
    for ( i =  0; i < 16 && i < l; ++i ) x->x = x->x<<4|(_rev[s[l-i-1]]&0xf);
    return x;
}
struct encode32 *str_encode32_rev(const char *s, int l)
{
    int i;
    struct encode32 *x = malloc(sizeof(*x));
    x->x[0] = x->x[1] = 0;
    x->l = 0;
    for ( i =  0; i < 16 && i < l; ++i ) x->x[0] = x->x[0]<<4|(_rev[s[l-i-1]]&0xf);
    for ( i = 16; i < 32 && i < l; ++i ) x->x[1] = x->x[1]<<4|(_rev[s[l-i-1]]&0xf);
    return x;
}
struct encode64 *str_encode64_rev(const char *s, int l)
{
    int i;
    struct encode64 *x = malloc(sizeof(*x));
    x->x[0] = x->x[1] = x->x[2] = x->x[3] = 0;
    x->l = 0;
    for ( i =  0; i < 16 && i < l; ++i ) x->x[0] = x->x[0]<<4|(_rev[s[l-i-1]]&0xf);
    for ( i = 16; i < 32 && i < l; ++i ) x->x[1] = x->x[1]<<4|(_rev[s[l-i-1]]&0xf);
    for ( i = 32; i < 48 && i < l; ++i ) x->x[2] = x->x[2]<<4|(_rev[s[l-i-1]]&0xf);
    for ( i = 48; i < 64 && i < l; ++i ) x->x[3] = x->x[3]<<4|(_rev[s[l-i-1]]&0xf);
    return x;
}
static struct encode* encode_init()
{
    struct encode *x = malloc(sizeof(*x));
    x->x = NULL;
    x->l = 0;
    return x;
}
struct encode *str_encode(const char *s)
{
    int l;
    l = strlen(s);
    struct encode *x = encode_init();
    if ( l == 0 ) return x; // empty struct
    if ( l <= 16 ) {
        x->l = l;
        x->x = (void*)str_encode16(s,l);
    }
    else if ( l <= 32 ) {
        x->l = l;
        x->x = (void*)str_encode32(s,l);
    }
    else if ( l <= 64 ) {
        x->l = l;
        x->x = (void*)str_encode64(s,l);
    }
    else {
        warnings("Trying to encode %d bases. Now only support 64 bases at maximal.", l);
        x->l = 64;
        x->x = (void*)str_encode64(s, 64);
    }
    return x;
}
struct encode *str_encode_rev(const char *s)
{
    int l;
    l = strlen(s);
    struct encode *x = encode_init();
    if ( l == 0 ) return x; // empty struct
    if ( l <= 16 ) {
        x->l = l;
        x->x = (void*)str_encode16_rev(s,l);
    }
    else if ( l <= 32 ) {
        x->l = l;
        x->x = (void*)str_encode32_rev(s,l);
    }
    else if ( l <= 64 ) {
        x->l = l;
        x->x = (void*)str_encode64_rev(s,l);
    }
    else {
        warnings("Trying to encode %d bases. Now only support 64 bases at maximal.", l);
        x->l = 64;
        x->x = (void*)str_encode64_rev(s, 64);
    }
    return x;
}

static void _encode_align_bits32(struct encode32 *x)
{
    int offset = (32 - x->l)*4;
    if ( offset > 0 ) {
        x->x[0] = x->x[0] << offset;
        x->x[0] = x->x[0] | (x->x[1] >> (64-offset));
        x->x[1] = (x->x[1] << offset) >> offset;
    }
}
static void _encode_align_bits64(struct encode64 *x)
{
    int offset = 64 - x->l;
    if ( offset >= 16 ) { // bits48
        offset = (offset-16)*4;
        x->x[0] = x->x[1];
        x->x[1] = x->x[2];
        x->x[2] = x->x[3];
        x->x[3] = 0;

        x->x[0] = x->x[0] << offset;
        x->x[0] = x->x[0] | (x->x[1] >> (64-offset));
        x->x[1] = x->x[1] << offset;
        x->x[1] = x->x[1] | (x->x[2] >> (64-offset));
        x->x[2] = (x->x[3] << offset) >> offset;        
    }
    else if ( offset == 16 ) { // bits48
        x->x[0] = x->x[1];
        x->x[1] = x->x[2];
        x->x[2] = x->x[3];
        x->x[3] = 0;
    }
    else if ( offset > 0 ) {
        offset = offset *4;
        x->x[0] = x->x[0] << offset;
        x->x[0] = x->x[0] | (x->x[1] >> (64-offset));
        x->x[1] = x->x[1] << offset;
        x->x[1] = x->x[1] | (x->x[2] >> (64-offset));
        x->x[2] = x->x[2] << offset;
        x->x[2] = x->x[2] | (x->x[3] >> (64-offset));
        x->x[3] = (x->x[3] << offset) >> offset;
    }    
}

#define BRANCH(_y, _x, _l, _i) do {             \
        _y = _y | ((_x>>_i)&0x1)<<(_l-_i-1);    \
    } while(0)

static void *_encode16_rev(struct encode16 *x)
{
    struct encode16 *r = malloc(sizeof(*r));
    r->l = x->l;
    r->x = 0;
    int i;
    for ( i = 0; i < x->l; ++i ) BRANCH(r->x, x->x, x->l, i);
    return (void*)r;
}
static void *_encode32_rev(struct encode32 *x)
{
    struct encode32 *r = malloc(sizeof(*r));
    r->l = x->l;
    r->x[0] = r->x[1] = 0;
    int i;
    for ( i = 0; i < 16 && i      < x->l; ++i ) BRANCH(r->x[0], x->x[1], 16, i);
    for ( i = 0; i < 16 && i + 16 < x->l; ++i ) BRANCH(r->x[1], x->x[0], 16, i);    
    
    _encode_align_bits32(r);
    
    return (void*)r;
}
static void *_encode64_rev(struct encode64 *x)
{
    struct encode64 *r = malloc(sizeof(*r));
    r->l = x->l;
    r->x[0] = r->x[1] = r->x[2] = r->x[3] = 0;
    int i;
    for ( i = 0; i < 16 && i      < x->l; ++i ) BRANCH(r->x[0], x->x[3], 16, i);
    for ( i = 0; i < 16 && i + 16 < x->l; ++i ) BRANCH(r->x[1], x->x[2], 16, i);
    for ( i = 0; i < 16 && i + 32 < x->l; ++i ) BRANCH(r->x[2], x->x[1], 16, i);
    for ( i = 0; i < 16 && i + 48 < x->l; ++i ) BRANCH(r->x[3], x->x[0], 16, i);    
    
    _encode_align_bits64(r);
    
    return (void*)r;
}

#undef BRANCH

struct encode *encode_rev(struct encode *x)
{
    struct encode *r = malloc(sizeof(*r));
    r->l = x->l;    
    if (x->l <= 16) r->x = _encode16_rev((struct encode16*)x->x);
    else if (x->l <= 32) r->x = _encode32_rev((struct encode32*)x->x);
    else r->x = _encode64_rev((struct encode64*)x);
    return r;
}
static int encode_count_bits(uint64_t x, uint64_t y, int l)
{
    int i, c = 0;
    x = x & y;    
    for ( i = 0; i < l; i++ ) {
        if ( ((x>>(i*4)) &0xf) > 0) c++;
    }
    return c;
}
/*
static uint64_t _shift_bits16(uint64_t x1, int offset, int l)
{
    uint64_t x = x1<<(offset*4);
    x = x>>((16-l)*4);
    return x;
}
static uint64_t _shift_bits32(uint64_t x1, uint64_t x2, int offset, int l)
{
    assert(l <= 16);
    uint64_t x = x1 << (offset*4);
    if ( offset <= 16 - l ) {
        x = x >> ((16-l)*4);
    }
    else {
        //int a1 = 16- offset;
        int a2 = offset+l-16;
        // 16 - a1 - a2 = 16 - 16 + offset - offset -l + 16 = 16 -l
        x = x >> (16-l);
        x = x | (x2 >>(16-a2));
    }
    return x;
}
static int _enc16_enc16_query(struct encode16 *q, struct encode16 *r, int m)
{
    uint64_t y = 0;
    int i;
    int l = r->l - q->l;
    for ( i = 0; i < l; ++i ) {
        y = _shift_bits16(r->x, i, q->l);
        if ( encode_count_bits(x, y) >= q->l -m ) return i;
    }
    return -1;
}                            
static int _enc16_enc32_query(struct encode16 *q, struct encode32 *r, int m)
{
    uint64_t y = 0;
    int i;
    int l = r->l - q->l;

    for ( i = 0; i < 16 - q->l; ++i ) {
        y = _shift_bits(r->x, i, q->l);
        if (encode_count_bits(x, u) >= q->l -m ) return i;
    }

    uint64_t x2 = r->x[1] << ((32-r->l)*4);    
    for ( i = 16-q->l; i < l; ++i ) {        
        y = _shift_bits(r->x[0], x2, i, q->l);
        if ( encode_count_bits(q->x, y) >= q->l -m ) return i;
    }
    return -1;
}
static int _enc16_enc64_query(struct encode16 *q, struct encode64 *r, int m)
{
    uint64_t y = 0;
}
static int _enc32_enc32_query(struct encode32 *q, struct encode32 *r, int m)
{
    uint64_t y[2] = {0,0};
}
static int _enc32_enc64_query(struct encode32 *q, struct encode64 *r, int m)
{
    uint64_t y[2] = {0,0};
}
static int _enc64_enc64_query(struct encode64 *q, struct encode64 *r, int m)
{
    uint64_t y[4] = {0,0,0,0};
}

// return matched location on ref, -1 on unfound
int encode_query(struct encode *q, struct encode *r, int m)
{
    if ( q->l > r->l )         return -1;
    if ( q->l <= 16 ) {
        if ( r->l <= 16 )      return _enc16_enc16_query(q, r, m);        
        else if ( r->l <= 32 ) return _enc16_enc32_query(q, r, m);
        else                   return _enc16_enc64_query(q, r, m);            
    }
    else if (q->l <= 32) {
        if (r->l <= 32 )       return _enc32_enc32_query(q, r, m);
        else                   return _enc32_enc64_query(q, r, m);
            
    }
    else                       return _enc64_enc64_query(q, r, m);
}
*/
static int _enc16_seq_query(struct encode16 *q, const char *s, int l, int m)
{
    int i;
    uint64_t x = 0;
    for ( i = 0; i < l; ++i ) {
        x = (x<<4)|(_enc[s[i]]&0xf);
        if ( i >= q->l) {
            int b = encode_count_bits(x, q->x, q->l);
            if ( b >= q->l-m) return i - q->l;
        }
    }
    return -1;
}
static int _enc32_seq_query(struct encode32 *q, const char *s, int l, int m)
{
    uint64_t x[2] = {0, 0};
    int i;
    // init
    for ( i =  0; i < 16;   ++i ) x[0] = (x[0]<<4)|(_enc[s[i]]&0xf);
    for ( i = 16; i < q->l; ++i ) x[1] = (x[1]<<4)|(_enc[s[i]]&0xf);

    
    // check
    for ( ; i < l; ++i ) {
        int b = encode_count_bits(x[0], q->x[0], 16) + encode_count_bits(x[1], q->x[1], q->l-16);
        if ( b > q->l -m ) return i-q->l;        
        x[0] = (x[0]<<4)|((x[1]>>((q->l-17)*4))&0xf);
        x[1] = (x[1]<<4)|(_enc[s[i]]&0xf);
        //x[1] = x[1]<<((32-q->l)*4)>>((32-q->l)*4);
    }
    
    return -1;
}
static int _enc64_seq_query(struct encode64 *q, const char *s, int l, int m)
{
    return -1;
}

// return matched location of seq, treat seq as streamed, return -1 on unfound
int encode_query_seq(struct encode *q, const char *s, int m)
{
    int l;
    l = strlen(s);
    if (l < q->l )          return -1;
    if ( q->l <= 16 )       return _enc16_seq_query(q->x, s, l, m);
    else if ( q->l <= 32 )  return _enc32_seq_query(q->x, s, l, m);
    else                    return _enc64_seq_query(q->x, s, l, m); 
}
int encode_query_seq_n(struct encode *q, const char *s, int l, int m)
{    
    if (l < q->l )          return -1;
    if ( q->l <= 16 )       return _enc16_seq_query(q->x, s, l, m);
    else if ( q->l <= 32 )  return _enc32_seq_query(q->x, s, l, m);
    else                    return _enc64_seq_query(q->x, s, l, m); 
}
struct motif *motif_init(char *seq)
{
    struct motif *m = malloc(sizeof(*m));
    m->n = m->m = 0;
    return m;
}
void motif_destroy(struct motif *m)
{
    free(m->name);
    encode_destory(m->enc);
    encode_destory(m->rev);
    free(m->map[0]);
    free(m->map[1]);
    free(m->map[2]);
    free(m->map[3]);
    free(m);    
}
static void motif_put_value_line(struct motif *m, float a_v, float c_v, float g_v, float t_v)
{
    if ( m->m == m->n ) {
        if ( m->n == 0 ) {
            m->m = 6;
            m->map[0] = malloc(sizeof(float)*m->m);
            m->map[1] = malloc(sizeof(float)*m->m);
            m->map[2] = malloc(sizeof(float)*m->m);
            m->map[3] = malloc(sizeof(float)*m->m);
        }
        else {
            m->m += 2;
            m->map[0] = realloc(m->map[0], sizeof(float)*m->m);
            m->map[1] = realloc(m->map[1], sizeof(float)*m->m);
            m->map[2] = realloc(m->map[2], sizeof(float)*m->m);
            m->map[3] = realloc(m->map[3], sizeof(float)*m->m);
        }
    }    
    m->map[0][m->n] = a_v;
    m->map[1][m->n] = c_v;
    m->map[2][m->n] = g_v;
    m->map[3][m->n] = t_v;
    m->n++;
}

#define BRANCH(_x,_m,_i) do {                                   \
        uint8_t b;                                              \
        int j;                                                  \
        for ( j = 0; j < 4; ++j ) {                             \
            b = b<<1;                                           \
            if (_m->map[j][i] > MOTIF_INF_PWM) b = b|0x1;       \
        }                                                       \
        _x = _x<<4|(b&0xf);                                     \
    } while(0)                                                  \

static void _motif_sync_enc16(struct motif *m)
{
    struct encode16 *x = malloc(sizeof(*x));
    x->l = m->n;
    m->enc = encode_init();
    m->enc->l = x->l;

    int i;
    for ( i = 0; i < x->l; ++i ) BRANCH(x->x, m, i);
    
    m->enc->x = (void*)x;
}
static void _motif_sync_enc32(struct motif *m)
{
    struct encode32 *x = malloc(sizeof(*x));
    x->l = m->n;
    m->enc = encode_init();
    m->enc->l = x->l;
        
    int i;
    for ( i =  0; i < 16 && i < x->l; ++i ) BRANCH(x->x[0], m, i);
    for ( i = 16; i < 32 && i < x->l; ++i ) BRANCH(x->x[1], m, i);

    m->enc->x = (void*)x;    
}
static void _motif_sync_enc64(struct motif *m)
{
    struct encode64 *x = malloc(sizeof(*x));
    x->l = m->n;
    m->enc = encode_init();
    m->enc->l = x->l;

    int i;
    for ( i =  0; i < 16 && i < x->l; ++i ) BRANCH(x->x[0], m, i);
    for ( i = 16; i < 32 && i < x->l; ++i ) BRANCH(x->x[1], m, i);
    for ( i = 32; i < 48 && i < x->l; ++i ) BRANCH(x->x[2], m, i);
    for ( i = 48; i < 64 && i < x->l; ++i ) BRANCH(x->x[3], m, i);

    m->enc->x = (void*)x;        
}
#undef BRANCH

static void motif_sync(struct motif *m)
{
    if ( m->n <= 16 )     _motif_sync_enc16(m);
    else if (m->n <= 32 ) _motif_sync_enc32(m);
    else if (m->n <= 64 ) _motif_sync_enc64(m);
    else                  error("Motif %s is longer than 64 bases, you sure about this?", m->name);

    // reversed motif sequence
    m->rev = encode_rev(m->enc);
}
struct motif **motif_read(const char *fname, int *_n)
{
    BGZF *fp = bgzf_open(fname, "r");    
    if ( fp == NULL ) error("%s: %s.", fname, strerror(errno));
    kstream_t *ks = ks_init(fp);
    kstring_t str = {0,0,0};
    int ret;
    int n;
    int n_motif = 0, m_motif = 1;
    struct motif **mm = malloc(1*sizeof(void*));
    
    while ( ks_getuntil(ks, 2, &str, &ret) >= 0 ) {

        if ( str.l == 0 ) continue; // empty line

        int *s = ksplit(&str, 0, &n);
        char *name = str.s + s[0];
        if ( *name == '>') { // motif name            
            struct motif *m = motif_init(name);
            m->name = strdup(name+1);
            mm[n_motif] = m;
            n_motif++;
            if (n_motif == m_motif) {
                m_motif = m_motif*2;
                mm = realloc(mm, m_motif*sizeof(void*));
            }
        }
        else {
            if ( n != 4 ) error("Line %s like to be truncated.", str.s);
            float a_v = atof(str.s+s[0]);
            float c_v = atof(str.s+s[1]);
            float g_v = atof(str.s+s[2]);
            float t_v = atof(str.s+s[3]);
            motif_put_value_line(mm[n_motif-1], a_v, c_v, g_v, t_v);
        }
    }
    int i;    
    for ( i = 0; i < n_motif; ++i ) motif_sync(mm[i]);
    bgzf_close(fp);
    *_n = n_motif;
    if ( n_motif == 0 ) { free(mm); return NULL; }
    return mm;
}
static int test_motif(const char *fasta_fname, const char *motif_fname)
{
    int i, nm;
    struct motif **mm = motif_read(motif_fname, &nm);
    LOG_print("Trying to read MOTIF file.");
    for (i = 0; i < nm; ++i ) {
        struct motif *m = mm[i];
        int j;
        printf("%s\n", m->name);
        for ( j =0; j < m->n; ++j) {
            LOG_print("%f\t%f\t%f\t%f", m->map[0][j], m->map[1][j], m->map[2][j], m->map[3][j]);
        }
    }

    BGZF *bgzf = bgzf_open(fasta_fname, "r");
    if ( bgzf == NULL ) error("%s: %s.", fasta_fname, strerror(errno));

    int c;
    int pos = -1;
    uint64_t x = 0;

    int n = 0;
    int l = 0, m = 2;
    char *name = malloc(2);
    // uint64_t mask = (1LL<<(enc->l*4)) -1;
    /*
    while (1) {
        c = bgzf_getc(bgzf);
        if (c < 0 ) break;
        if ( c == '>') {
            l = 0; // reset chromosome name
            while ( (c = bgzf_getc(bgzf))>=0 && c != '\n') {
                if (l == m) {
                    m += 10;
                    name = realloc(name, m);
                }
                name[l++] = c;
            }
            name[l] = '\0';
            if (c != '\n') while((c=bgzf_getc(bgzf))>0 && c != '\n');
            pos = -1; // reset from -1
        }
        else {
            if ( c == '\n'||isspace(c)) continue;
            pos++; // start from 0
            if ( c== 'N') { x = 0; n = 0; }
            else {
                x = (x<<4)|(tab[c]&0xf);
                
                for (j = 0; j
                if ( ++n >= enc->l ) {
                    int b = countbits(x, enc);
                    if (b >= enc->l-args.mis) {
                        printf("%s\t%d\t%d\n", name, pos-enc->l+1, pos+1);
                    }
                }
            }
        }
    }
    */
    free(name);
    bgzf_close(bgzf);
    

    
    return 0;
}
float motif_pwm_score(struct motif *m, char *s, int r)
{
    int l;
    l = strlen(s);
    if ( l < m->n ) error("Sequence is shorter than motif length.");
    float pwm = 0;
    int i;
    for ( i = 0; i < m->n; ++i ) {
        int b = r == 0 ? B4[s[i]] : 5-B4[s[m->n-i-1]];
        if ( b == 0 ) error("Unknown base at %s.",s);
        b -= 1;
        pwm += m->map[b][i];
    }
    return pwm;
}
struct MTF *MTF_init()
{
    struct MTF *MTF = malloc(sizeof(struct MTF));
    memset(MTF, 0, sizeof(struct MTF));
    MTF->tid = -1;
    // static struct plp_ref ref = {{NULL, NULL}, {-1, -1}, {0,0}};
    MTF->r = plp_ref_init();
    return MTF;
}
void MTF_destory(struct MTF *m)
{
    fai_destroy(m->fai);
    if (m->r->ref[0]) free(m->r->ref[0]);
    if (m->r->ref[1]) free(m->r->ref[1]);
    free(m);
}
//  PWM_score_change
//  CTCF_PWM_score_change
static int MTF_new_region_init(struct MTF *MTF, int tid, int pos)
{
    int ret;
    int has_ref;
    struct bedaux *bed = MTF->bed;
    char *seqname = (char*)MTF->bcf_hdr->id[BCF_DT_CTG][tid].key;
    ret = bed_position_covered(bed, seqname, pos, &MTF->start, &MTF->end);
    if (ret == 0) return 0;
    
    has_ref = plp_get_ref(MTF->r, seqname, MTF->fai, tid, &MTF->ref, &MTF->ref_len);
    if ( has_ref == 0 ) error("No such chromosome %s at reference.", seqname);

    return 1;
}
static int MTF_vcf_sync(struct MTF *MTF, bcf1_t *l)
{
    if ( MTF->tid != l->rid || MTF->end < l->pos ) 
        return MTF_new_region_init(MTF, l->rid, l->pos+1);
    return 1;
}
int anno_motif_setter_info_float(bcf_hdr_t *hdr, bcf1_t *line, struct anno_col *col, int n, float *v)
{
    return bcf_update_info_float_fixed(hdr, line, col->hdr_key, v, n);
}
/*
  annotate PWM changes in VCF
 */
static char *_construct_alt_seq(char *s, int loc, char *ref, char *alt, int l)
{    
    if (ref != NULL && _enc[s[loc]] != _enc[*ref]) {
        warnings("Inconsistant ref bases, %c vs %s", s[loc], ref);
        return NULL;
    }
        
    kstring_t str = {0,0,0};
    // put capped bases into string
    if (loc) kputsn(s, loc, &str);
    l = l - str.l;
    // put alternative sequences into string
    if (alt) kputs(alt, &str);
    // s point to the next position of ref
    s = s+loc+strlen(ref);
    // put enough sequences into string
    kputsn(s, l, &str);
    return str.s;
}
float find_best_match(struct motif *m, char *seq, int *loc, int mis, int l, int strand)
{
    *loc = encode_query_seq_n(m->enc, seq, l, mis);
    if ( *loc == -1 ) return -INT_MAX;
    int best_loc = *loc;
    float best_value = motif_pwm_score(m, seq+*loc, strand);
    int ll = *loc; // last offset
    char *new_seq = seq + ll + 1;
    int new_l = l - ll -1;

    for (;;) {

        if ( ll + m->n > l ) break;

        int lp = encode_query_seq_n(m->enc, new_seq, new_l, mis);
        if ( lp == -1 ) break;
        float s = motif_pwm_score(m, new_seq+lp, strand);

        
        new_seq = new_seq + lp + 1;
        new_l = new_l - lp -1;

        ll = new_seq - seq -1;

        if ( s > best_value ) { best_value = s; best_loc = ll; }
    }
    *loc = best_loc;
    return best_value;
}
int anno_vcf_motif_pwm(struct MTF *MTF, bcf1_t *line)
{
    extern int is_atcg(char *s);
    
    if ( MTF_vcf_sync(MTF, line) == 0 ) return 0;

    float pwm_change = 0.0;
    int i;
    // debug_print("%s\t%d",MTF->bcf_hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1);
    for ( i = 0; i < MTF->n; ++i ) { // for each motif, calculate PWM

        // motif struct
        struct motif *m = MTF->mm[i];
        //debug_print("%s", m->name);
        // variant must located inside the motif
        int start = line->pos - m->n;        
        int l = m->n*2; // length of scan region
        
        // in case at the end of chromosome
        if ( start +l > MTF->ref_len) l = MTF->ref_len - start;
        // if scan region smaller than motif length
        if ( l < m->n ) continue;

        int loc;
        char *ref = MTF->ref + start;
        int strand = 0;
        float ref_v = find_best_match(m, ref, &loc, 1, l, 0);
        
        if ( loc == -1) {
            float v = find_best_match(m, ref, &loc, 1, l, 1);            
            if ( loc != -1 ) {
                strand = 1;
                ref_v = v;
            }
            else continue;
        }
        /*
        loc = encode_query_seq_n(m->enc, ref, l, 1);
        if ( loc == -1 ) { // unfound
            loc = encode_query_seq_n(m->rev, ref, l, 1);
            if ( loc != -1) strand = 1;
            else continue;
        }
        */
        ref = ref + loc;
        
        // double ref_v = motif_pwm_score(m, ref, strand); // need test
        if (ref_v < MTF->min) continue;
        
        //kstring_t str1 = {0,0,0};
        //kputsn(ref, m->n, &str1);
        //debug_print("%s\t%d\t%s\t%f\t%s", MTF->bcf_hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1,m->name, ref_v, str1.s);
        

        int k;
        float *d = malloc((line->n_allele)*sizeof(float));
        d[0] = ref_v;
        int var_loc = line->pos - start - loc;
        for ( k =1; k < line->n_allele; ++k ) {
            d[k] = -99999;

            // for complex variants, skip
            if ( is_atcg(line->d.allele[k]) ) continue;
            
            char *alt = _construct_alt_seq(ref, var_loc, line->d.allele[0], line->d.allele[k], m->n);
            if ( alt == NULL ) error("Failed to construct alternative alleles. %s %d", MTF->bcf_hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1);
            //str1.l = 0;
            //kputsn(alt, m->n, &str1);
            //debug_print("%s", str1.s);
            d[k] = motif_pwm_score(m, alt, strand);
            float change = d[k] - d[0];
            if ( fabsf(change) > fabsf(pwm_change) ) pwm_change = change;
            free(alt);
            // looks like a regularory variants?
            if ( (d[0] >= 0 || d[k] >=0) && fabsf(d[k]) > 20.0 )
                LOG_print("Regulatory variants candidate: %s\t%d\t%s\t%s\t%.4f",
                          MTF->bcf_hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1,
                          line->d.allele[k] == NULL ? "." : line->d.allele[k], m->name, d[k-1]);
        }
        //free(str1.s);
        // update INFO MOTIF_PWMscore
        MTF->cols[i].func.pwm(MTF->bcf_hdr, line, &MTF->cols[i], line->n_allele, d);
    }

    // update INFO PWM_score_change
    if ( pwm_change != 0 ) 
        MTF->ccol->func.pwm(MTF->bcf_hdr, line, MTF->ccol, 1, &pwm_change);
    
    return 0;
}

int usage()
{
    fprintf(stderr, "PWM_change_score\n");
    fprintf(stderr, " -vcf          Variants in BCF/VCF format.\n");
    fprintf(stderr, " -bed          ATAC peak region in bed format.\n");
    fprintf(stderr, " -motif        Motif file.\n");    
    fprintf(stderr, " -ref          Reference in FASTA format.\n");
    fprintf(stderr, " -O <b|u|z|v>  Output format.\n");
    fprintf(stderr, " -o            Output file. Default is stdout.\n");
    fprintf(stderr, " -t            Threads.[5]\n");
    fprintf(stderr, " -min          Mininal PWM score, skip motifs below this value.\n");
    fprintf(stderr, " -record       Records per thread.\n");
    fprintf(stderr, "\n");
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
    const char *bed_fname;
    const char *motif_fname;
    const char *ref_fname;
    
    int n;
    struct motif **mm;    
    struct anno_col *pwm_cols;

    int motif_min;
    
    //faidx_t *fai;
    bcf_hdr_t *bcf_hdr;
    struct bedaux *bed;
    htsFile *fp_in;
    htsFile *fp_out;
    struct anno_col *pcs_col;
    int n_thread;
    int n_record;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .bed_fname = NULL,
    .motif_fname = NULL,
    .ref_fname = NULL,
    .n = 0,
    .mm = NULL,
    .pwm_cols = NULL,
    .motif_min = -20,
    //.fai = NULL,
    .bcf_hdr = NULL,
    .bed = NULL,
    .fp_in = NULL,
    .fp_out = NULL,
    .pcs_col = NULL,
    .n_thread = 5,
    .n_record = 1000,
};
int parse_args(int argc, char **argv)
{
    int i;
    const char *output_format = NULL;
    const char *record = NULL;
    const char *thread = NULL;
    const char *motif_min = NULL;
    
    for ( i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        if ( strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0 ) return usage();

        if ( strcmp(a, "-bed") == 0 )
            var = &args.bed_fname;
        else if ( strcmp(a, "-vcf") == 0 )            
            var = &args.input_fname;
        else if ( strcmp(a, "-motif") == 0 )
            var = &args.motif_fname;
        else if ( strcmp(a, "-O") == 0 )
            var = &output_format;
        else if ( strcmp(a, "-o") == 0 )
            var = &args.output_fname;
        else if ( strcmp(a, "-ref") == 0 )
            var = &args.ref_fname;
        else if ( strcmp(a, "-t") == 0 )
            var = &thread;
        else if ( strcmp(a, "-record") == 0 )
            var = &record;
        else if ( strcmp(a, "-min") == 0 )
            var = &motif_min;
        
        if ( var != 0 ) {
            if ( i == argc ) error("Missing an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        error("Unknown argument: %s, use -h see help information", a);
    }

    if ( args.bed_fname == NULL ) error("Parameter -bed is required.");
    if ( args.input_fname == NULL ) error("Parameter -vcf is required.");
    if ( args.motif_fname == NULL ) error("Parameter -motif is required.");
    if ( args.ref_fname == NULL ) error("Parameter -ref is required.");
    
    args.bed = bedaux_init();
    bed_read(args.bed, args.bed_fname);
    if (args.bed->flag & bed_bit_empty ) error("Cannot load BED file. %s", args.bed_fname);
    bed_merge(args.bed);

    args.fp_in = hts_open(args.input_fname, "r");
    if ( args.fp_in == NULL ) error("%s : %s.", args.input_fname, strerror(errno));

    args.mm = motif_read(args.motif_fname, &args.n);
    if ( args.n == 0 ) error("No motif records.");

    //args.fai = fai_load(args.ref_fname);
    //if ( args.fai == NULL ) error("Failed to load index of %s.", args.ref_fname);

    if ( thread ) args.n_thread = str2int((char*)thread);
    if ( record ) args.n_record = str2int((char*)record);
    if ( motif_min ) args.motif_min = str2int((char*)motif_min);
    int out_type = FT_VCF;
    if ( output_format != 0 ) {
	switch (output_format[0]) {
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

    args.fp_out = args.output_fname == NULL ? hts_open("-", hts_bcf_wmode(out_type)) :
        hts_open(args.output_fname, hts_bcf_wmode(out_type));
    
    args.bcf_hdr = bcf_hdr_read(args.fp_in);
    
#define BRANCH(_key, _description, _num) do {                           \
        int id;                                                         \
        id = bcf_hdr_id2int(args.bcf_hdr, BCF_DT_ID, _key);             \
        if (id == -1) {                                                 \
            bcf_hdr_append(args.bcf_hdr, _description);                 \
            bcf_hdr_sync(args.bcf_hdr);                                 \
            id = bcf_hdr_id2int(args.bcf_hdr, BCF_DT_ID, _key);         \
            assert(bcf_hdr_idinfo_exists(args.bcf_hdr, BCF_HL_INFO, id)); \
            _num = bcf_hdr_id2length(args.bcf_hdr, BCF_HL_INFO, id);    \
        }                                                               \
    } while(0)

    args.pwm_cols = malloc(args.n *sizeof(struct anno_col));
    
    for ( i = 0; i < args.n; ++i ) {
        struct motif *m = args.mm[i];
        kstring_t str = {0,0,0};
        kputs(m->name, &str);
        kputs("_PWM_score", &str);
        struct anno_col *col = &args.pwm_cols[i];
        col->hdr_key = strdup(str.s);
        col->func.pwm = anno_motif_setter_info_float;
        // col->replace = REPLACE_MISSING;
        str.l = 0;
        ksprintf(&str, "##INFO=<ID=%s,Number=R,Type=Float,Description=\"%s PWM score for each allele.\"", col->hdr_key, m->name);
        //debug_print("%s", str.s);
        BRANCH(col->hdr_key, str.s, col->number);
        free(str.s);
    }

    // Pwm Change Score column
    args.pcs_col = malloc(sizeof(struct anno_col));
    args.pcs_col->hdr_key = strdup("PWM_score_change");
    args.pcs_col->func.pwm = anno_motif_setter_info_float;
    // args.pcs_col->replace = REPLACE_MISSING;
    BRANCH(args.pcs_col->hdr_key, "##INFO=<ID=PWM_score_change,Number=1,Type=Float,Description=\"PWM score changes.\"", args.pcs_col->number);

    
#undef BRANCH
    
    
    bcf_hdr_write(args.fp_out, args.bcf_hdr);
    
    return 0;

}
void memory_release()
{
    hts_close(args.fp_in);
    hts_close(args.fp_out);
    int i;
    for ( i = 0; i < args.n; ++i ) {
        motif_destroy(args.mm[i]);
        free(args.pwm_cols[i].hdr_key);
    }
    free(args.mm);
    free(args.pcs_col->hdr_key);
    
    bcf_hdr_destroy(args.bcf_hdr);
    bed_destroy(args.bed);    
}

void *anno_pwm(void *arg, int idx)
{
    struct anno_pool *pool = (struct anno_pool*)arg;
    struct MTF **mm = (struct MTF**)pool->arg;
    struct MTF *m = mm[idx];
    int i;
    for ( i = 0; i < pool->n_reader; ++i ) {
        bcf1_t *l = pool->readers[i];
        if ( l->rid == -1 ) continue;
        if ( bcf_get_variant_types(l) == VCF_REF ) continue;
        if ( MTF_vcf_sync(m, l) == 0 ) continue;
        anno_vcf_motif_pwm(m, l);
    }
    return pool;
}

#include "anno_thread_pool.h"


int main(int argc, char **argv)
{

    if ( parse_args(argc, argv) ) return 1;
    if ( args.n_thread > 1 ) {
        struct MTF **M = calloc(args.n_thread, sizeof(void*));
        int i;
        // initize MTF
        for ( i = 0; i < args.n_thread; ++i ) {
            M[i] = MTF_init();
            M[i]->fai = fai_load(args.ref_fname);
            if (M[i]->fai == NULL ) error("Failed to load index of %s.", args.ref_fname);
            M[i]->n = args.n;
            M[i]->mm = args.mm;
            M[i]->bcf_hdr = args.bcf_hdr;
            M[i]->bed = args.bed;
            M[i]->cols = args.pwm_cols;
            M[i]->ccol = args.pcs_col;
            M[i]->min = args.motif_min;
        }
        
        struct thread_pool *p = thread_pool_init(args.n_thread);
        struct thread_pool_process *q = thread_pool_process_init(p, args.n_thread*2, 0);
        struct thread_pool_result *r;
        
        for (;;) {
            
            struct anno_pool *arg = anno_reader(args.fp_in, args.bcf_hdr, args.n_record);
            
            if ( arg->n_reader == 0 )
                break;
            
            arg->arg = M;
            
            int block;
            do {
                block = thread_pool_dispatch2(p, q, anno_pwm, arg, 1);
                if ( ( r = thread_pool_next_result(q) ) ) {
                    // generate output
                    struct anno_pool *data = (struct anno_pool*)r->data;
                    int i;
                    for ( i = 0; i < data->n_reader; ++i ) {
                        bcf_write1(args.fp_out, args.bcf_hdr, data->readers[i]);
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
                bcf_write1(args.fp_out, args.bcf_hdr, data->readers[i]);
                bcf_destroy(data->readers[i]);
            }
            free(data->readers);
            thread_pool_delete_result(r, 1);
        }
        thread_pool_process_destroy(q);
        thread_pool_destroy(p);
        for ( i = 0; i < args.n_thread; ++i ) MTF_destory(M[i]);
        free(M);
    }
    else {                           
        // initise MTF for each thread
        struct MTF *MTF = MTF_init();
        MTF->fai = fai_load(args.ref_fname);
        if (MTF->fai == NULL ) error("Failed to load index of %s.", args.ref_fname);

        MTF->n = args.n;
        MTF->mm = args.mm;
        MTF->bcf_hdr = args.bcf_hdr;
        MTF->bed = args.bed;
        MTF->cols = args.pwm_cols;
        MTF->ccol = args.pcs_col;
        MTF->min = args.motif_min;
        
        bcf1_t *line = bcf_init();

        for ( ;; ) {
            if ( bcf_read(args.fp_in, args.bcf_hdr, line) != 0 ) break;
            if ( line->rid == -1 ) continue;
            if ( bcf_get_variant_types(line) == VCF_REF ) continue;
            if ( MTF_vcf_sync(MTF, line) == 0 ) continue;
            anno_vcf_motif_pwm(MTF, line);
            bcf_write1(args.fp_out, args.bcf_hdr, line);
        }
        bcf_destroy(line);
        MTF_destory(MTF);
    }
    
    memory_release();
    return 0;
}
