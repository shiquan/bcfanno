#include "utils.h"
#include <string.h>
#include <math.h>

#ifndef MOTIF_H
#define MOTIF_H 

#include "htslib/faidx.h"
#include "htslib/vcf.h"
#include "wrap_pileup.h"
#include "htslib/kseq.h"
#include "bed_utils.h"
#include "anno_col.h"

#define BASEA 1
#define BASEC 2
#define BASEG 4
#define BASET 8

/*
  IUPAC nucleotide code
  R	A or G
  Y	C or T
  S	G or C
  W	A or T
  K	G or T
  M	A or C
  B	C or G or T
  D	A or G or T
  H	A or C or T
  V	A or C or G
*/

#define BASER 5
#define BASEY 10
#define BASES 6
#define BASEK 12
#define BASEM (1 | (1<<1))
#define BASEB ((1<<1) | (1<<2) | (1<<3))
#define BASED (1 | (1<<2) | (1<<3))
#define BASEH (1 | (1<<1) | (1<<3))
#define BASEV (1 | (1<<1) | (1<<3))

/*
static inline struct encode *revcode(struct encode *c)
{
    struct encode *r = malloc(sizeof(*r));
    r->l = c->l;
    r->x = 0;
    int i;
    for ( i =0; i < c->l; ++i ) r->x = r->x | ((c->x>>i)&0x1)<<(c->l-i-1);
    return r;
};

static inline struct encode *str_encode(const char *s)
{
    static const int motif_encode_tab[256] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 14, 2, 13, 0, 0, 4, 11, 0, 0, 12, 0, 3, 0, 0,
        0, 0, 5, 6, 8, 0, 11, 0, 0, 10, 0, 0, 0, 0, 0, 0,
        0, 1, 14, 2, 13, 0, 0, 4, 11, 0, 0, 12, 0, 3, 0, 0,
        0, 0, 5, 6, 8, 0, 11, 0, 0, 10, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };
    int i;
    int l = strlen(s);
    struct encode *x = malloc(sizeof(*x));
    x->x = 0;
    x->l = 0;
    // only encode first 16 bases or shorter
    x->l = l > 16 ? 16 : l;
    //uint64_t mask = (1ULL<<l*4)-1;
    for (i = 0; i < x->l; ++i) {
        x->x = x->x<<4 | (motif_encode_tab[s[i]]&0xf);
    }
    return x;
}

static inline struct encode *str_encode_rev(const char *s)
{
    static const int motif_encode_rev_tab[256] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 8, 11, 4, 11, 0, 0, 2, 13, 0, 0, 3, 0, 12, 0, 0,
        0, 0, 10, 6, 1, 0, 14, 0, 0, 5, 0, 0, 0, 0, 0, 0,
        0, 8, 11, 4, 11, 0, 0, 2, 13, 0, 0, 3, 0, 12, 0, 0,
        0, 0, 10, 6, 1, 0, 14, 0, 0, 5, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };

    int i;
    int l = strlen(s);
    struct encode *x = malloc(sizeof(*x));
    x->x = 0;
    x->l = 0;
    // only encode first 16 bases or shorter
    x->l = l > 16 ? 16 : l;
    //uint64_t mask = (1ULL<<l*4)-1;
    for (i = 0; i < x->l; ++i) {
        x->x = x->x<<4 | (motif_encode_rev_tab[s[i]]&0xf);
    }
    return x;
}
*/

/* encode bitcodes */

// encode at maximal 32 bases into bitcodes
struct encode16 {
    int l;
    uint64_t x;
};
// encode at maximal 32 bases into bitcodes
struct encode32 {
    int l;
    uint64_t x[2];
};
// encode at maximal 64 bases into bitcodes
struct encode64 {
    int l;
    uint64_t x[4];
};
// dynamic load encode bitcodes
struct encode {
    int l;
    void *x;
};

extern void encode_destory(struct encode *x);
extern struct encode *str_encode(const char *s);
extern struct encode *str_encode_rev(const char *s);
extern int encode_query(struct encode *q, struct encode *r, int m);
extern int encode_query_seq(struct encode *q, const char *s, int m);


#define MOTIF_INF 0.0000000001
#define MOTIF_INF_PWM -3 // <0.01
#define MOTIF_MAX_PWM 1.386
#define MOTIF_MIN_PWM -20


struct motif {
    char *name;
    struct encode *enc;
    struct encode *rev;
    int m; // allocated memory
    int n; // motif length
    float *map[4];
};

extern struct motif **motif_read(const char *fname, int *_n);

//
//  PWM_score_change
//  CTCF_PWM_score_change
//
struct MTF {
    int n;
    struct motif **mm; // point to args::mm
    int tid;
    int start, end;
    bcf_hdr_t *bcf_hdr; // point to args::bcf_hdr
    struct bedaux *bed; // point to args::bed
    int id; // maximal PWM_score_change
    int *motif_IDs; // point args::motif_IDs
    struct plp_ref *r;
    char *ref;
    int ref_len;
    faidx_t *fai;
    struct anno_col *cols; // point to args::pwm_cols
    struct anno_col *ccol; // point to args::pcs_col

    int min;
};

extern struct MTF *MTF_init();
extern void MTF_destory(struct MTF *m);
extern int anno_vcf_motif_pwm(struct MTF *MTF, bcf1_t *line);


#endif
