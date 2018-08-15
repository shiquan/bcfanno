#include "utils.h"
#include "number.h"
#include "htslib/bgzf.h"
#include <zlib.h>
#include <string.h>

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

int usage()
{
    fprintf(stderr, "motif_finder [options] ref.fa motif_sequence\n");
    fprintf(stderr, " -m [0]  Mismatches.\n");
    fprintf(stderr, " -rev    Count reverse strand also.\n");
    fprintf(stderr, " \n");
    return 1;
}
unsigned char *base_tab_init()
{
    unsigned char *t = malloc(sizeof(char)*256);
    memset(t, 0, sizeof(char)*256);
    t['A'] = BASEA;
    t['C'] = BASEC;
    t['G'] = BASEG;
    t['T'] = BASET; 
    t['a'] = BASEA;
    t['c'] = BASEC;
    t['g'] = BASEG;
    t['t'] = BASET;
    t['R'] = BASER;
    t['r'] = BASER;
    t['Y'] = BASEY;
    t['y'] = BASEY;
    t['S'] = BASES;
    t['s'] = BASES;
    t['K'] = BASEK;
    t['k'] = BASEK;
    t['M'] = BASEM;
    t['m'] = BASEM;
    t['B'] = BASEB;
    t['b'] = BASEB;
    t['D'] = BASED;
    t['d'] = BASED;
    t['H'] = BASEH;
    t['h'] = BASEH;
    t['V'] = BASEV;
    t['v'] = BASEV;
    return t;
}
unsigned char *base_tab_rev_init()
{
    unsigned char *t = malloc(sizeof(char)*256);
    memset(t, 0, sizeof(char)*256);
    t['A'] = BASET;
    t['C'] = BASEG;
    t['G'] = BASEC;
    t['T'] = BASEA; 
    t['a'] = BASET;
    t['c'] = BASEG;
    t['g'] = BASEC;
    t['t'] = BASEA;

    t['R'] = BASEY;
    t['r'] = BASEY;
    t['Y'] = BASER;
    t['y'] = BASER;
    t['S'] = BASES;
    t['s'] = BASES;
    t['K'] = BASEM;
    t['k'] = BASEM;
    t['M'] = BASEK;
    t['m'] = BASEK;
    t['B'] = BASEV;
    t['b'] = BASEV;
    t['D'] = BASEH;
    t['d'] = BASEH;
    t['H'] = BASED;
    t['h'] = BASED;
    t['V'] = BASEB;
    t['v'] = BASEB;

    return t;
}
struct encode {
    int l;
    uint64_t x;
};
int countbits(uint64_t x, struct encode *enc)
{
    int i;
    int l = 0;
    x = x & enc->x;
    for ( i = 0; i < enc->l; i++) {
        if ( ((x>>(i*4))&0xf) > 0 ) l++;
    }   
    return l;
}
struct encode *revcode(struct encode *c)
{
    struct encode *r = malloc(sizeof(*r));
    r->l = c->l;
    r->x = 0;
    int i;
    for ( i =0; i < c->l; ++i ) r->x = r->x | ((c->x>>i)&0x1)<<(c->l-i-1);
    return r;
};
struct encode *str_encode(const char *s, unsigned char *tab)
{
    int i;
    int l = strlen(s);
    struct encode *x = malloc(sizeof(*x));
    x->x = 0;
    x->l = 0;
    // only encode first 16 bases or shorter
    x->l = l > 16 ? 16 : l;
    //uint64_t mask = (1ULL<<l*4)-1;
    for (i = 0; i < x->l; ++i) {
        x->x = x->x<<4 | (tab[s[i]]&0xf);
    }
    return x;
}

struct args {
    const char *fasta_fname;
    const char *motif;
    int mis;
    int rev;
    unsigned char *tab_encode;
    unsigned char *rev_tab_encode;
    struct encode *motif_encode;
    struct encode *motif_rev_encode;
} args = {
    .fasta_fname = NULL,
    .motif = NULL,
    .mis = 0,
    .rev = 0,
    .tab_encode = NULL,
    .rev_tab_encode = NULL,
    .motif_encode = NULL,
    .motif_rev_encode = NULL,
};
int parse_args(int argc, char **argv)
{
    //if ( argc <= 3) return usage();
    int i;
    const char *mis = NULL;

    for ( i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        if ( strcmp(a, "-h") == 0 ) return usage();

        if ( strcmp(a, "-m") == 0 ) var = &mis;
        else if ( strcmp(a, "-rev") == 0 ) {
            args.rev = 1;
            continue;
        }

        if ( var != 0 ) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if ( args.fasta_fname == 0 ) {
            args.fasta_fname = a;
            continue;
        }

        if ( args.motif == 0) {
            args.motif = a;
            continue;
        }
        error("Unknown argument: %s, use -h see help information.", a);
    }

    if ( args.motif == NULL || args.fasta_fname == NULL ) error("Fasta and MOTIF sequence must be set!");
    if ( mis ) args.mis = str2int((char*)mis);
    
    args.tab_encode = base_tab_init();
    args.rev_tab_encode = base_tab_rev_init();

    args.motif_encode = str_encode(args.motif, args.tab_encode);
    args.motif_rev_encode = revcode(args.motif_encode);

    if ( args.motif_encode->l < 5 ) error("MOTIF is too short %s.", args.motif);
    
    if ( args.mis > 0 ) {
        int cf = args.motif_encode->l/10;
        if ( args.mis > cf ) {
            warnings("Cannot tolerate more mismatches than 10 percent of motifs. Trying to adapt mismatches to %d", cf);
            args.mis = cf;
        }
    }

    /* debug_print("%s\t%d\t%u", args.motif, args.motif_encode->l, args.motif_encode->x); */
    /* uint64_t x = args.motif_encode->x; */
    /* for ( i = 0; i < args.motif_encode->l; ++i ) { */
    /*     fprintf(stderr, "%c","NACNGNNNTNN"[x&0xf]); */
    /*     x = x>>4; */
    /* } */
    /* exit(1); */
    
    return 0;
}
void memory_release()
{
    free(args.tab_encode);
    free(args.rev_tab_encode);
    free(args.motif_encode);
    free(args.motif_rev_encode);
}
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv)) return 1;
    
    BGZF *bgzf = bgzf_open(args.fasta_fname, "r");
    if ( bgzf == NULL ) error("%s: %s.", args.fasta_fname, strerror(errno));

    int c;
    int pos = -1;
    uint64_t x = 0;
    unsigned char *tab = args.tab_encode;
    struct encode *enc = args.motif_encode;
    int n = 0;
    int l = 0, m = 2;
    char *name = malloc(2);
    uint64_t mask = (1LL<<(enc->l*4)) -1;
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
                // fprintf(stderr, "%c\t%c\n", c,"NACNGNNNTNN"[x&0xf]);
                if ( ++n >= enc->l ) {
                    int b = countbits(x, enc);
                    if (b >= enc->l-args.mis) {
                        printf("%s\t%d\t%d\n", name, pos-enc->l+1, pos+1);
                    }
                }
            }
        }
    }
    
    free(name);
    bgzf_close(bgzf);
    
    memory_release();
    return 0;
}
