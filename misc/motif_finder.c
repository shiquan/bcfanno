#include "utils.h"
#include "number.h"
#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "motif.h"
#include <zlib.h>
#include <string.h>

int usage()
{
    fprintf(stderr, "motif_finder [options] ref.fa motif_sequence\n");
    fprintf(stderr, " -m [0]  Mismatches.\n");
    fprintf(stderr, " -rev    Count reverse strand also.\n");
    fprintf(stderr, " -bed    Region to scan in BED format.\n");
    fprintf(stderr, " \n");
    return 1;
}
struct args {
    const char *fasta_fname;
    const char *motif;
    const char *bed_fname;
    int mis;
    int rev;
    unsigned char *tab_encode;
    unsigned char *rev_tab_encode;
    struct encode *motif_encode;
    struct encode *motif_rev_encode;
} args = {
    .fasta_fname = NULL,
    .motif = NULL,
    .bed_fname = NULL,
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
        if ( strcmp(a, "-bed") == 0 ) var = &args.bed_fname;
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
