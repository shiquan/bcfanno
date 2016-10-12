// sequence - a simple collection for handling sequences
//

#ifndef SEQUENCE_HEADER
#define SEQUENCE_HEADER

#include <stdio.h>
#include <stdlib.h>

#if defined(_MSC_VER) && !defined(__clang__)
# define inline __inline
#endif

#define C4_A  0
#define C4_C  1
#define C4_G  2
#define C4_T  3
#define C4_U  3
#define C4_N  4
#define seqarr    "ACGTN"
#define revseqarr "TGCAN"

#define SEQ_COMP(a, b) (a + b == 3)

static inline uint8_t seq2code4(uint8_t seq)
{
    static const uint8_t seq2num_table[256] = {
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    };
    return seq2num_table[seq];
}

static inline char *rev_seqs(const char *dna_seqs, unsigned long n)
{
    if ( n == 0 )
        return NULL;
    int i;
    char *rev = (char*)calloc(n+1, sizeof(char));    
    for ( i = 0; i < n; ++i) 
        rev[i] = revseqarr[seq2code4(dna_seqs[n-i-1])];
    rev[n] = '\0';
    return rev;
}


#define C4_Stop 0
#define C4_Phe  1
#define C4_Leu  2
#define C4_Ser  3
#define C4_Tyr  4
#define C4_Cys  5
#define C4_Trp  6
#define C4_Pro  7
#define C4_His  8
#define C4_Gln  9
#define C4_Arg 10
#define C4_Ile 11
#define C4_Met 12
#define C4_Thr 13
#define C4_Asn 14
#define C4_Lys 15
#define C4_Val 16
#define C4_Ala 17
#define C4_Asp 18
#define C4_Glu 19
#define C4_Gly 20

const static char *codon_names[] = {
    "Stop", "Phe", "Leu", "Ser", "Tyr", "Cys", "Trp", "Pro", "His", "Gln",
    "Arg", "Ile", "Met", "Thr", "Asn", "Lys", "Val", "Ala", "Asp", "Glu", "Gly",
};

const static int codon_matrix[4][4][4] = {
    {
        { C4_Lys, C4_Asn, C4_Lys, C4_Asn, },
        { C4_Thr, C4_Thr, C4_Thr, C4_Thr, },
        { C4_Arg, C4_Ser, C4_Arg, C4_Ser, },
        { C4_Ile, C4_Ile, C4_Met, C4_Ile, },
    },{    
        { C4_Gln, C4_His, C4_Gln, C4_His, },
        { C4_Pro, C4_Pro, C4_Pro, C4_Pro, },
        { C4_Arg, C4_Arg, C4_Arg, C4_Arg, },
        { C4_Leu, C4_Leu, C4_Leu, C4_Leu, },
    },{
                                     
        { C4_Glu, C4_Asp, C4_Glu, C4_Asp, },
        { C4_Ala, C4_Ala, C4_Ala, C4_Ala, },
        { C4_Gly, C4_Gly, C4_Gly, C4_Gly, },
        { C4_Val, C4_Val, C4_Val, C4_Val, },
    },{
        { C4_Stop, C4_Tyr, C4_Stop, C4_Tyr, },
        { C4_Ser, C4_Ser, C4_Ser, C4_Ser, },
        { C4_Stop, C4_Cys, C4_Trp, C4_Cys, },
        { C4_Leu, C4_Phe, C4_Leu, C4_Phe, },
    },
};
// no check the codon length for fast read
static inline int codon2aminoid(char *codon)
{
    return codon_matrix[seq2code4(codon[0])][seq2code4(codon[1])][seq2code4(codon[2])];
}

// check the variants type
enum var_type {
    _var_type_promoter_to_int = -1,
    var_is_unknown,
    var_is_reference,
    var_is_intron,
    var_is_noncoding,
    var_is_utr5,
    var_is_utr3,
    var_is_synonymous,
    var_is_nonsynonymous,
    var_is_inframe_insertion,
    var_is_inframe_deletion,
    var_is_stop_gained,
    var_is_frameshift,
    var_is_stop_lost,
    var_is_stop_retained,
    var_is_splice_site,
    var_is_splice_donor,
    var_is_splice_acceptor,
    var_is_complex,
};

static inline const char *var_type_string(enum var_type type)
{
    static const char* vartypes[20] = {
        "unknown",
        "reference",
        "intron",
        "noncoding",
        "utr5",
        "utr3",
        "synonymous",
        "nonsynonymous",
        "inframe insertion",
        "inframe deletion",
        "stop gained",
        "frameshift",
        "stop lost",
        "stop retained",
        "splice site",
        "splice donor",
        "splice acceptor",
        "complex",
        NULL,
        NULL,
    };
    assert(type >= 0);
    return vartypes[type];
}
// 1 on yes, 0 on no
static inline int check_is_stop(char *codon)
{
    return codon_matrix[seq2code4(codon[0])][seq2code4(codon[1])][seq2code4(codon[2])] == C4_Stop;
}

extern enum var_type check_var_type(char *block, int block_length, int start, char *ref, int ref_length, char *alt, int alt_length );


#endif
