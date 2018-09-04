/*  
    Copyright (C) 2016,2017  BGI Research

    Author: Shi Quan (shiquan@genomics.cn)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE. 
*/

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

#define SEQ_COMP(a,b) (a + b == 3)

typedef char * (*func_dup_seq)(const char *, unsigned long );

extern char *rev_seqs(const char *dna_seqs, unsigned long n);
extern int check_stop_codon(char *seq, char *p_end, int mito);
extern void compl_seq(char *seq, int l);
extern int seq2code4(int seq);
extern int same_DNA_seqs(const char *a, const char *b, int l );
extern int codon2aminoid(char *codon,int mito);
extern char *rev_seqs(const char *dna_seqs, unsigned long n);

#define X_CODO   0

#define C4_Stop  0
#define C4_Phe   1
#define C4_Leu   2
#define C4_Ser   3
#define C4_Tyr   4
#define C4_Cys   5
#define C4_Trp   6
#define C4_Pro   7
#define C4_His   8
#define C4_Gln   9
#define C4_Arg  10
#define C4_Ile  11
#define C4_Met  12
#define C4_Thr  13
#define C4_Asn  14
#define C4_Lys  15
#define C4_Val  16
#define C4_Ala  17
#define C4_Asp  18
#define C4_Glu  19
#define C4_Gly  20

const static char *codon_names[] = {
    "*", "Phe", "Leu", "Ser", "Tyr", "Cys", "Trp", "Pro", "His", "Gln", "Arg", "Ile", "Met", "Thr", "Asn", "Lys", "Val", "Ala", "Asp", "Glu", "Gly",
};

const static char *codon_short_names[] = { "*", "F", "L", "S", "Y", "C", "W", "P", "H", "Q", "R", "I", "M", "T", "N", "K", "V", "A", "D", "E", "G", };

const static int codon_matrix[4][4][4] = {
    { { C4_Lys, C4_Asn, C4_Lys, C4_Asn, },
      { C4_Thr, C4_Thr, C4_Thr, C4_Thr, },
      { C4_Arg, C4_Ser, C4_Arg, C4_Ser, },
      { C4_Ile, C4_Ile, C4_Met, C4_Ile, }, },
    
    { { C4_Gln, C4_His, C4_Gln, C4_His, },
      { C4_Pro, C4_Pro, C4_Pro, C4_Pro, },
      { C4_Arg, C4_Arg, C4_Arg, C4_Arg, },
      { C4_Leu, C4_Leu, C4_Leu, C4_Leu, }, },
    
    { { C4_Glu, C4_Asp, C4_Glu, C4_Asp, },
      { C4_Ala, C4_Ala, C4_Ala, C4_Ala, },
      { C4_Gly, C4_Gly, C4_Gly, C4_Gly, },
      { C4_Val, C4_Val, C4_Val, C4_Val, }, },
    
    { { C4_Stop, C4_Tyr, C4_Stop, C4_Tyr, },
      { C4_Ser, C4_Ser, C4_Ser, C4_Ser, },
      { C4_Stop, C4_Cys, C4_Trp, C4_Cys, },
      { C4_Leu, C4_Phe, C4_Leu, C4_Phe, }, },
};

/*
  https://www.mitomap.org/foswiki/bin/view/MITOMAP/HumanMitoCode

  For human mitochondrial genes, unlike the universal code, UGA codes for tryptophan instead of termination and AUA codes for methionine instead of isoleucine.
 */

const static int mitomap_codon_matrix[4][4][4] = {
    { { C4_Lys, C4_Asn, C4_Lys, C4_Asn, },
      { C4_Thr, C4_Thr, C4_Thr, C4_Thr, },
      { C4_Arg, C4_Ser, C4_Arg, C4_Ser, },
      { C4_Met, C4_Ile, C4_Met, C4_Ile, }, },
    
    { { C4_Gln, C4_His, C4_Gln, C4_His, },
      { C4_Pro, C4_Pro, C4_Pro, C4_Pro, },
      { C4_Arg, C4_Arg, C4_Arg, C4_Arg, },
      { C4_Leu, C4_Leu, C4_Leu, C4_Leu, }, },
    
    { { C4_Glu, C4_Asp, C4_Glu, C4_Asp, },
      { C4_Ala, C4_Ala, C4_Ala, C4_Ala, },
      { C4_Gly, C4_Gly, C4_Gly, C4_Gly, },
      { C4_Val, C4_Val, C4_Val, C4_Val, }, },
    
    { { C4_Stop, C4_Tyr, C4_Stop, C4_Tyr, },
      { C4_Ser, C4_Ser, C4_Ser, C4_Ser, },
      { C4_Trp, C4_Cys, C4_Trp, C4_Cys, },
      { C4_Leu, C4_Phe, C4_Leu, C4_Phe, }, },
};


// DNA level
enum variant_type {
    var_type_nonref = -1, // for gatk <NONREF> allele
    var_type_unknow = 0,
    var_type_ref,
    var_type_snp,
    var_type_del,
    var_type_ins,
    var_type_delins,
    var_type_copy, // dup for ins
    var_type_complex, 
};

// RNA level
// check the variants type
enum var_func_type {
    _var_type_promoter_to_int = -1,
    var_is_unknown,
    var_is_reference,
    var_is_intron,
    var_is_noncoding,
    var_is_utr5,
    var_is_utr3,
    var_is_synonymous,
    var_is_missense,
    var_is_nonsense, // stop gained
    var_is_inframe_insertion,
    var_is_inframe_deletion,
    var_is_inframe_delins,
    var_is_frameshift,
    var_is_stop_lost,
    var_is_stop_retained,
    var_is_complex,
    var_is_no_call,
    var_is_transcript_ablation, // whole exome deletion
    var_is_start_lost,
};

enum var_type_splice {
    var_is_not_splice = 0,
    var_is_splice_site,
    var_is_splice_donor,
    var_is_splice_acceptor,    
};

static inline const char *var_func_type_string(enum var_func_type type)
{
    static const char* vartypes[21] = {
        "Unknown",
        "Reference",
        "Intron",
        "Noncoding",
        "Utr5",
        "Utr3",
        "Synonymous",
        "Missense",
        "Nonsense",
        "InframeInsertion",
        "InframeDeletion",
        "InframeDelins",
        "Frameshift",
        "StopLost",
        "StopRetained",
        "Complex",
        "NoCall",
        "TranscriptAblation",
        "StartLost",
        NULL,
        NULL,
    };
    assert(type >= 0);
    return vartypes[type];
}

static inline const char *var_type_splice_string(enum var_type_splice type)
{
    static const char *splicetypes[5] = {
        "NotSplice",
        "SpliceSite",
        "SpliceDonor",
        "SpliceAcceptor",
        NULL,
    };
    assert(type >= 0 );
    return splicetypes[type];
}

// 1 on yes, 0 on no
static inline int check_is_stop(char *codon, int mito)
{
    if ( mito == 0 ) 
        return codon_matrix[seq2code4(codon[0])][seq2code4(codon[1])][seq2code4(codon[2])] == C4_Stop;
    else
        return mitomap_codon_matrix[seq2code4(codon[0])][seq2code4(codon[1])][seq2code4(codon[2])] == C4_Stop;
}

#endif
