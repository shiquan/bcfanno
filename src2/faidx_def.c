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

#include "utils.h"
#include "number.h"
#include "htslib/bgzf.h"
#include "htslib/khash.h"

struct faidx_val {
    int32_t line_len, line_blen;
    int64_t len;
    uint64_t offset;
};

KHASH_MAP_INIT_STR(s, struct faidx_val)

struct faidx {
    BGZF *bgzf;
    int n, m;
    char **name;
    khash_t(s) *hash;
};

// trans_retrieve designed to retrieve sequences and version from refMrna.fa.gz database.
// the header format of each sequence in refMrna.fa.gz should be defined as ">trans version"
int trans_retrieve_version(void *_fai, const char *trans)
{
    int l;
    int length;
    int version;

    // if transcript name already have version number, parse version and return it
    length = strlen(trans);
    for ( l = 0; l < length; ++l )
        if ( trans[l] == '.') break;
    l++;
    if ( l < length) {
        version = str2int_l((char*)trans+l, length-l);
        if ( version > 0 )
            return version;
    }
    
    // if no version number in the transcript name try to find version in the transcript fasta
    struct faidx *fai = (struct faidx*)_fai;
    khint_t iter;

    iter = kh_get(s, fai->hash, trans);
    if ( iter == kh_end(fai->hash) ) {
        warnings("Reference %s not fount in FASTA file, returing -1 version.", trans);
        return -1;
    }

    struct faidx_val *val = &kh_val(fai->hash, iter);
    char s[100];
    uint64_t offset = 0;
    if ( val->offset > 100 )
        offset = val->offset - 100;

    int ret;
    ret = bgzf_useek(fai->bgzf, offset, SEEK_SET);
    if ( ret < 0 )
        error("Seeking is a compressed, .gzi unindexed, file?");
    uint64_t i;
    int j = 0;
    int c;
    for ( i = offset; i < val->offset && (c=bgzf_getc(fai->bgzf))>= 0; i++ ) {
        s[j++] = c;
    }
    if ( c < 0 )
        error("Error reading file.");

    s[j] = '\0';
    for ( l = j-2; l >= 0; --l ) {
        if ( s[l] == '>') {
            int k;
            for ( k = l+1; k < j-2; k++)
                if ( s[k] == ' ' )
                    break;
            s[k] = '\0';
            if ( strncmp(s+l+1, trans, k - l -1) != 0 )
                error("Transcript inconsistant. %s vs %s", s+i+1, trans);
            for ( l = k + 1; l < j && s[l] != '\n'; ++l);
            version = str2int_l(s+k+1, l-k-1);
            break;
        }
    }

    if ( version > 0 )
        return version;

    // failed to parse version
    return -1;
}

#ifdef FAIDX_DEF_MAIN

#include "htslib/faidx.h"

struct args {
    faidx_t *fai;
    const char *name;
} args = {
    .fai = NULL,
    .name = 0,
};

int parse_args(int ac, char **av)
{
    if ( ac != 3 )
        error("Usage: retrieve_version refMrna.fa.gz transcript");

    args.fai = fai_load(av[1]);
    if ( args.fai == NULL )
        error("Failed to load index file of %s", av[1]);

    args.name = av[2];
    return 0;
}

int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    int version;
    version = trans_retrieve_version((void*)args.fai, args.name);
    printf("%s : %d\n", args.name, version);
    fai_destroy(args.fai);
    return 0;
}

#endif
