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

#ifndef ANNO_VCF_HEADER
#define ANNO_VCF_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>
#include <htslib/tbx.h>
#include "anno.h"
#include "vcmp.h"

struct anno_vcf_file {
    int id; // index
    char *fname;
    htsFile *fp;
    // header of vcf
    bcf_hdr_t *hdr;
    // index for bcf file, pre-indexed by bcftools
    hts_idx_t *bcf_idx;
    // index for vcf file, pre-indexed by tabix
    tbx_t *tbx_idx;
    // iterator
    hts_itr_t *itr;
    int no_such_chrom; // flag to set no this chromosome in the database
    int last_rid;
    // char *fname;   
    int cached, max;
    bcf1_t **buffer;
    // char *columns;
    int ncols;
    struct anno_col *cols;
};

struct vcfs_options {
    // this flag should be set 1 if vcf databases are inited, else set 0
    int vcfs_is_inited;
    // hdr is header struct of input vcf file, all the annotated tags should be inited in the hdr_out before export bcf lines
    bcf_hdr_t *hdr_out;
    // pre-indexed vcf/bcf databases
    // columns structure of vcf databases
    // struct anno_cols *cols;    
    // int i_data;
    int n_files;
    int m_files;
    struct anno_vcf_file *files;
    // for matching annotation and VCF lines by allele
    vcmp_t *vcmp;  
    int *sample_map, nsample_map, sample_is_file;   
    int mtmpi, mtmpf, mtmps;
    int mtmpi2, mtmpf2, mtmps2;
    int mtmpi3, mtmpf3, mtmps3;
    int32_t *tmpi, *tmpi2, *tmpi3;
    float *tmpf, *tmpf2, *tmpf3;
    char *tmps, *tmps2, **tmpp, **tmpp2;
    kstring_t tmpks;
};
// vcfs_options_init will init options but will not allocate memory for vcfs_options 
extern int vcfs_options_init(struct vcfs_options *opts);
// clean all options in vcfs_options, but will not free the memory of vcfs_options
extern int vcfs_options_destroy(struct vcfs_options *opts);
// add new vcf/bcf database and init columns
extern int vcfs_database_add(struct vcfs_options *opts, const char *fname, char *columns);
// core function to annotate
extern int anno_vcfs_core(struct vcfs_options *opts, bcf1_t *line);

#endif
