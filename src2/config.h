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

#ifndef BCFANNO_CONFIG_H
#define BCFANNO_CONFIG_H

#include <stdio.h>
#include <stdlib.h>

struct refgene_config {
    // if refgene_is_set == 1, genepred_fname is mandatory.
    int refgene_is_set;
    char *genepred_fname;
    char *refseq_fname;
    char *trans_list_fname;
    char *gene_list_fname;
    char *columns;
};

struct file_config {
    // file path
    char *fname;
    // columns string
    char *columns;
};
struct vcf_config {
    // vcf files number
    int n_vcf;
    struct file_config *files;
};
struct bed_config {
    // bed files number
    int n_bed;
    struct file_config *files;
};

struct module_config {
    int n_module;
    struct file_config *files;
};

// skip other keys except author, config_id and reference_version
struct bcfanno_config {
    char *author;
    char *config_id;
    char *reference_version;
    char *reference_path;
    struct vcf_config vcf;
    struct bed_config bed;
    struct refgene_config refgene;
    struct module_config module;
};

extern struct bcfanno_config *bcfanno_config_init(void);

extern void bcfanno_config_destroy(struct bcfanno_config *);

extern int bcfanno_load_config(struct bcfanno_config *, const char *);

extern int bcfanno_config_debug(struct bcfanno_config *config);

#endif
