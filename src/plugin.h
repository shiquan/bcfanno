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

#ifndef BCFANNO_PLUGIN_HEADER
#define BCFANNO_PLUGIN_HEADER

#include <stdio.h>
#include <stdlib.h>

struct pl_header {
    int n_records;
    char **records;
};

// replace mark for VCF tag
#ifndef REPLACE_MISSING
#define REPLACE_MISSING  0 // replace only missing values
#endif

#ifndef REPLACE_ALL
#define REPLACE_ALL      1 // replace both missing and existing values
#endif

#ifndef REPLACE_EXISTING
#define REPLACE_EXISTING 2 // replace only if tgt is not missing
#endif

#ifndef SET_OR_APPEND
#define SET_OR_APPEND    3 // set new value if missing or non-existent, append otherwise
#endif

enum info_type {
    info_val_is_flag = 0,
    info_val_is_float,
    info_val_is_int,
    info_val_is_str,
};

struct pl_column {
    int number;
    int replace;
    char *key;
    enum info_type type;
    int n;
    union {
        int32_t *i;
        float *f;
        char *str;
    } value;
};

struct module_data {
    int n_cols;
    struct pl_column *data;
};

#endif
