/*  
    Copyright (C) 2016,2017,2018  BGI Research

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

#ifndef UTILS_COMMON_HEADER
#define UTILS_COMMON_HEADER
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>
#include <time.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_GRAY    "\x1b[37m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define check_double_free(p) do {\
        void **_pp = (void**)&(p);                                      \
        if (_pp==NULL || *_pp == NULL) {				\
	    fprintf(stderr,"[error] double free %s, %d", __FUNCTION__, __LINE__); \
	    exit(EXIT_FAILURE);						\
	}								\
    } while(0)

#define ignore_free(p) do \
	void **_pp = (void**)&(p);					\
	if (*_pp!=NULL && _pp != NULL) {				\
	    free(*_pp);							\
	    *_pp = NULL;						\
	}								\
    } while(0)

#define safe_free(p) do				\
    {						\
	void **_pp = (void**)&(p);					\
        if (_pp==NULL || *_pp == NULL) {				\
	    fprintf(stderr,"[error] double free %s, %d", __FUNCTION__, __LINE__); \
	    exit(EXIT_FAILURE);						\
	}								\
	free(*_pp);							\
        *_pp = NULL;                                                    \
    } while(0)

#define check_mem(p) do				\
    {						\
	void **_pp = (void**)&p;		\
	if (_pp == NULL || *_pp == NULL) {				\
	    fprintf(stderr, "[memory out] func: %s, line: %d\n", __FUNCTION__, __LINE__);\
	    exit(EXIT_FAILURE);						\
	}								\
    }while(0)

#define str_errno() (errno == 0 ? "None" : strerror(errno))

#define clear_errno() do \
    {\
	fprintf(stderr, "%s\n", str_errno());\
	errno = 0;\
    }while(0)

#define error(line, ...) do						\
    {									\
	fprintf(stderr, ANSI_COLOR_RED "[error] [func: %s, line: %d] " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA line ANSI_COLOR_RESET "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__); \
	errno = 0;							\
	exit(EXIT_FAILURE);						\
    }while(0)

#define error_return(line, ...) do						\
    {									\
	fprintf(stderr, ANSI_COLOR_RED "[error] [func: %s, line: %d] " ANSI_COLOR_RESET ANSI_COLOR_MAGENTA line ANSI_COLOR_RESET "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__); \
	errno = 0;							\
    }while(0)

#define error_print(line, ...) do						\
    {									\
	fprintf(stderr, "[error] [func: %s, line: %d] " line "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__); \
    }while(0)

#define warnings(line, ...) do						\
    {									\
	fprintf(stderr, "[warnings] " line "\n", ##__VA_ARGS__);	\
    }while(0)

#define debug_print(line, ...) do {\
	fprintf(stderr, "[ ** DEBUG ** func: %s, line: %d ] " line "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__); \
    } while(0)

#define LOG_print(line, ...) do {\
	time_t second;\
	time(&second);\
	char _time_buff[100];							\
	strftime (_time_buff, 100, "%Y-%m-%d %H:%M:%S", localtime (&second));	\
	fprintf(stderr, "[%s] " ANSI_COLOR_GREEN line ANSI_COLOR_RESET"\n", _time_buff, ##__VA_ARGS__); \
    } while(0)

//#define DONOT_POST_ERR_STRING "Do NOT post this error message on forums or emails. Please read our online manual. Thank you!"

static inline void *bcfanno_realloc(void *x, size_t size)
{
    void *y = NULL;
    if ( x == NULL ) 
        y = malloc(size);
    else
        y = realloc(x, size);
    
    if ( y == NULL ) {
        free(x);
        return NULL;
    }    
    return y;
}

#endif
