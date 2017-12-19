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
#include "sort_list.h"

int count_list(const void *list)
{
    struct list_lite *pl = (struct list_lite *)list;
    int l;
    for ( l = 0; pl; l++ )
        pl = pl->next;
    return l;
}

int sort_list(void *plist, comp_func *func)
{
    struct list_lite **pp = (struct list_lite **)plist;
    struct list_lite *list = *pp;
    int l;
    l = count_list(list);
    // Too short to sort
    if ( l <= 1 )
        return 1;
    struct list_lite *el;
    struct list_lite **array;
    array = (struct list_lite**)malloc(l * sizeof(struct list_lite*));
    int i;
    for ( el = list, i= 0; el != NULL; el = el->next, i++ )
        array[i] = el;
    
    qsort(array, l, sizeof(array[0]), func);

    *pp = array[0];
    for ( i =0; i < l-1; ++i ) {
        array[i]->next = array[i+1];
        array[i+1]->next = NULL;
    }
    free(array);
    return 0;      
}

int sort_rmdup_list(void *plist, comp_func *func, del_func del_func)
{
    struct list_lite **pp = (struct list_lite**)plist;
    struct list_lite *list = *pp;
    int l;
    l = count_list(list);
    if (l <= 1 )
        return 1;
    struct list_lite *el;
    struct list_lite **array;
    array = (struct list_lite**)malloc(l*sizeof(struct list_lite*));
    int i;
    for ( el = list, i = 0; el != NULL; el = el->next, i++)
        array[i] = el;
    qsort(array, l, sizeof(array[0]), func);
    *pp = array[0];
    for ( el = *pp, i = 1; i < l; i++ ) {
        if ( func(array[i-1], array[i]) == 0 ) {
            del_func(array[i]);
        } else {
            el->next = array[i];
            el = array[i];
            el->next = NULL;
        }
    }
    free(array);
    return 0;        
}
void list_lite_del(void *plist, del_func del_func)
{
    struct list_lite **pp = (struct list_lite**)plist;
    struct list_lite *list = *pp;

    while ( *pp ) {
        list = *pp;
        *pp = list->next;
        del_func(list);
    }
}

