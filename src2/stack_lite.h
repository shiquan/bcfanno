#ifndef STACK_LITE_H
#define STACK_LITE_H

#include "utils.h"

struct anno_stack {
    int l, m;
    char **a;
};

static inline struct anno_stack *anno_stack_init()
{
    struct anno_stack *s = malloc(sizeof(*s));
    s->m = 0;
    s->l = 0;
    s->a = NULL;
    return s;
}
// return -1 on duplicate
//         0 on top
static inline int anno_stack_push(struct anno_stack *s, char *name)
{
    assert(name);
    int i;
    for ( i = 0; i < s->l; ++i ) {
        if ( strcmp(s->a[i], name) == 0 )
            return -1;
    }
    if ( s->m == s->l ) {
        s->m = s->m == 0 ? 2 : s->m+2;
        s->a = realloc(s->a, sizeof(void*)*s->m);
    }
    s->a[s->l++] = strdup(name);
    return 0;
}
static inline void anno_stack_destroy(struct anno_stack *s)
{
    int i;
    for ( i = 0; i < s->l; ++i )
        free(s->a[i]);
    free(s->a);
    free(s);
}

#endif
