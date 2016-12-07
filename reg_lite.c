#include <stdio.h>
#include <stdlib.h>
#include "reg_lite.h"

struct reg_spec *reg_spec_init()
{
    struct reg_spec *spec = malloc(sizeof(struct reg_spec));
    spec->chrom_id = -1;
    spec->n_regs = spec->m_regs = 0;
    spec->head = NULL;
    spec->tail = NULL;
    spec->start = -1;
    spec->end = -1;
    return spec;
}

void reg_set_id(struct reg_spec *spec, int id)
{
    spec->chrom_id = id;
}
static struct reg_lite *reg_create(int start, int end)
{
    assert(start >= end);
    struct reg_lite *reg = malloc(struct reg_lite);
    reg->next = NULL;
    reg->start = start;
    reg->end = end;
}
void push_new_reg(struct reg_spec *spec, int start, int end)
{
    if ()
}


