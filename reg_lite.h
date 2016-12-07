#ifndef REG_LITE_HEADER
#define REG_LITE_HEADER


struct reg_lite {
    int start;
    int end;
    struct reg_lite *next;
};

struct reg_spec {
    int chrom_id;
    int n_regs;
    int m_regs;    
    int start;
    int end;
    struct reg_lite *head;
    struct reg_lite *tail;
};

extern struct reg_spec * reg_spec_init();
extern void push_new_reg(struct reg_spec *spec, int start, int end);
extern void reg_set_id(struct reg_spec *spec, int id);
extern void reg_spec_destroy();

#endif
