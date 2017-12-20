PROG=       bcfanno
DEBUG_PROG= bcfanno_debug

all: $(PROG)

debug: $(DEBUG_PROG)

HTSDIR = ../htslib-1.6
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a

CC       = gcc
CFLAGS   = -Wall -O3 
DEBUG_CFLAGS   = -g -Wall -O0 -DDEBUG_MODE
DFLAGS   =
INCLUDES = -I. -I$(HTSDIR)/


all:$(PROG)

PACKAGE_VERSION := $(shell git describe --tags)
version.h:
	echo '#define BCFANNO_VERSION "$(PACKAGE_VERSION)"' > $@


.SUFFIXES:.c .o
.PHONY:all clean clean-all clean-plugins distclean install lib tags test testclean force plugins docs

force:

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@


bcfanno: $(HTSLIB) version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -lz -pthread -o $@  anno_bed.c anno_col.c anno_hgvs.c anno_pool.c anno_thread_pool.c anno_vcf.c bcfanno_main.c config.c genepred.c hgvs.c json_config.c kson.c name_list.c number.c sort_list.c variant_type.c vcf_annos.c vcmp.c flank_seq.c $(HTSLIB)

bcfanno_debug: $(HTSLIB) version.h
	$(CC) -DDEBUG_MODE $(DEBUG_CFLAGS) $(INCLUDES) -lz -pthread -o $@  anno_bed.c anno_col.c anno_hgvs.c anno_pool.c anno_thread_pool.c anno_vcf.c bcfanno_main.c config.c genepred.c hgvs.c json_config.c kson.c name_list.c number.c sort_list.c variant_type.c vcf_annos.c vcmp.c flank_seq.c $(HTSLIB)

test: $(HTSLIB) version.h

clean: testclean
	-rm -f gmon.out *.o *~ $(PROG) version.h 
	-rm -rf *.dSYM plugins/*.dSYM test/*.dSYM
	-rm -f anno_vcf bedadd vcfadd bcfanno anno_bed hgvs_generate hgvs_vcf
	-rm -f config bcfanno_debug vcf2tsv tsv2vcf vcf_rename_tags

testclean:
	-rm -f test/*.o test/*~ $(TEST_PROG)

distclean: clean
	-rm -f TAGS

clean-all: clean clean-htslib

tags:
	ctags -f TAGS *.[ch] plugins/*.[ch]
