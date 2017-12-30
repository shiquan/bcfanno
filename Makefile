PROG=       bcfanno vcf2tsv tsv2vcf vcf_rename_tags
DEBUG_PROG= bcfanno_debug

all: $(PROG)

debug: $(DEBUG_PROG)

HTSDIR = htslib-1.6
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a

CC       = gcc
CFLAGS   = -Wall -O3 
DEBUG_CFLAGS   = -g -Wall -O0 -DDEBUG_MODE
DFLAGS   =
INCLUDES = -I src2/ -I. -I$(HTSDIR)/ -I misc/


all:$(PROG)

# See htslib/Makefile
PACKAGE_VERSION := $(shell git describe --tags)

version.h:
	echo '#define BCFANNO_VERSION "$(PACKAGE_VERSION)"' > $@


.SUFFIXES:.c .o
.PHONY:all clean clean-all clean-plugins distclean install lib tags test testclean force plugins docs

force:

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@


GenePredExtGen: $(HTSLIB)
	$(CC)  $(DEBUG_CFLAGS) $(INCLUDES) -pthread -lz -o $@ misc/genepred_ext_gen.c misc/ksw.c src2/anno_thread_pool.c src2/genepred.c src2/number.c src2/faidx_def.c $(HTSLIB)

vcf2tsv: $(HTSLIB) version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ misc/vcf2tsv.c $(HTSLIB)

tsv2vcf: $(HTSLIB) version.h
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ misc/tsv2vcf.c misc/table2hash.c $(HTSLIB)

vcf_rename_tags: $(HTSLIB) version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ misc/vcf_rename_tags.c $(HTSLIB)

bcfanno: $(HTSLIB) version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -lz -pthread -o $@ src2/anno_bed.c src2/anno_col.c src2/anno_hgvs.c src2/anno_pool.c src2/anno_thread_pool.c src2/anno_vcf.c src2/bcfanno_main.c src2/config.c src2/flank_seq.c src2/genepred.c src2/hgvs.c src2/json_config.c src2/kson.c src2/name_list.c src2/number.c src2/sort_list.c src2/variant_type.c src2/vcf_annos.c src2/vcmp.c $(HTSLIB)

bcfanno_debug: $(HTSLIB) version.h
	$(CC) -DDEBUG_MODE $(DEBUG_CFLAGS) $(INCLUDES) -lz -pthread -o $@ src2/anno_bed.c src2/anno_col.c src2/anno_hgvs.c src2/anno_pool.c src2/anno_thread_pool.c src2/anno_vcf.c src2/bcfanno_main.c src2/config.c src2/flank_seq.c src2/genepred.c src2/hgvs.c src2/json_config.c src2/kson.c src2/name_list.c src2/number.c src2/sort_list.c src2/variant_type.c src2/vcf_annos.c src2/vcmp.c $(HTSLIB)


test: $(HTSLIB) version.h

clean: testclean
	-rm -f gmon.out *.o *~ $(PROG) version.h 
	-rm -rf *.dSYM plugins/*.dSYM test/*.dSYM
	-rm -f anno_vcf bedadd vcfadd bcfanno anno_bed hgvs_generate hgvs_vcf GenePredExtGen
	-rm -f config bcfanno_debug vcf2tsv tsv2vcf vcf_rename_tags

testclean:
	-rm -f test/*.o test/*~ $(TEST_PROG)

distclean: clean
	-rm -f TAGS

clean-all: clean clean-htslib

tags:
	ctags -f TAGS *.[ch] plugins/*.[ch]
