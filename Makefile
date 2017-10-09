PROG=       bcfanno vcf2tsv tsv2vcf vcf_rename_tags
DEBUG_PROG= bcfanno_debug

all: $(PROG)

debug: $(DEBUG_PROG)

HTSDIR = htslib-1.3
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a

CC       = gcc
CFLAGS   = -Wall -O3 -DHTS3
DEBUG_CFLAGS   = -g -Wall -O0 -DHTS3 -DDEBUG_MODE
DFLAGS   =
INCLUDES = -I src/ -I. -I$(HTSDIR)/


all:$(PROG)

# See htslib/Makefile
PACKAGE_VERSION = 0.01
ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --tags)
DOC_VERSION :=  $(shell git describe --tags)+
DOC_DATE := $(shell date +'%Y-%m-%d %R %Z')
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif
version.h:
	echo '#define BCFANNO_VERSION "$(PACKAGE_VERSION)"' > $@


.SUFFIXES:.c .o
.PHONY:all clean clean-all clean-plugins distclean install lib tags test testclean force plugins docs

force:

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@


# Plugin rules
#PLUGINC = $(foreach dir, plugins, $(wildcard $(dir)/*.c))
#PLUGINS = $(PLUGINC:.c=.so)
#PLUGINM = $(PLUGINC:.c=.mk)

#%.so: %.c version.h version.c $(HTSDIR)/libhts.so
#	$(CC) $(CFLAGS) $(INCLUDES) -fPIC -shared -o $@ version.c $< -L$(HTSDIR) -lhts

#-include $(PLUGINM)

#plugins: $(PLUGINS)

#hgvs_generate: $(HTSLIB) 
#	$(CC) -D_HGVS_MAIN $(DEBUG_CFLAGS) $(INCLUDES) -pthread -lz -o $@ sequence.c genepred.c hgvs_generate.c $(HTSLIB)

hgvs_vcf: $(HTSLIB)
	$(CC) -DHGVS_VCF_MAIN $(DEBUG_CFLAGS) $(INCLUDES) -pthread -lz -o $@ src/hgvs_vcf.c src/hgvs.c src/sequence.c src/genepred.c src/number.c src/sort_list.c $(HTSLIB)

vcfadd: $(HTSLIB) 
	$(CC) -D_VCF_ANNOS_MAIN $(DEBUG_CFLAGS) $(INCLUDES) -pthread -lz -o $@ src/vcf_annos.c src/json_config.c src/config.c src/kson.c src/vcmp.c $(HTSLIB)

bedadd: $(HTSLIB) 
	$(CC) -D_BED_ANNOS_MAIN $(DEBUG_CFLAGS) $(INCLUDES) -pthread -lz -o $@ src/anno_bed.c src/json_config.c src/config.c src/kson.c $(HTSLIB)

vcf2tsv: $(HTSLIB) version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ misc/vcf2tsv.c $(HTSLIB)

tsv2vcf: $(HTSLIB) version.h
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ misc/tsv2vcf.c src/table2hash.c $(HTSLIB)

vcf_rename_tags: $(HTSLIB) version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ misc/vcf_rename_tags.c $(HTSLIB)

bcfanno: $(HTSLIB) version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -lz -pthread -o $@ src/bcfanno_main.c src/vcf_annos.c src/anno_bed.c src/sequence.c src/genepred.c src/hgvs.c src/hgvs_vcf.c src/number.c src/vcmp.c src/json_config.c src/config.c src/kson.c src/file.c src/sort_list.c src/flank_seq.c $(HTSLIB)

bcfanno_debug: $(HTSLIB) version.h
	$(CC) -DDEBUG_MODE $(DEBUG_CFLAGS) $(INCLUDES) -lz -pthread -o $@ src/bcfanno_main.c src/vcf_annos.c src/anno_bed.c src/sequence.c src/genepred.c src/hgvs.c src/hgvs_vcf.c src/number.c src/vcmp.c src/json_config.c src/config.c src/kson.c src/file.c src/sort_list.c src/flank_seq.c $(HTSLIB)

test: $(HTSLIB) version.h vcfadd bedadd

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
