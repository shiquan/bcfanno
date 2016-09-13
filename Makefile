PROG=       vcfanno

all: $(PROG) $(TEST_PROG)

# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = htslib-1.3
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a

CC       = gcc
CFLAGS   = -g -Wall -Wc++-compat -O0
DFLAGS   =
#OBJS     = main.o vcfindex.o tabix.o \
#           vcfstats.o vcfisec.o vcfmerge.o vcfquery.o vcffilter.o filter.o vcfsom.o \
#           vcfnorm.o vcfgtcheck.o vcfview.o vcfannotate.o vcfroh.o vcfconcat.o \
#           vcfcall.o mcall.o vcmp.o gvcf.o reheader.o convert.o vcfconvert.o tsv2vcf.o \
#           vcfcnv.o HMM.o vcfplugin.o consensus.o ploidy.o version.o \
#           ccall.o em.o prob1.o kmin.o # the original samtools calling
INCLUDES = -I. -I$(HTSDIR)/


all:$(PROG)

# See htslib/Makefile
PACKAGE_VERSION = 0.01
ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)
DOC_VERSION :=  $(shell git describe --always)+
DOC_DATE := $(shell date +'%Y-%m-%d %R %Z')
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif
version.h:
	echo '#define VCFANNO_VERSION "$(PACKAGE_VERSION)"' > $@


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

hgvs_generate: $(HTSLIB) 
	$(CC) -D_HGVS_MAIN $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ hgvs_generate.c $(HTSLIB)

vcfadd: $(HTSLIB) 
	$(CC) -D_VCF_ANNOS_MAIN $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ vcf_annos.c config.c kson.c vcmp.c $(HTSLIB)

bedadd: $(HTSLIB) 
	$(CC) -D_BED_ANNOS_MAIN $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ anno_bed.c config.c kson.c $(HTSLIB)

vcf2tsv: $(HTSLIB) version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ misc/vcf2tsv.c $(HTSLIB)

vcf_rename_tags: $(HTSLIB) version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ misc/vcf_rename_tags.c $(HTSLIB)

vcfanno: $(HTSLIB) version.h vcf2tsv vcf_rename_tags
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ anno_core.c vcmp.c config.c kson.c vcf_annos.c anno_bed.c hgvs_generate.c $(HTSLIB)

vcfanno_debug: $(HTSLIB) version.h vcf2tsv vcf_rename_tags
	$(CC) -DDEBUG_MODE $(CFLAGS) $(INCLUDES) -pthread -lz -o $@ anno_core.c vcmp.c config.c kson.c vcf_annos.c anno_bed.c hgvs_generate.c $(HTSLIB)

test: $(HTSLIB) version.h hgvs_generate vcfadd bedadd

clean: testclean
	-rm -f gmon.out *.o *~ $(PROG) version.h 
	-rm -rf *.dSYM plugins/*.dSYM test/*.dSYM
	-rm -f anno_vcf bedadd vcfadd vcfanno anno_bed hgvs_generate

testclean:
	-rm -f test/*.o test/*~ $(TEST_PROG)

distclean: clean
	-rm -f TAGS

clean-all: clean clean-htslib

tags:
	ctags -f TAGS *.[ch] plugins/*.[ch]
