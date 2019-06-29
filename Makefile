PROG=       bcfanno vcf2tsv tsv2vcf vcf_rename_tags #GenePredExtGen
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
LIBS = -lz

all:$(PROG)

# See htslib/Makefile
PACKAGE_VERSION := $(shell git describe --tags)

bcfanno_version.h:
	echo '#define BCFANNO_VERSION "$(PACKAGE_VERSION)"' > $@


.SUFFIXES:.c .o
.PHONY:all clean clean-all clean-plugins distclean install lib tags test testclean force plugins docs

force:

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@


LIB_OBJ = src2/bed_utils.o src2/anno_bed.o src2/anno_col.o src2/anno_pool.o src2/sequence.o src2/gea.o src2/config.o src2/flank_seq.o src2/json_config.o src2/kson.o src2/name_list.o src2/number.o src2/sort_list.o src2/variant_type.o src2/vcf_annos.o src2/vcmp.o src2/anno_seqon.o src2/anno_vcf.o src2/anno_thread_pool.o

src2/bed_utils.o: src2/bed_utils.c
#src2/table2hash.o: src2/table2hash.c
src2/anno_bed.o: src2/anno_bed.c
src2/anno_col.o: src2/anno_col.c
src2/anno_pool.o: src2/anno_pool.c
src2/sequence.o: src2/sequence.c
src2/gea.o: src2/gea.c
src2/config.o: src2/config.c
src2/flank_seq.o: src2/flank_seq.c
src2/json_config.o: src2/json_config.c
src2/kson.o: src2/kson.c
src2/name_list.o: src2/name_list.c
src2/number.o: src2/number.c
src2/sort_list.o: src2/sort_list.c
src2/variant_type.o: src2/variant_type.c
src2/vcf_annos.o: src2/vcf_annos.c
src2/vcmp.o: src2/vcmp.c
src2/anno_seqon.o: src2/anno_seqon.c
src2/anno_vcf.o: src2/anno_vcf.c
src2/anno_thread_pool.o: src2/anno_thread_pool.c


liba.a: $(LIB_OBJ)
	@-rm -f src2/$@
	$(AR) -rcs src2/$@ $(LIB_OBJ)

#GenePredExtGen: $(HTSLIB)
#	$(CC) $(INCLUDES) -pthread -o $@ misc/genepred_ext_gen.c misc/ksw.c src2/anno_thread_pool.c src2/genepred.c src2/number.c src2/faidx_def.c $(HTSLIB) $(LIBS)

vcf2tsv: $(HTSLIB) bcfanno_version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -o $@ misc/vcf2tsv.c $(HTSLIB) $(LIBS)

tsv2vcf: $(HTSLIB) bcfanno_version.h
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -o $@ misc/tsv2vcf.c misc/table2hash.c $(HTSLIB) $(LIBS)

vcf_rename_tags: $(HTSLIB) bcfanno_version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -o $@ misc/vcf_rename_tags.c $(HTSLIB) $(LIBS)

#hgvs: $(HTSLIB) bcfanno_version.h
#	$(CC) $(CFLAGS) $(INCLUDES) -pthread -o bcfanno_hgvs -DANNO_HGVS_MAIN  src2/anno_col.c src2/anno_hgvs.c src2/hgvs.c src2/name_list.c src2/anno_thread_pool.c src2/anno_pool.c src2/number.c src2/vcmp.c src2/genepred.c src2/sort_list.c src2/variant_type.c $(HTSLIB) $(LIBS)

#bcfanno_atac: $(HTSLIB) bcfanno_version.h
#	$(CC) $(CFLAGS) $(INCLUDES) -pthread -o $@ src2/bed_utils.c src2/atac.c src2/number.c  $(HTSLIB) $(LIBS)

#bcfanno_pwm: $(HTSLIB) bcfanno_version.h
#	$(CC) $(DEBUG_CFLAGS) $(INCLUDES) -pthread -o $@ src2/bed_utils.c src2/motif.c src2/number.c src2/wrap_pileup.c src2/anno_col.c src2/anno_thread_pool.c src2/anno_pool.c src2/sequence.c $(HTSLIB) $(LIBS)

bcfanno: $(HTSLIB) liba.a bcfanno_version.h 
	$(CC) $(CFLAGS) $(INCLUDES) -pthread -o $@ src2/bcfanno_main.c src2/liba.a $(HTSLIB) $(LIBS)


bcfanno_debug: $(HTSLIB) bcfanno_version.h
	$(CC) -DDEBUG_MODE $(DEBUG_CFLAGS) $(INCLUDES) -pthread -o $@ src2/bcfanno_main.c src2/liba.a $(HTSLIB) $(LIBS)

test: $(HTSLIB) bcfanno_version.h

clean: testclean
	-rm -f gmon.out *.o *~ $(PROG) bcfanno_version.h 
	-rm -rf *.dSYM plugins/*.dSYM test/*.dSYM
	-rm -f anno_vcf bedadd vcfadd bcfanno anno_bed hgvs_generate hgvs_vcf GenePredExtGen bcfanno_hgvs
	-rm -f config bcfanno_debug vcf2tsv tsv2vcf vcf_rename_tags

testclean:
	-rm -f test/*.o test/*~ $(TEST_PROG)

distclean: clean
	-rm -f TAGS

clean-all: clean clean-htslib

tags:
	ctags -f TAGS *.[ch] plugins/*.[ch]
