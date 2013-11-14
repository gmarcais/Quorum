ACLOCAL_AMFLAGS=-I m4

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir) $(JELLYFISH_CFLAGS)
AM_CXXFLAGS = -Wall -g -O2 -std=c++0x
AM_LDFLAGS = $(STATIC_FLAGS)
LDADD = $(JELLYFISH_LIBS)

# Yaggo automatic rules with silencing
V_YAGGO = $(V_YAGGO_$(V))
V_YAGGO_ = $(V_YAGGO_$(AM_DEFAULT_VERBOSITY))
V_YAGGO_0 = @echo "  YAGGO " $@;
.yaggo.hpp:
	@echo yaggo $< $@
	$(V_YAGGO)$(YAGGO) --license $(srcdir)/header-license -o $@ $<
YAGGO_SOURCES = src/error_correct_reads_cmdline.hpp src/create_database_cmdline.hpp

BUILT_SOURCES = $(YAGGO_SOURCES)
noinst_HEADERS = $(YAGGO_SOURCES)

JFLIB_I = include/jflib

bin_PROGRAMS = quorum_error_correct_reads quorum_create_database
#bin_PROGRAMS = quorum_create_database

quorum_error_correct_reads_SOURCES = src/error_correct_reads.cc	\
                                     src/err_log.cc

# quorum_error_correct_reads_SOURCES += src/kmer.hpp			\
#                                       src/error_correct_reads.hpp	\
#                                       src/err_log.hpp			\
#                                       include/gzip_stream.hpp		\
#                                       include/misc.hpp			\
#                                       $(JFLIB_I)/atomic_field.hpp	\
#                                       $(JFLIB_I)/circular_buffer.hpp	\
#                                       $(JFLIB_I)/compare_and_swap.hpp	\
#                                       $(JFLIB_I)/locks_pthread.hpp	\
#                                       $(JFLIB_I)/multiplexed_io.hpp	\
#                                       $(JFLIB_I)/pool.hpp		\
#                                       include/divisor64.hpp

quorum_create_database_SOURCES = src/create_database.cc

###############################
# Adapter database generation #
###############################
dist_data_DATA = data/adapter_0
DISTCLEANFILES = data/adapter_0
SUFFIXES = _0 .fa

.fa_0:
	$(AM_V_GEN)mkdir -p $(@D); $(JELLYFISH) count -m 24 -s 5k -C -t 10 -r -o $$(dirname $@)/$$(basename $@ _0) $<


bin_SCRIPTS = quorum
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = src/quorum.in

if RELATIVE_PATHS
quorum_bindir=$$dir
else
quorum_bindir=$(bindir)
endif

if HAVE_RELATIVE_JF_PATH
quorum_jellyfish=$$dir/$(RELATIVE_JF_PATH)
else
quorum_jellyfish=$(JELLYFISH)
endif

do_subst = sed -e 's|[@]JELLYFISH[@]|$(quorum_jellyfish)|g' \
               -e 's|[@]bindir[@]|$(quorum_bindir)|g' \
	       -e 's|[@]JF_LIB_PATH[@]|$(JF_LIB_PATH)|g' \
	       -e 's|[@]PACKAGE_VERSION[@]|$(PACKAGE_VERSION)|g' \
	       -e 's|[@]JELLYFISH_VERSION[@]|$(JELLYFISH_VERSION)|g'

quorum: src/quorum.in Makefile
	$(do_subst) < $< > $@
	chmod +x $@

-include $(top_srcdir)/local.mk
