AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir) $(JELLYFISH_CFLAGS)
AM_CXXFLAGS = -Wall -Werror -g -O2 -std=c++0x
AM_LDFLAGS = $(JELLYFISH_LIBS)

# Yaggo automatic rules with silencing
V_YAGGO = $(V_YAGGO_$(V))
V_YAGGO_ = $(V_YAGGO_$(AM_DEFAULT_VERBOSITY))
V_YAGGO_0 = @echo "  YAGGO " $@;
.yaggo.hpp:
	$(V_YAGGO)$(YAGGO) --license $(srcdir)/header-license -p src/ -o $@ $<

JFLIB_I = include/jflib

bin_PROGRAMS = quorum_error_correct_reads quorum_combine_jf_dbs

quorum_error_correct_reads_SOURCES = src/error_correct_reads.cc src/kmer.cc	\
                                     src/err_log.cc src/yaggo.cpp
quorum_error_correct_reads_SOURCES += src/kmer.hpp				\
                                      src/error_correct_reads.hpp		\
                                      src/err_log.hpp				\
                                      src/error_correct_reads_cmdline.hpp	\
                                      src/yaggo.hpp				\
                                      include/gzip_stream.hpp			\
                                      include/misc.hpp				\
                                      $(JFLIB_I)/atomic_field.hpp		\
                                      $(JFLIB_I)/circular_buffer.hpp		\
                                      $(JFLIB_I)/compare_and_swap.hpp		\
                                      $(JFLIB_I)/locks_pthread.hpp		\
                                      $(JFLIB_I)/multiplexed_io.hpp		\
                                      $(JFLIB_I)/pool.hpp

quorum_combine_jf_dbs_SOURCES = src/combine_jf_dbs.cc src/yaggo.cpp
quorum_combine_jf_dbs_SOURCES += src/combine_jf_dbs.hpp src/yaggo.hpp

bin_SCRIPTS = quorum
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = src/quorum.in
do_subst = sed -e 's|[@]JELLYFISH[@]|$(JELLYFISH)|g' \
               -e 's|[@]bindir[@]|$(bindir)|g' \
	       -e 's|[@]JF_LIB_PATH[@]|$(JF_LIB_PATH)|g' \
	       -e 's|[@]PACKAGE_VERSION[@]|$(PACKAGE_VERSION)|g'

quorum: src/quorum.in Makefile
	$(do_subst) < $< > $@
	chmod +x $@
