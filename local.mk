# Build the easy install script
bin_SCRIPTS += quorum_easy_install
quorum_easy_install: quorum_easy_install.in Makefile
	$(do_subst) < $< > $@
	chmod +x $@
distcheck: quorum_easy_install
dist: quorum_easy_install

AM_CXXFLAGS += -Werror

# Count lines of code
.PHONY: cloc
cloc:
	cloc --force-lang="Ruby,yaggo" --force-lang="make,am" --force-lang="make,mk" \
	  --exclude-dir="gtest" --ignored=cloc_ignored_src_files \
	  $(top_srcdir)/src $(top_srcdir)/include \
	  $(top_srcdir)/Makefile.am $(top_srcdir)/*.mk
