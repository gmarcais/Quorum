# Build the easy install script
bin_SCRIPTS += quorum_easy_install
quorum_easy_install: quorum_easy_install.in Makefile
	$(do_subst) < $< > $@
	chmod +x $@
distcheck: quorum_easy_install
dist: quorum_easy_install
