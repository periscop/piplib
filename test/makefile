#
# Copyright 2002 (c) Leiden University, LIACS,
# All Rights Reserved.
#
# Permission is granted to copy, use, and distribute
# for any commercial or noncommercial purpose under the terms
# of the GNU General Public license, version 2, June 1991
# (see file : LICENSING).
#
# @author: Bart Kienhuis
# $Id: makefile,v 1.26 2002/09/26 13:03:53 kienhuis Exp $
#

ROOT = ..

# pbs with boulet bouleti and kirk 
PIPTEST =  \
	boulet.dat \
	crescat.dat \
	discr.dat \
	invert.dat \
	lineri.dat \
	loz.dat \
	maxb.dat \
	max.dat \
	test12i.dat

BITS = 64

test:
	failedtest=0; \
	for x in $(PIPTEST) ; do \
		echo "Verify file $$x" ; \
		$(ROOT)/obj$(BITS)_$(TARGET)/pip$(BITS)$(EXEC_SUFFIX) < $$x > xyz ;\
		diff -w xyz `basename $$x .dat`.ll ; \
		result=$$?; \
		if [ "$$result" -eq "1" ]; then \
			echo "Error: $$x is not the same"; \
			failedtest=`expr $$failedtest + 1`; \
		else \
			echo "$$x passed"; \
	        fi; \
	done ; \
	if [ $$failedtest != 0 ]; then \
		echo "$$failedtest tests failed"; \
	else \
	        echo "Test pip works correctly"; \
        fi

# Include the shared compilation rules
include $(ROOT)/Makefile
