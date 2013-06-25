#============================================================================#
# This file is part of the GrUMPy library                                    #
#                                                                            #
# GrUMPy is distributed under the Eclipse Public License as part of the      #
# COIN-OR repository (http://www.coin-or.org).                               #
#                                                                            #
# Authors: Brady Hunsaker, Google (bhunsaker@google.com)                     #
#          Osman Ozaltin, University of Waterloo (oozaltin@waterloo.edu)     #
#          Aykut Bulut (ayb211@lehigh.edu), Lehigh Univesity                 #
#          Ted Ralphs (ted@lehigh.edu), Lehigh University                    #
#                                                                            #
#                                                                            #
# Copyright (C) 2009-2013, Google, Lehigh University, University of Waterloo,#
# Brady Hunsaker, Osman Ozaltin, Aykut Bulut, and Ted Ralphs                 #
# All Rights Reserved.                                                       #
#============================================================================#

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS  = src/grumpy

# We don't want to compile the test subdirectory, unless the test target is
# specified.  But we need to list it as subdirectory to make sure that it is
# included in the tarball

if ALWAYS_FALSE
  SUBDIRS += test
endif

########################################################################
#             Additional files to be included in tarball               #
########################################################################

# Here we need include all files that are not mentioned in other Makefiles
EXTRA_DIST = 

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test; $(MAKE) test

unitTest: test

clean-local:
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

.PHONY: test unitTest 
