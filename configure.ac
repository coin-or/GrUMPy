#============================================================================#
# This file is part of the GrUMPy graph library                               #
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

#############################################################################
#                       Names and other basic things                        #
#############################################################################
# autoconf version used by COIN
AC_PREREQ(2.59)

# name of project, version number, contact info
AC_INIT([GrUMPy],[0.6.1],[coin-bazaar@lists.coin-or.org])

# copyright for configure script
AC_COPYRIGHT([
#============================================================================#
# This file is part of the GrUMPy graph library                               #
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
])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/grumpy/BB.py)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################
# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT(GrUMPy)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

# Get the EXEXT variable for CYGWIN
AC_EXEEXT

#############################################################################
#                       Check for Python                                    #
#############################################################################

AM_PATH_PYTHON([2.6],,[:])
if test "$PYTHON" != "" ; then
  PY_PREFIX=`$PYTHON -c 'import sys ; print sys.prefix'`
fi
AM_CONDITIONAL([HAVE_PYTHON], [test "$PYTHON" != "" ])

#############################################################################
#                       Generate the Makefiles                              #
#############################################################################

AC_CONFIG_FILES([Makefile
		src/grumpy/Makefile
		test/Makefile])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
