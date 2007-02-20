##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2006 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
## License, or any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You find a copy of the GNU General Public License in the file
## license.txt along with this program; if not, write to the Free
## Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
##
##

## last $Author$
## last $Date$
## $Revision$

import unittest as U
import glob
import os.path

from Biskit.LogFile import StdLog
import Biskit.tools as T

## system-wide log file
TESTLOG = StdLog()

## categories
NORMAL = 0
LONG   = 1
PVM    = 2
EXE    = 3

CORE   = 10
DOCK   = 11
MOD    = 12


class BiskitTest( U.TestCase):
    """
    Base class for Biskit test cases.
    BiskitTest adds some functionality over the standard L{unittest.TestCase}:

    * self.local reflects whether the Test is performed in the __main__ scope
      of the module it belongs to -- rather than as part of the whole
      test suite. 

    * Each test case should be classified by assigning flags to
      its static GROUPS field -- this will be used by the test runner
      to filter out, e.g. very long tests or tests that require PVM.

    * self.log should capture all the output, by default to the TESTLOG
      instance defined at the module level.

    Usage:
       Biskit test cases should be created by sub-classing BiskitTest and
       adding one or more test_* methods. setUp and tearDown can optionally
       be overriden but should be called in the overriding method.
    """
    ## categories for which this test case qualifies
    GROUPS = [ NORMAL, CORE ] 

    def setUp( self ):
        self.local =  self.__module__ == '__main__'
        self.log = TESTLOG

    def tearDown( self ):
        if self.local:
            globals().update( self.__dict__ )



class FilteredTestSuite( U.TestSuite ):
    """
    Collection of BiskitTests filtered by categories.

    FilteredTestSuite silently ignores Test cases that are either

    * classified into any of the forbidden groups

    * not classified into any of the allowed groups

    The groups are defined in the L{Biskit.BiskitTest} module.

    By default (if initialized without parameters), all groups are allowed
    and no groups are forbidden.
    """

    def __init__( self, tests=(), allowed=[], forbidden=[] ):
        """
        @param tests: iterable of TestCases
        @type  tests: ( BiskitTest, )
        @param allowed: list of allowed groups
        @type  allowed: [ int ]
        @param forbidden : list of forbidden groups
        @type  forbidden : [ int ]
        """
        self._allowed   = allowed
        self._forbidden = forbidden
        super(self.__class__).__init__( self, tests=tests )


    def addTest( self, test ):
        """
        Add Biskit test case if it is matching the allowed and disallowed
        groups.
        @param test: test case
        @type  test: BiskitTest
        """
        assert isinstance( test, BiskitTest ), \
               'FilteredTestSuite only accepts BiskitTest instances'

        matches = [ g for g in test.GROUPS if g in self._forbidden ]
        if len( matches ) > 0:
            return

        matches = [ g for g in test.GROUPS if g in self._allowed ]

        if not self._allowed or (self._allowed and len(matches) > 0):
            super(self.__class__).addTest( test )
            

class BiskitTestLoader( object ):
    """
    A replacement for the unittest TestLoaders. It automatically
    collects all BiskitTests from a whole package (that means a
    folder with python files) rather than a single module.
    """

    def collectModules( self, path=T.projectRoot(), module='Biskit' ):
        """
        Import all python files of a package as modules. Sub-packages
        are ignored and have to be collected separately.
        @param path:  single search path for a package
        @type  path:  str
        @param module: name of the python package
        @type  module: str
        @raise ImportError, if a python file cannot be imported
        """
        module_folder = module.replace('.', os.path.sep)
        
        files = glob.glob( os.path.join( path, module_folder,'*.py' ) )

        files = map( T.stripFilename, files )
        files = [ f for f in files if f[0] != '_' ]

        r = []

        for f in files:
            try:
                r += [ __import__( '.'.join([module, f]), fromlist=[module]) ]
            except:
                pass  ## temporary // remove after testing

        return r

    def collectTests( self, modules ):
        pass

