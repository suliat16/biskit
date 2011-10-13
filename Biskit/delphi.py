##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2011 Raik Gruenberg
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
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

## last $Author: graik $
## last $Date: 2009-05-09 14:17:28 +0200 (Sat, 09 May 2009) $
## $Revision: $
"""
Wrapper for Delphi -- Poisson-Boltzman electrostatic potential calculation.
"""

import tempfile, os
import numpy as N
import re

from Biskit import Executor, PDBModel, Reduce
from Biskit import AmberPrepParser
import Biskit.tools as T
import Biskit.mathUtils as U

class DelphiError( Exception ):
    pass

class DelphiCharges( object ):
    """
    Helper class for Delphi to write out a (custom) Delphi charge file. Atomic
    partial charges are taken from a dictionary of AmberResidue instances. 
    This dictionary can be created with the AmberResidue class from standard
    Amber topology files (bundled with Biskit in Biskit/data/amber/residues).
    
    Residues in the input PDB are matched against Amber residues with exactly 
    the same atom content (same set of atom names). That means, the name of 
    the residue can differ (as is the case for e.g. CALA, CTHR or NALA, etc.).
    
    DelphiCharges.customCharges() can be applied iteratively with different 
    sets of residue descriptions until all residues in the input model have been 
    matched against a topology entry.
    """
    
    def __init__(self, restypes={} ):
        """
        @param fname: output file name
        @type  fname: str
        @param restypes: dict of Residue type definitions indexed by 3-letter 
                         code
        @type  restypes: { str: AmberResidueType }
        """
        self.restypes = restypes
    
    def res2delphi(self, restype, resnumber=None, comment='' ):
        """
        Generate the Delphi charge record for a single residue (residue type).
        @type restype: AmberResidueType
        @param resnumber: residue number to use in the Delphi record
        @type  resnumber: int
        @return: several lines formatted as in Delphi charge files
        @rtype: str
        """
        resnumber = resnumber or ''

        r = ''
        for atom in restype.atoms.iterDicts():
 
            atom['resnumber'] = str(resnumber)
            atom['rescode'] = restype.code
            r += '%(name)-5s %(rescode)3s %(resnumber)-5s %(charge)6.3f'%\
                 atom

            if comment:
                r += '  ! %s' % comment
                comment = ''  # only print once
            r += '\n'
        return r
    
    def tofile( self, fname ):
        fname = T.absfile( fname )
        f = open( fname, 'w' )
        f.write('! charge file generated by Biskit.delphi.DelphiCharges\n')
        f.write('atom__resnumbc_charge_\n')
        try:
            for res in self.restypes:
                f.write( self.res2delphi( res ) )
        finally:
            f.close()
        
    
    def checkmodel( self, m ):
        """
        Verify that all residues and atoms in model m are covered by charges.
        @type m: PDBModel
        @return: residues not described or with missing atoms, atom names
        @rtype: ([int], [str])
        """
        assert( isinstance(m,PDBModel) )
        missing = {}

        for i, resm in enumerate( m.resModels() ):
            
            resname = resm['residue_name'][0]
            missing_atm = []
            missing_res = not resname in self.restypes
            
            if not missing_res:
                resatoms = resm.atomNames()
                standard = self.restypes[resname].atomNames()
                missing_atm = U.difference( resatoms, standard )
            
            if missing_res or missing_atm:
                missing[i] = missing_atm
        
        return missing
    
    def atomkey( self, residue ):
        """
        Create a string key encoding the atom content of residue.
        @param residue: model or AmberResidue
        @type  residue: PDBModel or AmberResidue
        @return: key formed from alphabetically sorted atom content of residue
        @rtype: str
        """
        r = residue.atomNames()
        r.sort()
        r = ''.join( r )
        return r
        
    
    def resindex( self, *resdics):
        """
        Build an index of residue types indexed by ordered atom content.
        @param *resdics: additional restype dictionaries indexed by res. name
        @type  *resdics: {str_resname: AmberResidue}
        @return: dict with alphabetically ordered atom content as key
        @rtype : { str_atoms : AmberResidue }
        """
        r = {}
        ## default residues treated last to give them priority
        resdics = resdics + ( self.restypes, )  

        for rdic in resdics:
            for resname, restype in rdic.items():
                r[ self.atomkey(restype) ] = restype

        return r
    
    def customCharges( self, model, comment='' ):
        """
        @param model: structure for which Delphi charge file should be created
        @type  model: PDBModel
        @return: Delphi charge table, list of non-matching residues
        @rtype: str
        """
        index = self.resindex()
        missing = []
        r = ''
        
        for res in model.resModels():
            rtype = index.get( self.atomkey( res ), None)
            resnumber = res['residue_number'][0]

            if rtype is not None:
                
                rtype.code = res['residue_name'][0]
                    
                r += self.res2delphi( rtype, resnumber=resnumber, 
                                      comment=comment )
                comment = ''  ## erase comment
            else:
                missing += [ res ]
    
        return r, missing
        
   
    
class Delphi( Executor ):
    """
    Calculate electrostatic potentials and potential maps with Delphi.
    [Work in Progress]
    
    The current workflow of this wrapper is:
    1. take input model/structure, remove hydrogens
    2. add and optimize hydrogens with the reduce program
    3. adapt residue and atom names to Amber conventions
    4. match each residue (by atom content) to a residue from a list of 
       Amber residue topology files
    5. Create custom delphi charge file with Amber partial charges
    6. Run Delphi in temporary folder 
    7. parse result energies into result dictionary
    
    
    Usage
    =====

    >>> D = Delphi( inputmodel )
    >>> result = D.run()
    >>> result
    {'scharge' :  1.4266   # surface charge
     'egrid' :  9105.51    # total grid energy
     'ecoul' :  -9849.664  # couloumb energy
     'erxn'  :  -664.7469   # corrected reaction field energy
     'erxnt' :  -21048.13  # total reaction field energy
     'eself' :  -20383.39 }  # self reaction field energy
    
    Note
    ====
    
    All energy values are in units of kT. The important terms for the
    calculation of (free) energy differences are 'egrid', 'erxn' and 'ecoul'.
    
    Grid dimensions
    ===============
    
    If no further options are given, the Delphi wrapper will determine grid
    dimensions that are centered on the molecule and lead to at most 60%
    "filling" of the grid in any of the x, y or z dimensions. In other words,
    the longest dimension of the solute along x, y or z axis will still have a
    20% margin to the grid boundary. This fill-factor can be overriden with
    the 'perfil' parameter.
    
    The density of the grid is, by default, 1.2 points per Angstroem
    and can be adjusted with the 'scale' parameter. The number of grid points 
    (gsize) will be calculated accordingly. 
    
    For a given combination of scale and perfil, this calculation will give
    approximately (but not exactly) the same grid dimensions as the one the
    delphi program would determine by itself. Minor differences may lead to
    two grid points more or less. This can also lead to some changes in
    resulting energies.
    
    For more control over grid dimensions, you should call the method
    Delphi.setGrid() *before* Delphi.run(). For example:
    
    >>> D = Delphi( inputmodel )
    >>> D.setGrid( scale=1.2, perfil=80 )
    
    ... will fix a grid centered on the molecule, with 1/1.2 A grid spacing and 
    at least 10% distance margin to the boundary. Another example:
    
    >>> D = Delphi( inputmodel )
    >>> D.setGrid( acenter=(0,0,0), scale=1.0, gsize=100 )
    
    ... circumvents any calculation and fixes a 100 x 100 x 100 grid centered 
    on 0 and with 1/1.0 A grid spacing -- a box of 100 A x 100 A x 100 A. If
    the scale parameter is not given, it will default to whatever was specified
    at the Delphi() constructor (and from there default to 1.2). 
    
    Delphi.setGrid() returns the three grid parameters as a dictionary for 
    later re-use. So if you want to use the same grid dimensions for two 
    different structures you can first calculate the dimensions based on one
    structure and then apply the same dimensions to another Delphi run:
    
    >>> D1 = Delphi( complex )
    >>> grid1 = D1.setGrid( scale=2.4, perfil=60 )
    >>> print grid1
    {'acenter': (0.1, 10., -0.5), 'scale': 2.400, 'gsize': 99 }
    
    >>> ligand = complex.takeChains( [1] ) ## extract part of structure
    >>> D2 = Delphi( ligand )
    >>> D2.setGrid( **grid1 )
    
    Note: The '**' syntax unpacks the dictionary as keyword parameters into the 
    method.
    
    
    Customization
    =============
    
    The default delphi parameter file can be replaced (parameter template).
    Note though that you should keep the place holders for input and output
    files. The default parameter file is taken from:

        Biskit/data/delphi/delphi_simple.prm
        
    As always, the place holders (e.g. %(salt)f ) in this file are replaced
    by the value of a variable of the same name (e.g. salt) within the name
    space of the Executor instance. In other words:
    
    >>> D = Delphi( inputmodel, salt=0.2, ionrad=2.5 )
    
    ...will override the default values for salt and ionradius for which there
    are the place hoders %(salt)f and %(ionrad)f in the template file.
    
    
    The default handling and matching of atomic partial charges can be 
    modified by:
    
        * providing a ready-made Delphi charge file (parameter f_charges)
        * or providing an alternative list of Amber topology files from which 
          residues are looked up by their atom content (parameter topologies)

    @note: Command configuration: biskit/Biskit/data/defaults/exe_delphi.dat
    """

    F_RADII = 'default.siz'           ## default Delphi atom radius file
    F_PARAMS = 'delphi_simple.prm'    ## default delphi parameter file
    
    ## list of Amber topology files in decending priority
    ## will be mined for matching residues to assign charges
    F_RESTYPES = ['all_amino03.in',
                  'all_aminoct03.in',
                  'all_aminont03.in',
                  'all_nuc02.in' ]
    
    RE_E_GRID = r'total grid energy\s+:\s+(?P<egrid>[0-9\-\.]+)\s+kt'
    RE_E_COUL = r'coulombic energy\s+:\s+(?P<ecoul>[0-9\-\.]+)\s+kt'
    RE_E_SELF = r'self-reaction field energy\s+:\s+(?P<eself>[0-9\-\.]+)\s+kt'
    RE_E_RXN  = r'corrected reaction field energy\s*:\s+(?P<erxn>[0-9\-\.]+)\s+kt'
    RE_E_RXNT = r'total reaction field energy\s*:\s+(?P<erxnt>[0-9\-\.]+)\s+kt'
    RE_SURFCH = r'total s\.charge\,no epsin carrying\s+:\s+(?P<scharge>[0-9\-\.]+)'
    
    
    def __init__( self, model, template=None, topologies=None,
                  f_charges=None,
                  indi=4.0, exdi=80.0, salt=0.15, ionrad=2, prbrad=1.4, 
                  bndcon=4, scale=1.2, perfil=60, 
                  **kw ):
        """
        @param model: structure for which potential should be calculated
        @type  model: PDBModel
        @param template: delphi command file template [None=use default]
        @type  template: str
        @param f_radii: alternative delphi atom radii file [None=use default]
        @type  f_radii: str
        @param topologies: alternative list of residue charge/topology files
                           [default: amber/residues/all*]
        @type  topologies: [ str ]
        @param f_charges: alternative delphi charge file 
                          [default: create custom]
        @type  f_charges: str

        @param indi: interior dilectric (4.0)
        @param exdi: exterior dielectric (80.0)
        @param salt: salt conc. in M (0.15)
        @param ionrad: ion radius (2)
        @param prbrad: probe radius (1.4) 
        @param bndcon: boundary condition (4, delphi default is 2)
        @param scale:  grid spacing (1.2)
        @param perfil: grid fill factor in % (for automatic grid, 60) 
        
        @param kw: additional key=value parameters for Executor:
        @type  kw: key=value pairs
        ::
          debug    - 0|1, keep all temporary files (default: 0)
          verbose  - 0|1, print progress messages to log (log != STDOUT)
          node     - str, host for calculation (None->local) NOT TESTED
                          (default: None)
          nice     - int, nice level (default: 0)
          log      - Biskit.LogFile, program log (None->STOUT) (default: None)
        """
        template = template or T.dataRoot() + '/delphi/' + self.F_PARAMS
        
        tempdir = self.newtempfolder( tempdir=True )  ## create new temp folder
        f_in = tempfile.mktemp( '.inp', 'delphi_', dir=tempdir )
        
        self.f_pdb = tempfile.mktemp( '.pdb', 'delphi_', dir=tempdir)
        self.f_grid = tempfile.mktemp( '_grid.phi', 'delphi_', dir=tempdir )
        self.f_map = None
        self.f_radii = None
        self.topologies = topologies or self.F_RESTYPES
        self.f_charges = f_charges or tempfile.mktemp( '.crg', 'delphi_',
                                                       dir=tempdir )
        
        ## DELPHI run parameters
        self.indi=indi  # interior dilectric(4.0)
        self.exdi=exdi  # exterior dielectric(80.0)
        self.salt=salt  # salt conc. in M (0.15)
        self.ionrad=ionrad # ion radius (2)
        self.prbrad=prbrad # probe radius (1.4) 
        self.bndcon=bndcon # boundary condition (4, delphi default is 2)
        
        ## DELPHI parameters for custom grid
        self.scale=scale   # grid spacing (1.2)
        self.perfil=perfil # grid fill factor in % (for automatic grid, 60)
        self.gsize = None
        self.acenter = None
        self.strcenter = '(0.0,0.0,0.0)'
        
        Executor.__init__( self, 'delphi', 
                           template=template,
                           f_in=f_in,
                           args=f_in,
                           catch_err=True,
                           tempdir=tempdir,
                           cwd=tempdir,
                           **kw )
        
        self.model = model
        self.delphimodel = None
    

    def delphiDimensions( self, model ):
        """
        Calculate "geometric" center and molecular dimensions as defined by 
        Delphi (the delphi geometric center is NOT exactly what a geometric
        center is commonly defined as). 
        @param model: PDBModel for which center and dimensions should be
                      calculated
        @type model:  PDBModel
        @return: center and dimensions
        @rtype : N.array([x,y,z] of float), N.array([x,y,z] of float)
        """
        m = model.compress( model.maskHeavy() )
        xyz = m.getXyz()
        
        ## largest protein length in x, y, z direction
        dimensions = N.max( xyz, 0 ) - N.min( xyz, 0 )
        
        center = N.min( xyz, 0) + dimensions / 2
        
        ## + 1 C diameter
        dimensions += 2*1.7
        
        return center, dimensions
        

    def setGrid( self, acenter=None, gsize=None, scale=None, perfil=None ):
        """
        Specify or calculate Delphi grid. There are two options:
        (1) specify the actual grid by giving acenter, scale and gsize.
        (2) calculate new grid dimensions from acenter, scale and perfil.
        If not given, acenter defaults to the geometric center of the structure
        model.
        @param acenter: center coordinates for the grid
        @type  acenter: [float, float, float]
        @param gsize: number of grid points in x, y, and z direction
        @type  gsize: int
        @param scale: distance between grid points
        @type scale : float
        @param perfil: percent fill factor 
        """
        center, dimensions = self.delphiDimensions( self.model )

        if acenter is not None:
            self.acenter = acenter
        else:
            if self.acenter is None:
                self.acenter = center
        
        self.strcenter = str( tuple( self.acenter) )
        self.scale= scale or self.scale
        
        if not gsize:
            self.perfil = perfil or self.perfil
              
            ## grid size in number of points at self.scale per Angstrom density
            gsize = (N.max( dimensions ) * 100. / self.perfil) * self.scale
            gsize = int( round( gsize ) )
            
            ## grid size must be an uneven number
            if not gsize % 2:
                gsize += 1
                
            self.gsize = gsize
        
        return {'acenter':self.acenter, 'scale':self.scale, 'gsize':self.gsize}
        
        
    def version(self):
        return 'Delphi $Revision: $'


    def __prepareFolder( self ):
        """
        Link default parameter files into working directory.
        """
        try:
            f_radii = self.f_radii or T.dataRoot() + '/delphi/' + self.F_RADII

            target = os.path.join(self.cwd, 'radii.siz')
            if not os.path.exists( target ):
                os.symlink( f_radii, target )

        except OSError, error:
            raise DelphiError, \
                  'Error preparing temporary folder for Delphi\n'+\
                  'Error: %r\n' % error +\
                  'folder: %r\n' % self.cwd
        

    def __prepareCharges(self, f_out ):

        unhandled = self.delphimodel
        
        try:
            f = open( f_out, 'w' )
            f.write('! custom charge file generated by Biskit.delphi.Delphi\n')
            f.write('atom__resnumbc_charge_\n')
            
            for topo in self.topologies:
                if unhandled:
                    dc = DelphiCharges( AmberPrepParser(topo).residueDict() )

                    r, missing = dc.customCharges( unhandled, 
                                                   comment= 'from %s:'% topo)
                    f.write( r )
                    
                    if missing: ## create new model with only unhandled residues
                        unhandled = missing[0].concat( *missing[1:] )
                    else:
                        unhandled = None
            
            f.close()
          
        except IOError, why: 
            raise IOError, 'Error creating custom delphi charge file '+f_out+\
                  '( '+str(why)+' )'

    
    def prepare( self ):
        """
        Overrides Executor method.
        """
        Executor.prepare( self )
        
        self.__prepareFolder()
        
        ## if setGrid hasn't been called yet, create automatic grid
        if not self.gsize:
            self.setGrid()
        
        reducer = Reduce( self.model, verbose=self.verbose,
                          tempdir=self.tempdir, cwd=self.cwd,
                          log=self.log, debug=self.debug )
        if self.verbose: 
            self.log.add('adding hydrogen atoms to input structure')

        self.delphimodel = reducer.run()
        self.delphimodel.xplor2amber()
        self.delphimodel.writePdb( self.f_pdb )
        
        if not os.path.exists( self.f_charges ):
            self.__prepareCharges( self.f_charges )
        

    def cleanup( self ):
        """
        Tidy up the mess you created.
        """        
        if not self.debug:
            T.tryRemove( self.f_pdb )

        Executor.cleanup( self )

    def isFailed( self ):
        """
        Overrides Executor method
        """
        return self.output is None or \
               not 'energy calculations done' in self.output

    def fail( self ):
        """
        Overrides Executor method. Called when execution fails.
        """
        s = 'Delphi failed. Please check the program output in the '+\
          'field `output` of this Delphi instance (e.g. `print x.output`)!'
        self.log.add( s )

        raise DelphiError, s

    def postProcess( self ):
        """
        Called directly after execution. Read delphi output.
        """
        try:
            f = open( self.f_out, 'r')
            self.output = f.read()
            f.close()
        except IOError:
            self.output = None
    

    def parseOutput( self ):
        """
        Assumes output file has been parsed into self.output
        """
        r = {}
        for pattern in [self.RE_E_COUL, self.RE_E_GRID, self.RE_E_RXN, 
                        self.RE_E_SELF, self.RE_E_RXNT, self.RE_SURFCH]:
            ex = re.compile( pattern )
            hit = ex.search( self.output )
            try:
                r.update( hit.groupdict() )
            except:
                self.log.writeln('Warning, no match for: ' + pattern)
        
        for k, v in r.items():
            r[k] = float( v )
            
        return r
    
    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        self.result = self.parseOutput()
            


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test class"""

    TAGS = [ BT.EXE ]
    MODEL= None

    def test_delphi( self ):
        """Delphi test"""
        if self.local: print 'Loading PDB...'

        self.m1 = self.MODEL or \
            PDBModel( T.testRoot( 'lig/1A19_dry.model' ) )
        Test.MODEL = self.m1


        if self.local: print 'Starting Delphi'
        self.x = Delphi( self.m1, debug=self.DEBUG,
                         verbose=self.local )

        if self.local:
            print 'Running'

        self.r = self.x.run()

        if self.local:
            print "Result: "
            print self.r
            
        expect = {'scharge': 1.427, 'egrid': 9075., 'ecoul': -9849.70, 
                  'eself': -20383.4, 'erxn': -666.7}
        
        for k, v in expect.items():
            self.assertAlmostEqual( expect[k], self.r[k], 0, 
                                    'missmatch in energy values: '+k)
        
            
    def test_delphiCharges( self ):
        """DelphiCharges test"""
        resnormal = AmberPrepParser( 'all_amino03.in' ).residueDict()
        resnterm  = AmberPrepParser( 'all_aminont03.in').residueDict()
        rescterm  = AmberPrepParser( 'all_aminoct03.in').residueDict()
        
        if self.local:
            T.errWrite( 'loading PDB...' )

        self.m1 = self.MODEL or PDBModel( T.testRoot( 'lig/1A19_dry.model' ) )
        Test.MODEL = self.m1
        if self.local:
            T.errWriteln( 'Done.' )
        
        if self.local:
            T.errWrite( 'Adding hydrogens to model (reduce)...' )

        self.rmodel = Reduce( self.m1 ).run()
        self.rmodel.xplor2amber()
        if self.local:
            T.errWriteln( 'Done.' )
        
        self.dc = DelphiCharges( resnormal )
        if self.local:
            print
            print self.dc.res2delphi( self.dc.restypes['CYS'] )

        missing = self.dc.checkmodel( self.rmodel )
        if self.local:
            print 'atom types identified as missing: ', missing
        
        self.assert_( 0 in missing.keys(), 'residue missmatch (0)' )
        self.assert_( 88 in missing.keys(), 'residue missmatch (88)' )
        self.assert_( U.difference(['H1', 'H2', 'H3'], missing[0]) == [],
                      'atom missmatch (0)' )
        self.assert_( 'OXT' in missing[88], 'atom missmatch (88)' )
        
        self.customCharges = self.dc.customCharges( self.rmodel )
        
        

if __name__ == '__main__':

    BT.localTest(debug=True)
    