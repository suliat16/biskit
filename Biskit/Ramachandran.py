##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
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
## last $Date$
## $Revision$
## last $Author:  
"""
Display a Ramachandran plot for a list of PDBModels
"""

from Biskit import ColorSpectrum as CS
from Biskit.MatrixPlot import Legend
import biggles
from Biskit.PDBDope import PDBDope 

import Biskit.mathUtils as MU
import Biskit.PDBModel as PDBModel
import Biskit.tools as T

import Numeric as N


class Ramachandran:

    def __init__( self, models, name=None, profileName='relAS',
                  noPlot=0):
        """
        @param model: List of models display a Ramachandran plot for
        @type  model: [ PDBModel ] OR PDBModel
        @param modName: model name, will show up in plot
        @type  modName: str
        @param profileName: name of profile to use for coloring
                            (default: 'relAS')
        @type  profileName: str
        """
        if type(models) != type([]):
            models = [ models ]

        self.psi = []
        self.phi = []

        self.gly = []
        self.pro = []

        self.prof=[]
        self.profileName =  profileName

        self.name=name

        # calculate angles, profiles ...
        self.calc( models, noPlot )
        
        self.prof = N.ravel(self.prof)
        self.gly = N.ravel(self.gly)
        self.pro = N.ravel(self.pro)


    def calc( self, models ):
        """
        Calculate angles, profiles and other things needed.
        
        @param models: List of models 
        @type  model: [ PDBModel ]
        """              
        res_count = 0
        for m in models:
            
            ## add profile if not there
            if self.profileName:
                self.prof += [ self.calcProfiles( m ) ]
                     
            ## calclate phi and psi angles for model
            self.phi_and_psi( m )

            ## get list with GLY and PRO residue indices
            gly_atomInd = m.indices(lambda a: a['residue_name']=='GLY')
            gly_resInd  = N.array( m.atom2resIndices( gly_atomInd ) )
            pro_atomInd = m.indices(lambda a: a['residue_name']=='PRO')
            pro_resInd  = N.array( m.atom2resIndices( pro_atomInd ) )    
            self.gly.append( gly_resInd + res_count )
            self.pro.append( pro_resInd + res_count )
            res_count += m.lenResidues()
            

    def calcProfiles( self, m ):
        """
        Calculate needed profiles.

        @param m: PDBModel to calculate data for
        @type  m: PDBModel
        """
        print "Initiating PDBDope..."
        d = PDBDope( m )
                
        if not self.profileName in m.aProfiles.keys():
            
            if self.profileName in ['MS', 'AS', 'curvature', 'relAS', 'relMS']:
                print "Adding SurfaceRacer profile...",
                d.addSurfaceRacer()
                            
            if self.profileName in ['relASA']:
                print "Adding WhatIf ASA...",
                d.addASA()
                        
            if self.profileName in ['density']:
                print "Adding surface density...",
                d.addDensity()
                                     
        if not self.profileName in m.rProfiles.keys():
                    
            if self.profileName in ['cons_abs', 'cons_max', 'cons_ent']:
                print "Adding conservation data...",
                d.addConservation()
                                               
            if self.profileName in ['ASA_total','ASA_sc', 'ASA_bb']:
                print "Adding WhatIf ASA...",
                d.addASA()
      
        print 'Done.'

        ## convert atom profiles to average residue profile
        if self.profileName in m.aProfiles.keys():
            prof = []
            aProfile = m.profile( self.profileName )
            resIdx =  m.resIndex()
            resIdx += [ m.lenAtoms()]
            for i in range(len(resIdx)-1):
                prof+=[ N.average( N.take(aProfile, range(resIdx[i],
                                                          resIdx[i+1]) ) )]
        else:
            prof = m.profile( self.profileName )

        return prof

    
    def phi_and_psi( self, model ):
        """
        Calculate phi and psi torsion angles for all
        residues in model::
        
          phi - rotation about the N-CA bond
              - last position in a chain = None
          psi - rotation about CA-C
              - first position in a chain = None          

        @param model: PDBModel
        @type  model: PDBModel 
        """
        for c in range( model.lenChains(breaks=1) ):
            cModel = model.takeChains( [c], breaks=1 )

            xyz = cModel.xyz

            xyz_CA =  N.compress( cModel.maskCA(), xyz,0 )
            xyz_N  =  N.compress( cModel.mask( ['N'] ), xyz,0 )
            xyz_C  =  N.compress( cModel.mask( ['C'] ), xyz,0 )

            ## phi: c1 - N
            ##      c2 - CA
            ##      c3 - C
            ##      c4 - N of next residue
            for i in range( len(xyz_N)-1 ):
                self.phi += [self.dihedral( xyz_N[i], xyz_CA[i],
                                            xyz_C[i], xyz_N[i+1] )]
            self.phi += [None]

            ## psi: c1 - C of previous residue  
            ##      c2 - N
            ##      c3 - CA
            ##      c4 - C
            self.psi += [None]
            for i in range( 1, len(xyz_N) ):
                self.psi += [self.dihedral( xyz_C[i-1], xyz_N[i],
                                            xyz_CA[i], xyz_C[i] )]

    def dihedral( self, c1, c2, c3, c4 ):
        """
        Calculates the torsion angle of a set of four atom coordinates.
        The dihedral angle returned is the angle between the projection
        of i1-i2 and the projection of i4-i3 onto a plane normal to i2-i3.

        @param c1: coordinates
        @type  c1: [float]
        @param c2: coordinates
        @type  c2: [float]
        @param c3: coordinates
        @type  c3: [float]
        @param c4: coordinates
        @type  c4: [float]        
        """
        v21 = c2 - c1
        v32 = c3 - c2
        L = MU.cross( v21, v32 )
        L_norm = N.sqrt(sum(L**2))

        v43 = c4 - c3
        v23 = c2 - c3
        R = MU.cross( v43, v23 )
        R_norm = N.sqrt(sum(R**2))

        S   = MU.cross( L, R )
        angle = sum( L*R ) / ( L_norm * R_norm )

        if angle > 1.0: angle=1.0
        if angle <-1.0: angle = -1.0

        angle = N.arccos(angle) *180/N.pi
        if sum(S*v32) < 0.0:
            angle = -angle

        return angle

    
    def ramachandran( self, property=None ):
        """
        Create all the ramachandran plot points.
        
        @param property: property to use for coloring of data
        @type  property: array

        @return: list of biggles.Point objects (all the points of the
                 plot)and a biggles.Inset object (property scale).
        @rtype: [ biggles.Point ], biggles.Inset
        """
        p = []

        ## calculate colors and create a legend if a property is given
        if self.profileName:
            palette=CS('plasma', 0,100)
            col = palette.color_array( self.prof )

            legend = Legend(  palette.legend() )
            inset = biggles.Inset((1.1, 0.60), (1.2, .97), legend)

        else:
            col = ['black']*len(self.phi)
            inset = None

        ## add data points to plot
        for i in range(len(self.phi)):
            ## don't add termini - has missing angles
            if self.phi[i] and self.psi[i]:
                if i in self.gly:
                    p += [biggles.Point( self.psi[i], self.phi[i],
                                         type="star", size=1, color=col[i] )]
                elif i in self.pro:
                    p += [biggles.Point( self.psi[i], self.phi[i],
                                         type="filled square", size=1,
                                         color=col[i] )]
                else:
                    p += [biggles.Point( self.psi[i], self.phi[i],
                                         type="filled circle", size=1,
                                         color=col[i] )]
        return p, inset


    def ramachandran_background( self ):
        """
        Creates a background (favoured regions) for a ramachandran plot.

        @return: list of biggles.Point objects
        @rtype: [ biggles.Point ]
        """
        bg = []
        mat = biggles.read_matrix( T.projectRoot() +
                                   '/external/biggles/ramachandran_bg.dat')
        x, y = N.shape(mat)
        for i in range(x):
            for j in range(y):
                if mat[i,j]<200:
                    a = (360./y)*j    - 180
                    b = (360./x)*(x-i)- 180
                    bg += [ biggles.Point( a, b, type="dot" )]
        return bg


    def show( self ):
        """
        Show ramachandran plot.
        """
        plot = biggles.FramedPlot()
        plot.xrange = (-180., 180.)
        plot.yrange = (-180., 180.)
        plot.xlabel = "$\Phi$"
        plot.ylabel = "$\Psi$"      
        
        if self.name:
            plot.title = self.name

        ## add allowed regions
        bg_plot = self.ramachandran_background( )
        for p in bg_plot:
            plot.add( p )

        ## add ramachandran phi, psi valies
        points, inset = self.ramachandran(  )
        for p in points:
            plot.add(p)
        if inset:
            plot.add( inset )
      
        plot.add( biggles.PlotLabel( 1.14, 0.55, self.profileName, size=2) )
        plot.add( biggles.PlotLabel( 1.1, 0.45, "GLY star", size=2) )
        plot.add( biggles.PlotLabel( 1.12, 0.40, "PRO square", size=2) )
      
        plot.show()


## TEST ##

if __name__ == '__main__':

    import glob

    f = glob.glob( T.testRoot()+'/lig_pcr_00/pcr_00/*_1_*pdb' )[:2]
    m = [ PDBModel(i) for i in f ]
    m = [ i.compress( i.maskProtein() ) for i in m ]

    r = Ramachandran( m , name='test')
    r.show()
    
