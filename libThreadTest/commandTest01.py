# -*- coding: utf-8 -*-
 
import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)

### GEOM component
import GEOM
from salome.geom import geomBuilder
import SALOMEDS
import datetime, math


# add lib path
import inspect, os
scriptFullPathName = inspect.getfile(inspect.currentframe())
scriptDir = os.path.abspath(os.path.dirname(scriptFullPathName))
sys.path.append(scriptDir)
import ebt_3d


geompy = geomBuilder.New(theStudy)

# erase all ------------
nb = theStudy.NewBuilder()
for compName in ["SMESH", "GEOM"]:
  comp = salome.myStudy.FindComponent(compName)
  if comp:
    iterator = salome.myStudy.NewChildIterator(comp)
    while iterator.More():
      sobj = iterator.Value()
      iterator.Next()
      nb.RemoveObjectWithChildren(sobj)

salome.sg. EraseAll()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )



# tests ----------------------------------------------------

MakeFuse    = True # False #


meanThreadRadius = 6.375
pitch = 1.5
height = 15.
bBolt = True # False #
nonThreadRadius = 0 #meanThreadRadius -2.
threadAngle = 90.
nbThreads = 1
spacetBetweenThread =  0 # pitch / float(nbThreads) / 2.
smoothStartLength = 0
smoothEndLength   = pitch
rightHanded = True

screw1 = ebt_3d.makeNutBolt2(geompy,
                     meanThreadRadius = meanThreadRadius,
                     pitch = pitch,
                     height = height,
                     bBolt = bBolt,
                     nonThreadRadius = nonThreadRadius,
                     threadAngle = threadAngle,
                     nbThreads = nbThreads,
                     spacetBetweenThread = spacetBetweenThread,
                     smoothStartLength = smoothStartLength,
                     smoothEndLength   = smoothEndLength,
                     rightHanded = rightHanded,
                     MakeFuse = MakeFuse)

geompy.addToStudy(screw1, "screw1")

#-----------------------------------------------------------

meanThreadRadius = 6.375
pitch = 15
height = 15.
bBolt = True # False #
nonThreadRadius = meanThreadRadius -2.
threadAngle = 120
nbThreads = 8
spacetBetweenThread =  pitch / float(nbThreads) / 2.
smoothStartLength = 2
smoothEndLength   = 2
rightHanded = False

screw2 = ebt_3d.makeNutBolt2(geompy,
                     meanThreadRadius = meanThreadRadius,
                     pitch = pitch,
                     height = height,
                     bBolt = bBolt,
                     nonThreadRadius = nonThreadRadius,
                     threadAngle = threadAngle,
                     nbThreads = nbThreads,
                     spacetBetweenThread = spacetBetweenThread,
                     smoothStartLength = smoothStartLength,
                     smoothEndLength   = smoothEndLength,
                     rightHanded = rightHanded,
                     MakeFuse = MakeFuse)

geompy.addToStudy(screw2, "screw2")


#-----------------------------------------------------------

meanThreadRadius = 6.375
pitch = 1.5
height = 5
bBolt = False #
nonThreadRadius = meanThreadRadius +3.
threadAngle = 60
nbThreads = 1
spacetBetweenThread =  0
smoothStartLength = 1
smoothEndLength   = 1
rightHanded = True

screw3 = ebt_3d.makeNutBolt2(geompy,
                     meanThreadRadius = meanThreadRadius,
                     pitch = pitch,
                     height = height,
                     bBolt = bBolt,
                     nonThreadRadius = nonThreadRadius,
                     threadAngle = threadAngle,
                     nbThreads = nbThreads,
                     spacetBetweenThread = spacetBetweenThread,
                     smoothStartLength = smoothStartLength,
                     smoothEndLength   = smoothEndLength,
                     rightHanded = rightHanded,
                     MakeFuse = MakeFuse)

geompy.addToStudy(screw3, "screw3")

#-----------------------------------------------------------

