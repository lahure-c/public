# -*- coding: utf-8 -*-

# pour recharger la bibliothèque depuis la console: reload(<nom_biblio>)
# exemple: reload(ebt_3d)

import sys
import salome
import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS
import datetime

def GetSolidList_OrderedFromDistToPoint(geompy, CompoundSolid, RefPoint):
  """  Measure the minimal distance of all solids indide CompoundSolid
       returns the list of these solids ordered from min to max distance
       Returns a list of solid GEOM objects
       this function was made to replace the buggy function geompy.MinDistanceComponents() with the "SOLID" type
       here the min distance is the min of the faces' min distances
  """
  UnorderedSolidList = geompy.ExtractShapes(CompoundSolid, geompy.ShapeType["SOLID"],isSorted = False)
  distancesS = []   # tupple: [ (solid object, distanceToPoint), ...]
  for solid in UnorderedSolidList:
    UnorderedFaceList = geompy.ExtractShapes(solid, geompy.ShapeType["FACE"],isSorted = False)
    distancesF = []  # list of distances of faces
    for face in UnorderedFaceList:    
      distance_aDistDxDyDz = geompy.MinDistanceComponents(RefPoint, face)  # distance min de chaque face du solide
      distancesF.append(distance_aDistDxDyDz[0])             
    distancesS.append((solid,min(distancesF)))                             # objet solide et sa distance min au point

  OrderedSolidTuppleList = sorted(distancesS, key=lambda id: id[1])        # tri suivant la colonne d'indice 1
  OrderedSolidList = [x[0] for x in OrderedSolidTuppleList]                # extraction de la colonne d'indice 0
  return OrderedSolidList







def makeNutBolt2(geompy, meanThreadRadius, pitch, height, bBolt=True, nonThreadRadius=0., threadAngle= 90., nbThreads=1., spacetBetweenThread=0., smoothStartLength=0, smoothEndLength=0, rightHanded=True, MakeFuse=False):
  """Make a bolt or a nut with a 90° thread angle
  meanThreadRadius is the mean radius of the thread = x2 (rayon moyen de la vis)
  pitch is the pitch of the screw, the advance length after one turn = z2 (pas de la vis)
  height is the length of the bolt = z1 (longueur de la vis)
  bBolt: boolean value,  True makes a bolt, False makes a nut  (True -> vis, False -> Ecrou)
  nonThreadRadius: outside radius of the nut, or radius of the hole inside bolt (0 for a plain bolt) = x4  
  threadAngle  in degrees = a1 (angle du filet en degrés
  spacetBetweenThread = vertical distance in between thread = x4 
  nbThreads: number of threads (nombre de filets)
  smoothEndTurns       fade off length of the end of the threads = z5 (longueur d'évanescence des filets)
  smoothStartLength    like smoothEndTurns
  rightHanded: True for right handed thread, False for left handed   (pas à droite où à gauche)
  MakeFuse if True fuses all the components of the screw in a single solid, otherwise return a compound of solids (True -> fusionne la vis)
  returns a SOLID GEOM_Object
  """

  dateStart = datetime.datetime.now()
  outputDebug = False
  #outputDebug = True
  if outputDebug:
    import inspect, os
    scriptFullPathName = inspect.getfile(inspect.currentframe())
    scriptDir = os.path.abspath(os.path.dirname(scriptFullPathName))
    DebugFileName = scriptDir + "/DebugFile_" + str(dateStart);
    print DebugFileName
    DebugFile = open(DebugFileName, "w")

  x2 = meanThreadRadius
  z2 = pitch
  z1 = height
  x4 = nonThreadRadius
  a1 = threadAngle /180.*math.pi
  z5start = smoothStartLength
  z5end   = smoothEndLength
  z3 = z2 / float(nbThreads)

  z4 = spacetBetweenThread

  z6 = z3 - z4

  # what is x1 and x3
  # x1 -x3 = z6/2 / tan(a1/2)
  x1Mx3 = z6/2. / math.tan(a1/2.)

  x1 = x2 + x1Mx3/2.
  x3 = x2 - x1Mx3/2.

  zeroDistance = 1e-7

  if outputDebug:
    print "x1= ", x1, "  x2= ", x2, "  x3= ", x3, "  x4= ", x4, "  z1= ", z1, "  z2= ", z2, "  z3= ", z3, "  z4= ", z4, "  z6= ", z6, "  nbThreads= ", nbThreads


  fuseBuggDistance1 = 0.001     # should be zero
  bFuseBuggTrick1  = True       # True # False #   should be False. True takes one more boolean operation
  bFuseBuggTrick2  = True       # False may bugg. True is safer and faster (the final fuse is splitted)

  '''
  Bugg (no solid output) with:
  meanThreadRadius = 6.375
  pitch = 1.5
  height = 15.
  bBolt = False
  nonThreadRadius = meanThreadRadius +2.
  threadAngle= 90.
  nbThreads=1
  smoothStartLength = 0
  smoothEndLength   = pitch
  rightHanded = True
  MakeFuse    = True
  spacetBetweenThread = 0.04  # or 0 ?
  fuseBuggDistance1 = 0.0001
  bFuseBuggTrick1  = True
  bFuseBuggTrick2  = True or False
  #
  solved with fuseBuggDistance1 = 0.001   time 18s
  solved with bFuseBuggTrick1  = False     time  6s
  solved with bFuseBuggTrick1=False  and  fuseBuggDistance1 = 0.001    time 5s
  '''

  '''
  Bugg (no thread in solid output) with:
  meanThreadRadius = 6.375
  pitch = 1.5
  height = 15.
  bBolt = False
  nonThreadRadius = meanThreadRadius +2.
  threadAngle= 90.
  nbThreads=1
  smoothStartLength = pitch
  smoothEndLength   = pitch
  rightHanded = True
  MakeFuse    = True
  spacetBetweenThread = 0.04  # or 0 ?
  fuseBuggDistance1 = 0.0001
  bFuseBuggTrick1  = False
  bFuseBuggTrick2  = True or False
  #
  solved with fuseBuggDistance1 = 0.001   time  5s
  solved with bFuseBuggTrick1  = True     time 19s
  '''

  '''
  Bugg (no thread in solid output  or Boolean operation aborted : non valid shape result) with:
  meanThreadRadius = 4.7940000000000005
  pitch = 1.5
  height = 8.
  bBolt = False
  nonThreadRadius = 7.2940000000000005
  threadAngle= 90.
  nbThreads=1
  smoothStartLength = 0
  smoothEndLength   = 0.8
  rightHanded = True
  MakeFuse    = True
  spacetBetweenThread = 0 # or 0.001 
  fuseBuggDistance1 = 0.001
  bFuseBuggTrick1  = False
  bFuseBuggTrick2  = True
  #
  solved with spacetBetweenThread = 0.1
  solved with bFuseBuggTrick1  = True
  '''

  # testing the input values -------------------------
  if x2 <0 or z2 <0 or z1<0 or nbThreads<1 or z5start<0 or z5end<0:
    print sys._getframe().f_code.co_name, ": error in the parameters: "
    return None

  # check x4
  if bBolt:
    if(x4 < 0):
      x4 = 0
    if(x4 > x2):
      x4 = x2 -1.
      print "hole radius of the nucleus cannot be greater than x2 it to: ", x2 -1

  else:  # Nut
    if(x4 < x1):
      x4 = x1 + 1.
      print "external radius of the nut cannot be lower than x1 Setting it to: ", x1 +1.



  # éléments de base -----------------------------
  O = geompy.MakeVertex(0, 0, 0)
  OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
  OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
  OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
  if outputDebug: geompy.addToStudy( O, 'O' )
  if outputDebug: geompy.addToStudy( OX, 'OX' )
  if outputDebug: geompy.addToStudy( OY, 'OY' )
  if outputDebug: geompy.addToStudy( OZ, 'OZ' )


  # if bFuseBuggTrick1 ==True, we will use dzNonThreadBooleanBugg
  dzNonThreadBooleanBugg = 1.
  coefDzNonThreadBooleanBugg = 1 + 2*dzNonThreadBooleanBugg/z1    # coefficient for scaling the nonThread part along the Z axis


  VertexBN_Bottom = geompy.MakeVertex(0, 0, 0)
  VertexBN_Middle = geompy.MakeVertex(0, 0, z1/2.)
  # if outputDebug: geompy.addToStudy(VertexBN_Middle, 'VertexBN_Middle')
  VertexBN_Top    = geompy.MakeVertex(0, 0, z1)
  VectBN_Axis1  = geompy.MakeVector(VertexBN_Bottom, VertexBN_Top)
  # if outputDebug: geompy.addToStudy(VectBN_Axis1, 'VectBN_Axis1')
  VectBN_Axis2  = geompy.MakeVector(VertexBN_Top, VertexBN_Bottom)
  # if outputDebug: geompy.addToStudy(VectBN_Axis2, 'VectBN_Axis2')
  # calcul de la taille du plan de coupe
  planeSize = 3 * (x2+z2)
  if not bBolt and x4 > x2+z2:
    planeSize = 3 * x4
  PlaneBN_Top = geompy.MakePlane(VertexBN_Top, VectBN_Axis1, planeSize)
  PlaneBN_Bottom = geompy.MakePlane(VertexBN_Bottom, VectBN_Axis1, planeSize)
  # if outputDebug: geompy.addToStudy( VertexBN_Middle, 'Vertex_Middle' )
  # if outputDebug: geompy.addToStudy( VertexBN_Top, 'Vertex_Top' )
  if outputDebug: geompy.addToStudy( PlaneBN_Top, 'Plane_Top' )
  if outputDebug: geompy.addToStudy( PlaneBN_Bottom, 'Plane_Bottom' )



  # génération de la spirale pour servir de chemin d'extrusion ------------------------------
  if outputDebug: print "making spiral path"
  nbPtsParTour = 512 # or 256.  Number of calculated points per turn
  P= []      # list of points of the path
  iCpt =0
  sensPas = 1.0    # sens du pas: 1 si à droite; -1 si à gauche
  if not rightHanded:
    sensPas = -1.0

  smoothingDeltaRadius = math.fabs(x1-x3)  # the size of the thread
  extraTurns = 2.                          # number of extra spiral turns outside before and after the normal thread, should be an integer

  while True:
    # height of the spiral point
    Pz = -extraTurns*z2 + iCpt*z2/float(nbPtsParTour)  # start extraTurns before in order to prepare for a smooth start
    #  Pz = iCpt*z2/float(nbPtsParTour)  # start @ 0
    if Pz  > z1 + extraTurns* z2:                # end + extraTurns
      break

    angle = sensPas * 2*math.pi * iCpt / float(nbPtsParTour)  # angle of the spiral point

    # radius of the spiral point, normally x2, but different if smoothing is required  
    # in case of a smooth thread, we want the thread size to increase from 0 to its nominal size
    nominalSpiralRadius = x2 # the nominal value (without smoothing)
    # we set the smoothed size of the spiral radius of x2 +/- z3
    if bBolt:
      smoothedSpiralRadius = x2 - smoothingDeltaRadius
    else:
      smoothedSpiralRadius = x2 + smoothingDeltaRadius

    spiralRadius = nominalSpiralRadius   # default radius is the nominal value

    if z5start > 0:
      # in that case the sketch needs to be translated to abscisse of smoothedSpiralRadius instead of nominalSpiralRadius
      if Pz < 0:
        spiralRadius = smoothedSpiralRadius
      elif Pz < z5start:
        alpha = math.pi * Pz / z5start  # varie de 0 à pi
        spiralRadius = smoothedSpiralRadius + (nominalSpiralRadius - smoothedSpiralRadius) *  0.5 * (1 - math.cos(alpha))
        #if outputDebug: print spiralRadius

    if z5end > 0:
      if Pz > z1:
        spiralRadius = smoothedSpiralRadius
      elif Pz > z1 - z5end:
        alpha = math.pi * (Pz - (z1 - z5end)) / z5end  # varie de 0 à pi
        spiralRadius = nominalSpiralRadius - (nominalSpiralRadius - smoothedSpiralRadius) *  0.5 * (1 - math.cos(alpha))
        # if outputDebug: print spiralRadius

    Px = spiralRadius * math.cos(angle)
    Py = spiralRadius * math.sin(angle)
    Pt = geompy.MakeVertex(Px, Py, Pz)
    # if outputDebug: geompy.addToStudy(Pt,  'P['+str(iCpt)+']')
    P.append(Pt)                               #  List of points
    iCpt+=1

  SpiralPath = geompy.MakeInterpol(P, False)  # Path with extra length
  if outputDebug: geompy.addToStudy(SpiralPath,'SpiralPath')
  #---------------------------------------------------


  #----- thread profile ----------------------------------------
  if outputDebug: print "making sketch"
  # Generation of bolt male thread profile
  sk = geompy.Sketcher3D()
  sk.addPointsAbsolute(x3, 0, -z6/2.)
  sk.addPointsAbsolute(x3, 0,  z6/2.)
  sk.addPointsAbsolute(x1, 0,  0.)
  sk.close()

  Sketch_Bolt = sk.wire()
  if outputDebug: geompy.addToStudy(Sketch_Bolt, 'Sketch_Bolt')

  # in the case of a smooth start the sketch needs to be at a the lower radius
  # then the spiral will grow in radius
  # if there's no smooth start (even with a smooth end) the sketch is at the nominal radius
  if z5start > 0:
    geompy.TranslateDXDYDZ(Sketch_Bolt, -smoothingDeltaRadius, 0, 0)
  # ---------------------------------------------------------------------------


  # génération du croquis d'extrusion ----------------------
  ExtrudedThreadList = []
  if bBolt:
    TheSketch = Sketch_Bolt
  else:  
    # the nut female profile is the bolt male thread profile symmetry along the x2 axis + translation
    Vertex_mR_Bottom = geompy.MakeVertex(x2, 0, 0)
    Vertex_mR_Top    = geompy.MakeVertex(x2, 0, z1)
    Vect_mR_Axis     = geompy.MakeVector(Vertex_mR_Bottom, Vertex_mR_Top)
    if outputDebug: geompy.addToStudy(Vect_mR_Axis, 'Vect_mR_Axis')
    Sketch_Nut_mirrored = geompy.MakeMirrorByAxis(Sketch_Bolt, Vect_mR_Axis)
    if outputDebug: geompy.addToStudy( Sketch_Nut_mirrored, "Sketch_Nut_mirrored")
    TheSketch = Sketch_Nut_mirrored

  if outputDebug: geompy.addToStudy(TheSketch,'TheSketch')

  # génération des extrusions ------------------
  for iThread in range(int(nbThreads)):
    if outputDebug: print "Extruding thread ", iThread
    TheSketch_ = geompy.MakeRotation(TheSketch, OZ, 2*math.pi*iThread/float(nbThreads))
    if outputDebug: geompy.addToStudyInFather(TheSketch, TheSketch_,'TheSketch_'+ str(iThread))

    ThreadFace_  = geompy.MakeFaceWires([TheSketch_], 1)
    if outputDebug: geompy.addToStudyInFather(TheSketch, ThreadFace_,'ThreadFace_'+ str(iThread))

    SpiralPath_ = geompy.MakeRotation(SpiralPath, OZ, 2*math.pi*iThread/float(nbThreads))
    if outputDebug: geompy.addToStudyInFather(SpiralPath, SpiralPath_,'SpiralPath_'+ str(iThread))

    ExtrudedThreadList.append(geompy.MakePipeBiNormalAlongVector(ThreadFace_, SpiralPath_, VectBN_Axis1))

  ExtrudedThreadCompound = geompy.MakeCompound(ExtrudedThreadList)
  if outputDebug: geompy.addToStudy(ExtrudedThreadCompound, 'ExtrudedThreadCompound')
  ThreadShellList = geompy.ExtractShapes(ExtrudedThreadCompound, geompy.ShapeType["SHELL"],isSorted = False)
  icpt=0
  for shell in ThreadShellList:
    if outputDebug: geompy.addToStudyInFather(ExtrudedThreadCompound, shell, 'shell-'+ str(icpt))
    icpt +=1
  #--------------------------------------------------- 


  # pour le rajout du flanc du thread. la partie cylindrique intérieure du filet pour la vis, extérieure pour l'écrou
  CylFF = geompy.MakeCylinderRH(x3, z1 + 6*z2)
  geompy.TranslateDXDYDZ(CylFF, 0, 0, -3*z2)  # décalage du cylindre pour accélérer le partionnement: absence de surfaces coplanaires accélère le partitionnement
  if outputDebug: geompy.addToStudy( CylFF, 'CylFF' )

  CylFF_revList = geompy.GetShapesOnCylinder(CylFF, geompy.ShapeType["FACE"], OZ, x3, GEOM.ST_ON)
  CylFF_rev = geompy.MakeCompound(CylFF_revList)
  if outputDebug: geompy.addToStudy( CylFF_rev, 'CylFF_rev' )
  #--------------------------------------------------------

  SolidThreadList = []
  for iThread in range(int(nbThreads)):    # for each thread

    if outputDebug: print "Partitionning the extrusion ", iThread
    # partionnement suivant les plans horizontaux
    ExtrudedThreadPartition = geompy.MakePartition([ExtrudedThreadList[iThread]], [PlaneBN_Bottom, PlaneBN_Top], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
    if outputDebug: geompy.addToStudyInFather(ExtrudedThreadCompound, ExtrudedThreadPartition, 'ExtrudedThreadPartition-'+ str(iThread))

    if outputDebug: print "Getting the middle solid of extrusion ", iThread
    ordSList = GetSolidList_OrderedFromDistToPoint(geompy, ExtrudedThreadPartition, VertexBN_Middle)

    if outputDebug:
      icpt=0
      for solid in ordSList:
        geompy.addToStudyInFather(ExtrudedThreadPartition, solid, 'orderedSolid-'+ str(iThread) + '-' + str(icpt))
        icpt +=1

    if len(ordSList) != 3 :
      print sys._getframe().f_code.co_name, " Error in ExtrudedThreadPartition, instead of 3 solids, got ",len(ordSList)
      return None 

    MiddleSolid = ordSList[0]   # le solide le plus proche du point
    if outputDebug: geompy.addToStudyInFather(ExtrudedThreadPartition, MiddleSolid, 'MidlleSolid-'+ str(iThread))
    #---------------------

    SolidThreadList.append(MiddleSolid)   # stores the object


  # ------------------------------------------------------

  if outputDebug:
    print "Making final thread "
    DebugFile.write("Making final thread\n")


  if bBolt:
    CylPlain = geompy.MakeCylinderRH(x3 +fuseBuggDistance1, z1)  # plain part of the screw
    if x4 > zeroDistance:                         # case of a hole in the nucleus
      CylCut = geompy.MakeCylinderRH(x4, z1)   # hollow part of the screw is CylCut
      CylNucleus = geompy.MakeCut(CylPlain,CylCut)
    else:
      CylNucleus = CylPlain

    if MakeFuse and bFuseBuggTrick1:   # we make it taller to avoid fuse buggs - we'll cut it after the fuse
      CylNucleus = geompy.MakeScaleAlongAxes(CylNucleus, VertexBN_Middle, 1, 1, coefDzNonThreadBooleanBugg)
    if outputDebug: geompy.addToStudy(CylNucleus, 'CylNucleus')
    SolidThreadList.insert(0, CylNucleus)   # in the first place for helping the use operation

  else:  # nut
    CylPlain = geompy.MakeCylinderRH(x4, z1)      
    if outputDebug: geompy.addToStudy(CylPlain, "CylPlain")

    CylCut = geompy.MakeCylinderRH(x1 -fuseBuggDistance1, z1)
    if outputDebug: geompy.addToStudy(CylCut, "CylCut")

    CylPipe = geompy.MakeCut(CylPlain, CylCut)
    if MakeFuse and bFuseBuggTrick1:   # we make it taller to avoid fuse buggs - we'll cut it after the fuse
      CylPipe = geompy.MakeScaleAlongAxes(CylPipe, VertexBN_Middle, 1, 1, coefDzNonThreadBooleanBugg)

    if outputDebug: geompy.addToStudy(CylPipe, 'CylPipe')
    SolidThreadList.insert(0, CylPipe)      # in the first place for helping the fuse operation

  FinalSolid = geompy.MakeCompound(SolidThreadList)       


  # fuse all the components together
  if MakeFuse:
    if outputDebug: print "Fusing threads with body"

    if len(SolidThreadList) <=2 or not bFuseBuggTrick2:   # when there are only two objects or we don't care about the bugg, we use the simplest fuse method
      # fuse method: only one fuse (simplest)
      if outputDebug:
        date01 = datetime.datetime.now()
        DebugFile.write("only one fuse\n")
        DebugFile.flush()
      FinalSolid = geompy.MakeFuseList(SolidThreadList)
      if outputDebug:
        date02 = datetime.datetime.now()
        DebugFile.write("  time spent= {}\n" .format(date02-date01))
        DebugFile.flush()

    else:   # whith multiple thread it may bugg, we need a different fuse method
      if True:
        # fuse method: even and odd (fastest)
        # fuse is made in 3 times: 1: nucleus + even threads; 2: nucleus + odd threads; 3: fuse of the the 2 fuses
        # this way the fuse avoid fusing contiguous threads
        evenList=[]
        oddList=[]
        evenList.append(SolidThreadList[0])
        oddList.append(SolidThreadList[0])
        ii = 1
        for ii in range(1, len(SolidThreadList)):
          if ii%2==0:
            evenList.append(SolidThreadList[ii])
          else:
            oddList.append(SolidThreadList[ii])

        if outputDebug:
          date01 = datetime.datetime.now()
          DebugFile.write("even fuse operation started, numb of threads={}\n" .format(len(evenList)))
          DebugFile.flush()
        evenFuse = geompy.MakeFuseList(evenList, False, False)
        if outputDebug:  geompy.addToStudy(evenFuse, "evenFuse")

        if outputDebug:
          date02 = datetime.datetime.now()
          DebugFile.write("  time spent= {}\n" .format(date02-date01))
          DebugFile.write("odd fuse operation started, numb of threads={}\n" .format(len(oddList)))
          DebugFile.flush()
        oddFuse  = geompy.MakeFuseList(oddList, False, False)
        if outputDebug:  geompy.addToStudy(oddFuse, "oddFuse")

        if outputDebug:
          date03 = datetime.datetime.now()
          DebugFile.write("  time spent= {}\n" .format(date03-date02))
          DebugFile.write("even and odd fuse operation started\n")
          DebugFile.flush()
        FinalSolid = geompy.MakeFuseList([evenFuse, oddFuse], False, False)
        if outputDebug:
          date04 = datetime.datetime.now()
          DebugFile.write("  time spent= {}\n" .format(date04-date03))
          DebugFile.flush()

      else:
        # fuse method: two by two (was the safest ?)
        # fuse is made progressively, two by two, in order to avoid the bugs
        FinalSolid = SolidThreadList[0]
        for ii in range(1, len(SolidThreadList)):
          if outputDebug:
            dateIterationStart = datetime.datetime.now()
            DebugFile.write("2 by 2 fuse operation iteraton num= {}\n" .format(ii))
            DebugFile.flush()

          FinalSolid = geompy.MakeFuse(FinalSolid, SolidThreadList[ii], False, False)

          if outputDebug:
            dateIterationStop = datetime.datetime.now()
            DebugFile.write("time spent= {}\n" .format(dateIterationStop-dateIterationStart))
            DebugFile.write("\n")
            DebugFile.flush()

  # fin bMakeFuse


  if MakeFuse and bFuseBuggTrick1:
    CylPlainReal = geompy.MakeCylinderRH(max([x4 +1., x2 + z3/4.0 +1.]), z1)  # cut along z axis only
    if outputDebug:  geompy.addToStudy(FinalSolid, 'FinalSolidbeforeFuse')
    if outputDebug:  geompy.addToStudy(CylPlainReal, 'CylPlainReal')      
    FinalSolid = geompy.MakeCommonList([CylPlainReal, FinalSolid])


  if outputDebug:
    geompy.addToStudy(FinalSolid, 'FinalSolid')
    ii = 0
    for solid in SolidThreadList:
      geompy.addToStudyInFather(FinalSolid, solid, 'solid_'+ str(ii))
      ii+=1


  dateStop = datetime.datetime.now()
  if outputDebug:
    DebugFile.write("Fin, temps de calcul = {0}\n" .format(dateStop - dateStart))
    DebugFile.close()

  return FinalSolid



