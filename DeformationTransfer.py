#
"""
DeformationTransfer
    "Deformation Transfer" for Maya 2018

    Reference:
        - Deformation Transfer, SIGGRAPH 2005
        
"""

__author__ = "Naoya Iwamoto <iwanao731@gmail.com>"
__status__ = "beta release"
__version_ = "0.1"
__date__ = "15 Feb 2019"

import sys
import maya.api.OpenMaya as om
import maya.cmds as cmds
import maya.mel
import numpy as np
import scipy.linalg as linalg

def maya_useNewAPI():
    pass

def initializePlugin(plugin):
    fnPlugin = om.MFnPlugin(plugin, vendor = 'Euclid Lab.', version = '0.1')
    try:
        fnPlugin.registerCommand(DeformationTransferCmd.commandName, DeformationTransferCmd.creator)
    except:
        sys.stderr.write('Failed to register command: {0}\n'.format(DeformationTransferCmd.commandName))
        raise

    cmds.setParent('MayaWindow')
    if not cmds.menu('EuclidLab', query=True, exists=True):
        cmds.menu('EuclidLab', label='EuclidLab')
    cmds.setParent('EuclidLab', menu=True)
    if not cmds.menu('Deformation Transfer', query=True, exists=True):
        cmds.menuItem('Deformation Transfer', label='Deformation Transfer', tearOff=True, command='import maya.mel;maya.mel.eval("DeformationTransferBuild")')

def uninitializePlugin(plugin):
    fnPlugin = om.MFnPlugin(plugin)
    try:
        fnPlugin.deregisterCommand(DeformationTransferCmd.commandName)
    except:
        sys.stderr.write('Failed to unregister command: {0}\n'.format(DeformationTransferCmd.commandName))
        raise
    cmds.deleteUI('MayaWindow|EuclidLab|Deformation Transfer', menuItem=True)

class DeformationTransferCmd(om.MPxCommand):
    commandName = 'DeformationTransferBuild'

    def __init__(self):
        om.MPxCommand.__init__(self)

    @staticmethod
    def creator():
        return DeformationTransferCmd()

    def doIt(self, args):
        print '---------------------'
        print 'Deformation Transfer'
        print '---------------------'
        DeformationTransfer()

def getSelectMesh():
    slist = om.MGlobal.getActiveSelectionList()
    itsl = om.MItSelectionList(slist)
    meshPaths = []
    while not itsl.isDone():
        dagPath = itsl.getDagPath()
        itsl.next()
        if dagPath is None:
            continue
        apiType = dagPath.apiType()
        if apiType != om.MFn.kTransform:
            continue
        for c in xrange(dagPath.childCount()):
            child = dagPath.child(c)
            if child.apiType() != om.MFn.kMesh:
                continue
            path = dagPath.getAPathTo(child)
            mesh = om.MFnMesh(path)
            if not mesh.findPlug('intermediateObject', True).asBool():
                meshPaths.append(path)
                break
    return meshPaths
    
def setMatrixCol(mat, vec, iCol):
    mat[0][iCol] = vec[0]
    mat[1][iCol] = vec[1]
    mat[2][iCol] = vec[2]

def computeMatrixA(meshPath):
    mesh = om.MFnMesh(meshPath)
    #numVertices, vertexList = mesh.getVertices()

    #print "print numVertices : ", mesh.numVertices

    # get triangle information
    ___, indices = mesh.getTriangles()

    neighbor = []
    offset = len(neighbor)
    neighbor = neighbor + [set() for v in xrange(mesh.numVertices)]

    MatA = np.zeros((len(indices), mesh.numVertices))

    for triID in xrange(len(indices) / 3):
        i0 = indices[triID * 3 + 0] + offset
        i1 = indices[triID * 3 + 1] + offset
        i2 = indices[triID * 3 + 2] + offset

        e0 = mesh.getPoint(i1) - mesh.getPoint(i0);
        e1 = mesh.getPoint(i2) - mesh.getPoint(i0);

        # Construct Va
        Va = np.zeros((3, 2));
        setMatrixCol(Va, e0, 0);
        setMatrixCol(Va, e1, 1);

        # QR Decomposition (TBD)        
        Q, R = np.linalg.qr(Va)

        invRQT = np.dot(np.linalg.inv(R), Q.transpose())

        # i0
        MatA[triID * 3 + 0][i0] = - invRQT[0][0] - invRQT[1][0] 
        MatA[triID * 3 + 1][i0] = - invRQT[0][1] - invRQT[1][1] 
        MatA[triID * 3 + 2][i0] = - invRQT[0][2] - invRQT[1][2] 

        # i1
        MatA[triID * 3 + 0][i1] = invRQT[0][0]
        MatA[triID * 3 + 1][i1] = invRQT[0][1]
        MatA[triID * 3 + 2][i1] = invRQT[0][2]

        # i2
        MatA[triID * 3 + 0][i2] = invRQT[1][0]
        MatA[triID * 3 + 1][i2] = invRQT[1][1]
        MatA[triID * 3 + 2][i2] = invRQT[1][2]

    return MatA

def calcNormal(v1, v2, v3):
    return np.cross((v2-v1),(v3-v1))

def setTriEdgeMatrix(mesh, i1, i2, i3):

    matV = np.zeros((3,3))

    v1 = mesh.getPoint(i1);
    v2 = mesh.getPoint(i2);
    v3 = mesh.getPoint(i3);

    e1 = v2 - v1;
    e2 = v3 - v1;
    e3 = calcNormal(v1, v2, v3);

    setMatrixCol(matV, e1, 0);
    setMatrixCol(matV, e2, 1);
    setMatrixCol(matV, e3, 2);

    return matV

def setMatrixBlock(matF, matBlock, irow, icol):

    rows, cols = matBlock.shape

    for r in range(rows):
        for c in range(cols):
            matF[irow+r][icol+c] = matBlock[r][c]

def computeMatrixF(srcRefMeshPath, srcDefMeshPath):

    src_ref_mesh = om.MFnMesh(srcRefMeshPath)
    src_def_mesh = om.MFnMesh(srcDefMeshPath)
    
    #numVertices, src_vertexList = src_ref_mesh.getVertices()

    # get triangle information
    ___, indices = src_ref_mesh.getTriangles()

    numTriangle = len(indices)/3

    # loop triangle
    neighbor = []
    offset = len(neighbor)
    neighbor = neighbor + [set() for v in xrange(src_ref_mesh.numVertices)]

    # set MatF
    matF = np.zeros((numTriangle*3, 3))
    
    for triID in xrange(numTriangle):

        i0 = indices[triID * 3 + 0] + offset
        i1 = indices[triID * 3 + 1] + offset
        i2 = indices[triID * 3 + 2] + offset

        # 3 x 3 matrix
        Va = setTriEdgeMatrix(src_ref_mesh, i0, i1, i2)
        Vb = setTriEdgeMatrix(src_def_mesh, i0, i1, i2)

        # QR Decomposition
        Q, R = np.linalg.qr(Va)

        invRQT = np.dot(np.linalg.inv(R), Q.transpose())

        Sa = np.dot(Vb, invRQT)
        SaT = Sa.transpose()

        setMatrixBlock(matF, SaT, triID*3, 0)

    return matF

def DeformationTransfer():

    meshPaths = getSelectMesh()

    src_ref_mesh_path = meshPaths[0]
    src_def_mesh_path = meshPaths[1]
    trg_ref_mesh_path = meshPaths[2]
    trg_def_mesh_path = meshPaths[3]

    # set target reference model
    print "analysis target reference model"
    MatA = computeMatrixA(trg_ref_mesh_path)
    print "MatA : ", MatA.shape

    MatAt = MatA.transpose()
    print "MatAt : ", MatAt.shape

    LU = linalg.lu_factor(np.dot(MatAt, MatA))   # LU decompose 

    # transfer to target model
    print "analysis transfer of source ref and def"
    matF = computeMatrixF(src_ref_mesh_path, src_def_mesh_path) # deformation gradient for source model
    
    print "matF : ", matF.shape

    # solve deformation transfer
    UtS = np.dot(MatAt, matF)

    print "UtS : ", UtS.shape

    trg_def_x = linalg.lu_solve(LU, UtS)

    print "trg_def_x : ", trg_def_x.shape

    # set result to target def model
    trg_def_mesh = om.MFnMesh(trg_def_mesh_path)
    for i in xrange(trg_def_mesh.numVertices):
        pos = trg_def_mesh.getPoint(i)
        pos[0] = trg_def_x[i][0]
        pos[1] = trg_def_x[i][1]
        pos[2] = trg_def_x[i][2]
        trg_def_mesh.setPoint(i, pos)

    print "done"

