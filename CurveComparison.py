import dxfgrabber
import numpy as np
from numpy.linalg import inv

#Getting Bezier (currently from whole polyline)
def GetBezier(PointCnt, Points):
    M = MMtrx(PointCnt)
    P = PMtrx(PointCnt, Points)
    T = TMtrx(PointCnt, P)
    C = CMtrx(PointCnt, M, P, T)
    for i in range(min(PointCnt, 8)):
        m=1
#Get M matrix, where number represents number of Bezier fit-points https://pomax.github.io/bezierinfo/#curvefitting
#Max 8 - degree Bezier to get close enough results
def MMtrx(PointCnt):
    M3 = [[1, 0, 0],[-2, 2, 0],[1, -2, 1]]
    M4 = [[1, 0, 0, 0],[-3, 3, 0, 0],[3, -6, 3, 0],[-1, 3, -3, 1]]
    M5 = [[1, 0, 0, 0, 0],[-4, 4, 0, 0, 0],[6, -12, 6, 0, 0],[-4, 12, -12, 4, 0],[1, -4, 6, -4, 1]]
    M6 = [[1, 0, 0, 0, 0, 0],[-5, 5, 0, 0, 0, 0],[10, -20, 10, 0, 0, 0],[-10, 30, -30, 10, 0, 0],[5, -20, 30, -20, 5, 0],[-1, 5, -10, 10, -5, 1]]
    M7 = [[1, 0, 0, 0, 0, 0, 0],[-6, 6, 0, 0, 0, 0, 0],[15, -30, 15, 0, 0, 0, 0],[-20, 60, -60, 20, 0, 0, 0],[15, -60, 90, -60, 15, 0, 0],[-6, 30, -60, 60, -30, 6, 0],[1, -6, 15, -20, 15, -6, 1]]
    M8 = [[1, 0, 0, 0, 0, 0, 0, 0],[-7, 7, 0, 0, 0, 0, 0, 0],[21, -42, 21, 0, 0, 0, 0, 0],[-35, 105, -105, 35, 0, 0, 0, 0],[35, -140, 210, -140, 35, 0, 0, 0],[-21, 105, -210, 210, -105, 21, 0, 0],[7, -42, 105, -140, 105, -42, 7, 0], [-1, 7, -21, 35, -35, 21, -7, 1]]
    M = [M3, M4, M5, M6, M7]
    if (PointCnt >= 8):
        return M8
    else:
        return M[PointCnt]
#Get fit-point matrix
def PMtrx(PointCnt, Points):
    it = PointCnt / 8   #choose points evenly from polyline point array
    if (it == 0):
        it = 1
    k = 0
    for i in range(min(PointCnt - 1, 8 - 1)):
        P[i] = Point[k]
        k += it
    return P
#get T matrix with parameter values
def TMtrx(PointCnt, P):
    for i in range(min(PointCnt, 8)):   #create point (t) distance array
        d[i] = d[i - 1] + PointDist(P[i-1], P[i])
    for i in range(min(PointCnt, 8)):   #scale points to interval [0..1]
        d[i] /= d[min[PointCnt,8] - 1]
    for i in range(min(PointCnt, 8)):   #fill T matrix
        for j in range(min(PointCnt, 8)):
            T[i][j] = d[i]^j
    return T
#get controlpoint matrix
def CMtrx(PointCnt, M, P, T):
    Tt = np.transpose(T)
    Mi = inv(M)
    C = np.matmul(Tt, T)
    C = inv(C)
    C = np.matmul(Mi, C)
    C = np.matmul(C, Tt)
    C = np.matmul(C, P)
    return C
#Returns distance between two points
def PointDist(ax, ay, bx, by):
    return sqrt((ax-bx)^2+(ay-by)^2)

dxf = dxfgrabber.readfile("svarki_002.dxf")

#type of objects in file
print("type of objects in file")
type=[entity.dxftype for entity in dxf.entities]
print(type)

output = [entity for entity in dxf.entities]

Color = [entity.color for entity in output]
print(Color)
Linetype = [entity.linetype for entity in output if entity.dxftype != 'POINT']
print(Linetype)

#get parameters of objects
#Point
PointCoord = [entity.point for entity in output if entity.dxftype == 'POINT']
#Line
LineStart = [entity.start for entity in output if entity.dxftype == 'LINE']
LineEnd = [entity.end for entity in output if entity.dxftype == 'LINE']
#Circle
CenterPoints = [entity.center for entity in output if entity.dxftype == 'CIRCLE']
Radius = [entity.radius for entity in output if entity.dxftype == 'CIRCLE']
#Ellipse
EllipseCenter = [entity.center for entity in output if entity.dxftype == 'ELLIPSE']
EllipseMajorAxis = [entity.major_axis for entity in output if entity.dxftype == 'ELLIPSE']
#Arc
ArcCenter = [entity.center for entity in output if entity.dxftype == 'ARC']
ArcRadius = [entity.radius for entity in output if entity.dxftype == 'ARC']
ArcStartAngle = [entity.start_angle for entity in output if entity.dxftype == 'ARC']
ArcEndAngle = [entity.end_angle for entity in output if entity.dxftype == 'ARC']
#Polyline
PolylineIsClosed = [entity.is_closed for entity in output if entity.dxftype == 'POLYLINE']
PolylineSplineType = [entity.spline_type for entity in output if entity.dxftype == 'POLYLINE']
PolylinePoints = [entity.points for entity in output if entity.dxftype == 'POLYLINE']
PolylineControlPoints = [entity.control_points for entity in output if entity.dxftype == 'POLYLINE']
PolylineBulge = [entity.bulge for entity in output if entity.dxftype == 'POLYLINE']
PolylineVertexCount = [entity.__len__() for entity in output if entity.dxftype == 'POLYLINE']
#LWPolyline
LWPolylinePoints = [entity.points for entity in output if entity.dxftype == 'LWPOLYLINE']
LWPolylineIsClosed = [entity.is_closed for entity in output if entity.dxftype == 'LWPOLYLINE']
LWPolylinePointsCount = [entity.__len__() for entity in output if entity.dxftype == 'LWPOLYLINE']
#Spline
SplineDegree = [entity.degree for entity in output if entity.dxftype == 'SPLINE']
SplineStartTangent = [entity.start_tangent for entity in output if entity.dxftype == 'SPLINE']
SplineEndTangent = [entity.end_tangent for entity in output if entity.dxftype == 'SPLINE']
SplineControlPoints = [entity.control_points for entity in output if entity.dxftype == 'SPLINE']
SplineFitPoints = [entity.fit_points for entity in output if entity.dxftype == 'SPLINE']
SplineKnots = [entity.knots for entity in output if entity.dxftype == 'SPLINE']
SplineIsClosed = [entity.is_closed for entity in output if entity.dxftype == 'SPLINE']

#Point
for x in range(len(PointCoord)):
    print(" POINT\n")
    print(len(PointCoord))
    print(x)
    print("^")
    print(PointCoord[x])
#Line
for x in range(len(LineStart)):
    print(" LINE\n")
    print("Line start")
    print(LineStart[x])
    print("Line end")
    print(LineEnd[x])
#Circle
for x in range(len(CenterPoints)):
    print(" CIRCLE \n")
    print("Center Point")
    print(CenterPoints[x])
    print("Radius")
    print(Radius[x])
#Polyline
for x in range(len(PolylineIsClosed)):
    print(" POLYLINE\n")
    print("is closed")
    print(PolylineIsClosed[x])
    print("Spline type")
    print(PolylineSplineType[x])
    print("line points")
    print(PolylinePoints[x])
    print("Control points")
    print(PolylineControlPoints[x])
    print("bulge")
    print(PolylineBulge[x])
    print("vertex count:")
    print(PolylineVertexCount[x])
    Bezier = GetBezier(PolylineVertexCount, PolylinePoints)
#LWPolyline
for x in range(len(LWPolylineIsClosed)):
    print(" LWPOLYLINE\n")
    print(LWPolylinePoints)
    print(LWPolylineIsClosed)
    print(LWPolylinePointsCount)
#Ellipse
for x in range(len(EllipseCenter)):
    print(" ELLIPSE\n")
    print(EllipseCenter)
    print(EllipseMajorAxis)
#Arc
for x in range(len(ArcCenter)):
    print(" ARC\n")
    print(ArcCenter)
    print(ArcRadius)
    print(ArcStartAngle)
    print(ArcEndAngle)
#Spline
for x in range(len(SplineIsClosed)):
    print(" SPLINE\n")
    print(SplineDegree)
    print(SplineStartTangent)
    print(SplineEndTangent)
    print(SplineControlPoints)
    print(SplineFitPoints)
    print(SplineKnots)
    print(SplineIsClosed)



#Block
#print("BLOCKS")
#BlockBasepoint= [Block.BlockBasepoint for Block in dxf.entities]
#print(BlockBasepoint)
print("=============")
