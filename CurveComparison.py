import dxfgrabber
import numpy as np
from numpy.linalg import inv
import array as arr

#get factorial
def Fact(n):
    f = 1
    for i in range(n):
        f = f * (i + 1)
    return f
#get binominal coefficients
def BinCoef(n, k):
    return Fact(n)/(Fact(k) * Fact(n - k))
#Getting Bezier (currently from whole polyline)
def GetBezier(PointCnt, Points):
    M = MMtrx(PointCnt)
    P = PMtrx(PointCnt, Points)
    T = TMtrx(PointCnt, P)
    C = CMtrx(PointCnt, M, P, T)
    Bx = ""
    By = ""
    n = min(PointCnt, 8)
    for i in range(1, n + 1):
        BinC = BinCoef(n - 1, i - 1)
        Bx = f"{Bx} + {BinC} * (1 - t)^{n - i}*t^{i - 1}({P[i - 1][0]}) " #Create Bx(t) formula as a string
        By = f"{By} + {BinC} * (1 - t)^{n - i}*t^{i - 1}({P[i - 1][1]}) " #Create Bx(t) formula as a string
    return (f"({Bx},{By})")
        #'+' BinC '* (1 - t) *' (1 - t)**(n - i) * t**i * P[i - 1][i - 1])
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
        return M[PointCnt - 3]
#Get fit-point matrix
def PMtrx(PointCnt, Points):
    it = (PointCnt - 1) / 7   #choose points evenly from polyline point array
    if (it == 0):
        it = 1
    P = []
    k = 0
    for i in range(min(PointCnt, 8) - 1):
        P.append(Points[int(round(k))])
        k += it
    P.append(Points[-1])
    return P
#get T matrix with parameter values
def TMtrx(PointCnt, P):
    d = []
    d.append(0)
    for i in range(1, min(PointCnt, 8)):   #create point (t) distance array
        dist = d[i - 1] + PointDist(P[i-1][0], P[i][0], P[i-1][1], P[i][1])
        d.append(dist)
    for i in range(min(PointCnt, 8)):   #scale points to interval [0..1]
        d[i] = d[i] / d[-1]
    T = []
    for i in range(min(PointCnt, 8)):   #fill T matrix
        T.append([])
        for j in range(min(PointCnt, 8)):
            T[i].append(d[i]**j)
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
    return np.sqrt((ax-bx)**2+(ay-by)**2)

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
for entity in output:
    #Point
    if entity.dxftype == 'POINT':
        PointCoord = entity.point
#        for x in range(len(PointCoord)):
#            print(" POINT\n")
#            print(len(PointCoord))
#            print(PointCoord[x])
    #Line
    if entity.dxftype == 'LINE':
        LineStart = entity.start
        LineEnd = entity.end
#        for x in range(len(LineStart)):
#            print(" LINE\n")
#            print("Line start")
#            print(LineStart[x])
#            print("Line end")
#            print(LineEnd[x])
    #Circle
    if entity.dxftype == 'CIRCLE':
        CenterPoints = entity.center
        Radius = entity.radius
        for x in range(len(CenterPoints)):
            print(" CIRCLE \n")
            print("Center Point")
            print(CenterPoints[x])
            print("Radius")
            print(Radius[x])
    #Ellipse
    if entity.dxftype == 'ELLIPSE':
        EllipseCenter = entity.center
        EllipseMajorAxis = entity.major_axis
        for x in range(len(EllipseCenter)):
            print(" ELLIPSE\n")
            print(EllipseCenter)
            print(EllipseMajorAxis)
    #Arc
    if entity.dxftype == 'ARC':
        ArcCenter = entity.center
        ArcRadius = entity.radius
        ArcStartAngle = entity.start_angle
        ArcEndAngle = entity.end_angle
        for x in range(len(ArcCenter)):
            print(" ARC\n")
            print(ArcCenter)
            print(ArcRadius)
            print(ArcStartAngle)
            print(ArcEndAngle)

    #Polyline
    if entity.dxftype == 'POLYLINE':
        PolylineIsClosed = entity.is_closed
        PolylineSplineType = entity.spline_type
        PolylinePoints = entity.points
        PolylineControlPoints = entity.control_points
        PolylineBulge = entity.bulge
        PolylineVertexCount = entity.__len__()

        print(" POLYLINE\n")
        print("is closed")
        print(PolylineIsClosed)
        print("Spline type")
        print(PolylineSplineType)
        print("line points")
        print(PolylinePoints)
#        print("Control points")
#        print(PolylineControlPoints)
#        print("bulge")
#        print(PolylineBulge)
        print("vertex count:")
        print(PolylineVertexCount)
        PointCnt = PolylineVertexCount
        Bezier = GetBezier(PointCnt, PolylinePoints)
        print("formula")
        print(Bezier)

    #LWPolyline
    if entity.dxftype == 'LWPOLYLINE':
        LWPolylinePoints = entity.points
        LWPolylineIsClosed = entity.is_closed
        LWPolylinePointsCount = entity.__len__()
        for x in range(len(LWPolylineIsClosed)):
            print(" LWPOLYLINE\n")
            print(LWPolylinePoints)
            print(LWPolylineIsClosed)
            print(LWPolylinePointsCount)
    #Spline
    if entity.dxftype == 'SPLINE':
        SplineDegree = entity.degree
        SplineStartTangent = entity.start_tangent
        SplineEndTangent = entity.end_tangent
        SplineControlPoints = entity.control_points
        SplineFitPoints = entity.fit_points
        SplineKnots = entity.knots
        SplineIsClosed = entity.is_closed
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
