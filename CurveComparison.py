import dxfgrabber
import numpy as np
from numpy.linalg import inv
import array as arr
import sympy
from sympy import init_printing, Symbol, UnevaluatedExpr, expand, pretty_print as pprint, latex

#Crete Composite Bezier from given points
def CompositeBezier(n, Points):
    P = PMtrx(n, Points)
    C = [[0],[1],[2],[3]]
    Bezier = []
    t = Symbol('t')
    print((P[0][0] + P[1][0]) / 2)
    MiddlePnt = (P[0][0] + P[1][0]) / 2
    LF = LineFormula(P[0], P[2])
    d = PointDist(P[0][0], P[0][1], P[1][0], P[1][1]) / 3
    dp = DistantPoint(P[1], - d, LineSlope(P[0], P[1]))
    PF = PerpFormula(P[0], P[2], dp)
    print("Middle")
    print(MiddlePnt[0])
    tmp = (LineIntersect(PerpFormula(P[0], P[1], MiddlePnt), LF))
    C[1] = 2 * MiddlePnt - tmp
    C[0].append(P[0])
    C[3].append(P[1])
    for i in (Points):
        LF = LineFormula(P[i + 1], P[i + 2])
        dp = PointDist(P[i], P[i + 1]) / 3
        C[2] = DistantPoint(P[i + 1], - d, LineSlope(P[i], P[i + 1]))
        Bezier.append(BezierFormulaComp(C))
        C[1] = DistantPoint(P[i + 1], d, LineSlope(P[i + 1], P[i + 2]))
        if i == len(Points) - 3:
            MiddlePnt = (P[n - 1][0] + P[n - 2][0]) / 2
            LF = LineFormula(P[i], P[i + 1])
            d = PointDist(P[i], P[i + 2]) / 3
            dp = DistantPoint(P[i + 2], - d, LineSlope(P[i + 1], P[i - 1]))
            PF = PerpFormula(P[i], P[i + 2], dp)
            tmp = (LineIntersect(PerpFormula(P[-1], P[-2], MiddlePnt), LineFormula(P[-1],P[-3])))
            C[1] = 2 * MiddlePnt - tmp
            Bezier.append(BezierFormulaComp(C))
            break
    return Bezier

def BezierFormulaComp(C):
    t = Symbol('t')
    for i in range(4):
        BinC = BinCoef(3 - i, i)
        Bxtmp = (UnevaluatedExpr(BinC))*(1 - t)**(4 - i - 1)*t**(i - 2)*(C[i][0])    #Create Bx(t) formula as a string
        Bx = Bx + Bxtmp
        Bytmp = (UnevaluatedExpr(BinC))*(1 - t)**(4 - i - 1)*t**(i - 2)*(C[i][1])  #Create Bx(t) formula as a string
        By = By + Bytmp
#        Bx = f"{Bx} + ({BinC})(1 - t)^{n - i}*t^{i - 1}({P[i - 1][0]}) " #Create Bx(t) formula as a string
#        By = f"{By} + ({BinC})(1 - t)^{n - i}*t^{i - 1}({P[i - 1][1]}) " #Create Bx(t) formula as a string

    return (Bx, By)

def LineIntersect(P11, P12, P21, P22):
    LF1 = LineFormula(P11, P12)
    LF2 = LineFormula(P21, P22)
    LFx = solve(LF1 - LF2, x)
    Px, Py = solve([LF1, LF2], [x, y])
    P = []
    P.append(x)
    P[x].append(Py)
    return P

#Returns slope of line trough two given points
def LineSlope(P1, P2):
    return (P1[1] - P2[1])/(P1[0] - P2[0])
#Returns distance between two points
def PointDist(ax, ay, bx, by):
    return np.sqrt((ax-bx)**2+(ay-by)**2)

#Point coordinates given distance from point on a line
def DistantPoint(P, d, slope):
    nP = []
    dist = P[0] + d * np.cos(np.arctan(slope))
    nP.append(dist)
    nP.append(P[1] + d * np.sin(np.arctan(slope)))
    return(P)
#Get parametric line Formula
def ParamLineFormula(P1, P2):
    t = Symbol('t')
    PLFx = ((1 - t) * P1[0] + t * P2[0])
    PLFy = ((1 - t) * P1[1] + t * P2[1])

#Get Perpendicular formula
def PerpFormula(P1, P2, P):
    x = Symbol('x')
    print(P[0] * ((P1[0] - P2[0]) / (P1[1] - P2[1])))
    PF = (- (P1[0] + P2[0]) / (P1[1] - P2[1]) * x + P[0] * ((P1[0] - P2[0]) / (P1[1] - P2[1])) + P[1])
    return(PF)

#Get line formula
def LineFormula(P1, P2):
    x = Symbol('x')
    LF = ((P1[1] - P2[1]) / (P1[0] - P2[0])) * x + (P1[0] * P2[1] - P2[0] * P1[1]) / (P1[0] - P2[0])
    return LF

#get factorial
def Fact(n):
    f = 1
    for i in range(n):
        f = f * (i + 1)
    return f
#get binominal coefficients
def BinCoef(n, k):
    return int(Fact(n)/(Fact(k) * Fact(n - k)))
#Getting Bezier (currently from whole polyline)
def GetBezier(n, Pn, Points):
    P = PMtrx(Pn, Points)
    Pn = len(P)
    T = TMtrx(n, Pn, P)
    M = MMtrx(n)
    C = CMtrx(Pn, M, P, T)
    return BezierFormula(n, C)


def BezierFormula(n, C):
    Bx = 0
    By = 0
    t = Symbol('t')
    for i in range(1, n + 1):
        BinC = BinCoef(n - 1, i - 1)
        Bxtmp = (UnevaluatedExpr(BinC))*(1 - t)**(n - i)*t**(i - 1)*(C[i - 1][0])    #Create Bx(t) formula as a string
        Bx = Bx + Bxtmp
        Bytmp = (UnevaluatedExpr(BinC))*(1 - t)**(n - i)*t**(i - 1)*(C[i - 1][1])  #Create Bx(t) formula as a string
        By = By + Bytmp
#        Bx = f"{Bx} + ({BinC})(1 - t)^{n - i}*t^{i - 1}({P[i - 1][0]}) " #Create Bx(t) formula as a string
#        By = f"{By} + ({BinC})(1 - t)^{n - i}*t^{i - 1}({P[i - 1][1]}) " #Create Bx(t) formula as a string

    return (Bx, By)

#def ToCubicBezier()

#Get M matrix, where number represents number of Bezier fit-points https://pomax.github.io/bezierinfo/#curvefitting
#Max 8 - degree Bezier to get close enough results
#def MMtrx(i, eq)
def MMtrx(n):
    t = Symbol('t')
    Mtrx = []
    for i in range(n):
        Mtrx.append([])
        BinC = BinCoef(n - 1 , i)
        Mtmp = expand((BinC)*(1 - t)**(n - i - 1)*t**(i))
        for j in range(n):
            Mtrx[i].append(float(Mtmp.coeff(t, n - j - 1)))
    return(Mtrx)

#Get fit-point matrix
def PMtrx(n, Points):
    it = (len(Points) - 1) / (n - 1)   #choose points evenly from polyline point array
    if (it == 0):
        it = 1
    P = []
    k = 0
    for i in range(n - 1):
        P.append(Points[int(round(k))])
        k += it
        if (int(round(k)) >= len(Points) - 1):
            if (P[-1] != Points[-1]):
                P.append(Points[-1])
            return P
    P.append(Points[-1])
    return P
#get T matrix with parameter values
def TMtrx(n, Pn, P):
    d = []
    d.append(0)
    for i in range(1, Pn):   #create point (t) distance array
        dist = d[i - 1] + PointDist(P[i-1][0], P[i][0], P[i-1][1], P[i][1])
        d.append(dist)
    for i in range(Pn):   #scale points to interval [0..1]
        d[i] = d[i] / d[-1]
    T = []
    for i in range(Pn):   #fill T matrix
        T.append([])
        for j in range(n):
            T[i].append(d[i]**j)
    return T
#get controlpoint matrix
def CMtrx(n, M, P, T):
    M = np.flip(M, 0)
    Tt = np.transpose(T)
    Mi = inv(M)
    C = np.matmul(Tt, T)
    C = inv(C)
    C = np.matmul(Mi, C)
    C = np.matmul(C, Tt)
    C = np.matmul(C, P)
    return C






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
        Bx, By = GetBezier(5, min(8, PointCnt), PolylinePoints)
        print("formula")
#        print(latex(Bx))
#        print(latex(By))
        CompositeBezier(min(8, PointCnt), PolylinePoints)
        for i in (Bezier):
            a = 1
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
