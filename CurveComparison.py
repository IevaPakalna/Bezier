#PF problem - (something with const y line)
import dxfgrabber
import numpy as np
from numpy.linalg import inv
import array as arr
import sympy
from sympy import init_printing, Symbol, UnevaluatedExpr, expand, pretty_print as pprint, latex
from sympy.solvers import solve

#Crete Composite Bezier from given points
def CompositeBezier(n, Points):
    P = PMtrx(len(Points), Points)
    C = [[0],[1],[2],[3]]
    Bezier = []
    t = Symbol('t')
    x = Symbol('x')
    LF = LineFormula(P[1], P[2])
    d = PointDist(P[0], P[1])
    C[2] = DistantPoint(P[1], - d / 3, LineSlope(P[1], P[2]))
    dp = DistantPoint(P[1], -d / 2, LineSlope(P[0], P[1]))
    PF1 = PerpFormula(P[0], P[1], dp)

    PF1tmp = []
    PF1tmp.append(0)
    PF1tmp.append((PF1).subs(x, 0))

    PF2 = PerpFormula(dp, PF1tmp, C[2])

    PF2tmp = []
    PF2tmp.append(0)
    PF2tmp.append((PF2).subs(x, 0))

    PF1xPF2 = LineIntersect(dp, PF1tmp, C[2], PF2tmp)
    C1x = 2 * PF1xPF2[0] - C[2][0]
    C1y = 2 * PF1xPF2[1] - C[2][1]

    C1 = []
    C1.append(C1x)
    C1.append(C1y)
    C[1] = C1

    C[0] = (P[0])
    C[3] = (P[1])
    Bezier.append(BezierFormulaComp(C))
    for i in range(len(P)):
        LF = LineFormula(P[i + 1], P[i + 2])
        dp = PointDist(P[i], P[i + 1]) / 3
        C[2] = DistantPoint(P[i + 1], - dp, LineSlope(P[i + 1], P[i + 2]))
        C[0] = P[i]
        C[3] = P[i + 1]
        Bezier.append(BezierFormulaComp(C))
        C[1] = DistantPoint(P[i], dp, LineSlope(P[i - 1], P[i]))
        if i == len(P) - 3:
            LF = LineFormula(P[-2], P[-3])
            d = PointDist(P[-1], P[-2])
            C[1] = DistantPoint(P[-2], - d / 3, LineSlope(P[-2], P[-3]))
            dp = DistantPoint(P[-2], -d / 2, LineSlope(P[-1], P[-2]))
            PF1 = PerpFormula(P[-1], P[-2], dp)

            PF1tmp = []
            PF1tmp.append(0)
            PF1tmp.append((PF1).subs(x, 0))

            PF2 = PerpFormula(dp, PF1tmp, C[1])

            PF2tmp = []
            PF2tmp.append(0)
            PF2tmp.append((PF2).subs(x, 0))

            PF1xPF2 = LineIntersect(dp, PF1tmp, C[1], PF2tmp)
            C2x = 2 * PF1xPF2[0] - C[1][0]
            C2y = 2 * PF1xPF2[1] - C[1][1]

            C2 = []
            C2.append(C2x)
            C2.append(C2y)
            C[2] = C2

            C[0] = (P[-1])
            C[3] = (P[-4])
            Bezier.append(BezierFormulaComp(C))
            break
    return Bezier

def BezierFormulaComp(C):
    t = Symbol('t')
    Bx = 0
    By = 0
    for i in range(4):
        BinC = BinCoef(3, i)
        Bxtmp = (UnevaluatedExpr(BinC))*(1 - t)**(4 - i - 1)*t**(i)*(C[i][0])    #Create Bx(t) formula
        Bx = UnevaluatedExpr(Bx) + Bxtmp
        Bytmp = (UnevaluatedExpr(BinC))*(1 - t)**(4 - i - 1)*t**(i)*(C[i][1])  #Create Bx(t) formula
        By = UnevaluatedExpr(By) + Bytmp
#        Bx = f"{Bx} + ({BinC})(1 - t)^{n - i}*t^{i - 1}({P[i - 1][0]}) " #Create Bx(t) formula as a string
#        By = f"{By} + ({BinC})(1 - t)^{n - i}*t^{i - 1}({P[i - 1][1]}) " #Create Bx(t) formula as a string
    print(Bx, '  ', By)
    return (Bx, By)

def LineIntersect(P11, P12, P21, P22):
    S1 = LineSlope(P11, P12)
    S2 = LineSlope(P21, P22)
    print(P11[0])
    print(P11[1])
    print(P12[0])
    print(P12[1])
    b1 = ((P11[0] * P12[1] - P12[0] * P11[1]) / (P11[0] - P12[0]))
    b2 = ((P21[0] * P22[1] - P22[0] * P21[1]) / (P21[0] - P22[0]))
    a = np.array([[S1, b1],[S2, b2]], dtype = 'float')
    b = np.array([1, 1], dtype = 'float')
    print(a)
    print(b)
    Ptmp = np.linalg.solve(a,b)
    P = []
    P.append(Ptmp[0])
    P.append(Ptmp[1])
    return P

#Returns slope of line trough two given points
def LineSlope(P1, P2):
    return (P1[1] - P2[1])/(P1[0] - P2[0])
#Returns distance between two points
def PointDist(a, b):
    return np.sqrt((a[0] - b[0])**2+(a[1] - b[1])**2)

#Point coordinates given distance from point on a line
def DistantPoint(P, d, slope):
    nP = []
    nP.append(P[0] + d * np.cos(np.arctan(slope)))
    nP.append(P[1] + d * np.sin(np.arctan(slope)))
    return(nP)
#Get parametric line Formula
def ParamLineFormula(P1, P2):
    t = Symbol('t')
    PLFx = ((1 - t) * P1[0] + t * P2[0])
    PLFy = ((1 - t) * P1[1] + t * P2[1])

#Get Perpendicular formula
def PerpFormula(P1, P2, P):
    x = Symbol('x')
    if P1[1] == P2[1]:
        PF = P[1] + 0 * x
        return PF
    if P1[0] == P2[0]:
        PF = P[0] + 0 * y
        return PF
    PF = - (P1[0] - P2[0]) / (P1[1] - P2[1]) * x + P[0] * ((P1[0] - P2[0]) / (P1[1] - P2[1])) + P[1]
    return(PF)

#Get line formula
def LineFormula(P1, P2):
    x = Symbol('x')
    LF = (UnevaluatedExpr((P1[1] - P2[1]) / (P1[0] - P2[0])) * x + (P1[0] * P2[1] - P2[0] * P1[1]) / (P1[0] - P2[0]))
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
        Bxtmp = (UnevaluatedExpr(BinC))*(1 - t)**(n - i)*t**(i - 1)*(C[i - 1][0])    #Create Bx(t) formula
        Bx = Bx + Bxtmp
        Bytmp = (UnevaluatedExpr(BinC))*(1 - t)**(n - i)*t**(i - 1)*(C[i - 1][1])  #Create Bx(t) formula
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
        dist = d[i - 1] + PointDist(P[i-1], P[i])
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
        Bezier = []
        Bezier.append(CompositeBezier(min(8, PointCnt), PolylinePoints))
        for i in (Bezier):
            print()
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
