#PF problem - (something with const y line)
import dxfgrabber
import numpy as np
from numpy.linalg import inv
import array as arr
import sympy
from sympy import init_printing, Symbol, UnevaluatedExpr, expand, pretty_print as pprint, latex
from sympy.solvers import solve
import matplotlib.pyplot as plt
import matplotlib.axes as axes


plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
major_ticks = np.arange(-10000, 10000, 100)
minor_ticks = np.arange(-10000, 10000, 10)
ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor = True)
ax.set_yticks(major_ticks)
ax.set_yticks(minor_ticks, minor = True)
ax.grid(which = 'both')
ax.grid(which = 'major', alpha = 0.5)
ax.grid(which = 'minor', alpha = 0.2)

#Crete Composite Bezier from given points
def CompositeBezier(Points, nr):
    P = PMtrx(len(Points), Points)
    C = [[0],[1],[2],[3]]
    CP = []
    Bezier = []
    t = Symbol('t')
    x = Symbol('x')
    dtmp = PointDist(P[0], P[1])
    d = dtmp
    C[2] = DistantPoint(P[1], P[2], P[0], d / 3) #Point on (P[1], P[2]) line in 1/3 distance of [P[0],P[1]] (C2 controlpoint)
    middlePnt = DistantPoint(P[0], P[0], P[1], d / 2) #Middle point of [P[0],P[1]]
    PF1 = PerpFormula(P[0], P[1], middlePnt)   #Middle perpendicular of [P[0],P[1]]

    PF1tmp = [] #Point on line PF1
    if P[0][0] == P[1][0]:
        PF1tmp.append(-1)
        PF1tmp.append(middlePnt[1])
    elif P[0][1] == P[1][1]:
        PF1tmp.append(middlePnt[0])
        PF1tmp.append(+1)
    else:
        PF1tmp.append(0)
        PF1tmp.append((PF1).subs(x, 0))

    PF2 = PerpFormula(middlePnt, PF1tmp, C[2]) #Perpendicular of PF1 that goes trough C[2]

    PF2tmp = [] #Point on line PF2
    #cases where x or y - constant
    if middlePnt[1] == PF1tmp[1]:
        PF2tmp.append(C[2][0])
        PF2tmp.append(+1)
    elif middlePnt[0] == PF1tmp[0]:
        PF2tmp.append(-1)
        PF2tmp.append(C[2][1])
    else:
        PF2tmp.append(0)
        PF2tmp.append((PF2).subs(x, 0))
    PF1xPF2 = LineIntersect(middlePnt, PF1tmp, C[2], PF2tmp)   #Intersection point of PF1 and PF2
    #C1 controlpoint coords
    C1x = 2 * PF1xPF2[0] - C[2][0]
    C1y = 2 * PF1xPF2[1] - C[2][1]

    C1 = []
    C1.append(C1x)
    C1.append(C1y)
    C[1] = C1

    C[0] = (P[0])
    C[3] = (P[1])
    CP.append([C[0], C[1], C[2], C[3]])
    Bx, By = BezierFormulaComp(C)
    Bezier.append([Bx, By])
    for i in range(1, len(P)): #Calculate controlpoints for P[i], P[i + 1] segment
        dtmp = PointDist(P[i], P[i + 1]) / 3
        d = dtmp
        C[2] = DistantPoint(P[i + 1], P[i + 2], P[i], d)
        C[0] = P[i]
        C[3] = P[i + 1]
        C[1] = DistantPoint(P[i], P[i - 1], P[i + 1], d)
        CP.append([C[0], C[1], C[2], C[3]])
        Bx, By = BezierFormulaComp(C)
        Bezier.append([Bx, By])
        if i == len(P) - 3:
            dtmp = PointDist(P[-1], P[-2]) / 3
            d = dtmp
            C[1] = DistantPoint(P[-2], P[-3], P[-1], d)
            middlePnt = DistantPoint(P[-1], P[-1], P[-2], d)
            PF1 = PerpFormula(P[-1], P[-2], middlePnt)

            PF1tmp = []
            if P[-1][0] == P[-2][0]:
                PF1tmp.append(-1)
                PF1tmp.append(middlePnt[-2])
            elif P[-1][1] == P[-2][1]:
                PF1tmp.append(middlePnt[-1])
                PF1tmp.append(+1)
            else:
                PF1tmp.append(0)
                PF1tmp.append((PF1).subs(x, 0))

            PF2 = PerpFormula(middlePnt, PF1tmp, C[1])

            PF2tmp = []
            #cases where x or y - constant
            if middlePnt[-2] == PF1tmp[-2]:
                PF2tmp.append(C[-3][0])
                PF2tmp.append(+1)
            elif middlePnt[-1] == PF1tmp[-1]:
                PF2tmp.append(-1)
                PF2tmp.append(C[-3][1])
            else:
                PF2tmp.append(0)
                PF2tmp.append((PF2).subs(x, 0))
            PF1xPF2 = LineIntersect(middlePnt, PF1tmp, C[1], PF2tmp)
            C2x = 2 * PF1xPF2[0] - C[1][0]
            C2y = 2 * PF1xPF2[1] - C[1][1]
            C2 = []
            C2.append(C2x)
            C2.append(C2y)
            C[2] = C2
            C[0] = (P[-1])
            C[3] = (P[-2])
            CP.append([C[0], C[1], C[2], C[3]])
            Bx, By = BezierFormulaComp(C)
            Bezier.append([Bx, By])
            break
    return Bezier, CP

#Coordinates of projected point on a line
def PointProjecOnLine(P, P1, P2):
    slope = LineSlope(P1, P2)
    b = 0
    if P1[0] == P2[0]:
        b = P1[0]
    else:
        b = ((P1[0] * P2[1] - P2[0] * P1[1]) / (P1[0] - P2[0]))

    x = (- (- P[0] - slope * P[1]) - slope * b) / (slope**2 + 1)
    y = (slope * (P[0] + slope * P[1]) + b) / (slope**2 + 1)
    Pp = []
    Pp.append(x)
    Pp.append(y)
    return Pp

def BezierFormulaComp(C):
    t = Symbol('t')
    Bx = 0
    By = 0
    for i in range(4):
        BinC = BinCoef(3, i)
        Bxtmp = BinC*(1 - t)**(4 - i - 1)*t**(i)*(C[i][0])    #Create Bx(t) formula
        Bx = Bx + Bxtmp
        Bytmp = BinC*(1 - t)**(4 - i - 1)*t**(i)*(C[i][1])  #Create Bx(t) formula
        By = By + Bytmp

#        Bx = f"{Bx} + ({BinC})(1 - t)^{n - i}*t^{i - 1}({P[i - 1][0]}) " #Create Bx(t) formula as a string
#        By = f"{By} + ({BinC})(1 - t)^{n - i}*t^{i - 1}({P[i - 1][1]}) " #Create Bx(t) formula as a string
    return Bx, By

def PlotBezier(C, clr, transp, lstyle) :
    for i in range(len(C)) :
        t1 = np.linspace(0, 1, 50)
        Bx1 = 0
        By1 = 0
        for j in range(4) :
            BinC = BinCoef(3, j)
            Bxtmp = BinC * (1 - t1)**(4 - j - 1) * t1**(j) * (C[i][j][0])    #Create Bx(t) formula
            Bx1 = (Bx1) + Bxtmp
            Bytmp = BinC * (1 - t1)**(4 - j - 1) * t1**(j) * (C[i][j][1])  #Create Bx(t) formula
            By1 = (By1) + Bytmp
        plt.plot(Bx1, By1, color = clr, alpha = transp, linestyle = lstyle)
    return


def LineIntersect(P11, P12, P21, P22):
    y1 = - 1
    y2 = - 1
    Ptmp = []
    if P11[0] == P12[0] and P21[1] == P22[1] :
        Ptmp.append(P11[0])
        Ptmp.append(P21[1])
        return Ptmp
    if P11[1] == P12[1] and P21[0] == P22[0] :
        Ptmp.append(P21[0])
        Ptmp.append(P11[1])
        return Ptmp
    if P11[1] == P12[1] :
        S1 = 0
    else:
        S1 = LineSlope(P11, P12)
    if P21[1] == P22[1] :
        S2 = 0
    else:
        S2 = LineSlope(P21, P22)
    if P11[0] == P12[0]:
        y1 = 0
        b1 = P11[0]
    else:
        b1 = ((P11[0] * P12[1] - P12[0] * P11[1]) / (P11[0] - P12[0]))
    if P21[0] == P22[0]:
        y2 = 0
        b2 = P21[0]
    else:
        b2 = ((P21[0] * P22[1] - P22[0] * P21[1]) / (P21[0] - P22[0]))
    a = np.array([[S1, y1],[S2, y2]], dtype = 'float')
    b = np.array([-b1, -b2], dtype = 'float')
    Ptmp = np.linalg.solve(a,b)
    P = []
    P.append(Ptmp[0])
    P.append(Ptmp[1])
    return P

#Returns slope of line trough two given points
def LineSlope(P1, P2):
    if P1[0] == P2[0]:
        return 0
    return (P1[1] - P2[1])/(P1[0] - P2[0])
#Returns distance between two points
def PointDist(a, b):
    dist = np.sqrt(round((a[0] - b[0])**2 + (a[1] - b[1])**2, 5))
    return dist

#Calculate distant point along a line a certain distance away given point
def DistantPoint(P, P1, P2, d):
    Ptmp = []
    dist = 0
    v = []
    P1n = []
    P2n = []
    if P != P1 :
        Ptmp = PointProjecOnLine(P, P1, P2) #Project point on line
        v.append(P[0] - Ptmp[0]) #Get distance vector
        v.append(P[1] - Ptmp[1])
        P1n.append(P1[0] + v[0])
        P1n.append(P1[1] + v[1])
        P2n.append(P2[0] + v[0])
        P2n.append(P2[1] + v[1])
        v1 = P2n[0] - P[0]
        v2 = P2n[1] - P[1]
    else:
        v1 = P2[0] - P[0]
        v2 = P2[1] - P[1]
    u1 = float(v1 / np.sqrt(v1**2 + v2**2))
    u2 = float(v2 / np.sqrt(v1**2 + v2**2))
    nP = []
    nP.append(P[0] + d * u1)
    nP.append(P[1] + d * u2)
    return(nP)

#Get parametric line Formula
def ParamLineFormula(P1, P2):
    t = Symbol('t')
    PLFx = ((1 - t) * P1[0] + t * P2[0])
    PLFy = ((1 - t) * P1[1] + t * P2[1])
    return PLFx, PLFy

#Get Perpendicular formula
def PerpFormula(P1, P2, P):
    x = Symbol('x')
    if P1[1] == P2[1]:
        PF = P[1] + 0 * x
        return PF
    PF = - (P1[0] - P2[0]) / (P1[1] - P2[1]) * x + P[0] * ((P1[0] - P2[0]) / (P1[1] - P2[1])) + P[1]
    return(PF)

#Get line formula
def LineFormula(P1, P2):
    x = Symbol('x')
    if P1[1] == P2[1]:
        LF = (P1[1] + 0 * x)
        return LF
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

#Binary sort (for abscissa)
def SortInsertPos(P, points, l, r): #l - left side, r - right side of segment in array
    if l == r :
        if points[l][0] == P[0] :
            if points[l][1] == P[1] :
                return -1
            if points[l][1] < P[1] :
                return l + 1
            return l
        if points[l][0] < P[0] :
            return l + 1
        return l
    med = l + (r - l) // 2
    if points[med][0] > P[0] :
        return SortInsertPos(P, points, l, med)
    if points[med][0] < P[0] :
        return SortInsertPos(P, points, med + 1, r)
    if points[med][0] == P[0] :
        if points[med][1] == P[1] :
            return -1 #there already exists identical point, therefore we will not save it
        if points[med][1] > P[1] :
            return SortInsertPos(P, points, l, med)
        if points[med][1] < P[1] :
            return SortInsertPos(P, points, med + 1, r)

#transform points
def pointTransform(P, Vx, Vy, rP, alpha) :
    for i in range(len(P) - 1) :
        tmpx = (P[i][0] - Vx)
        tmpy = (P[i][1] - Vy)
        tmpx = np.cos(alpha) * (tmpx - rP[0]) + np.sin(alpha) * (tmpy - rP[1]) + rP[0]
        tmpy = np.cos(alpha) * (tmpx - rP[0]) + np.sin(alpha) * (tmpy - rP[1]) + rP[1]
        P[i][0] = tmpx
        P[i][1] = tmpy
    return P

dxf1 = dxfgrabber.readfile("parastie_platgurnu_m2_p2_002.dxf")

#type of objects in file
type1 = [entity.dxftype for entity in dxf1.entities]

output1 = [entity for entity in dxf1.entities]

Color1 = [entity.color for entity in output1]
Linetype1 = [entity.linetype for entity in output1 if entity.dxftype != 'POINT']

CP1 = []

points1 = []
line1 = []
Bezier1 = []

#get parameters of objects
for entity in output1:
    #Point
    if entity.dxftype == 'POINT':
        point1 = entity.point
        if len(points1) == 0 :
            points1.append([point1[0], point1[1]])
            continue
        pos = SortInsertPos(point1, points1, 0, len(points1) - 1)
        if pos != -1 :
            points1.insert(pos, [point1[0], point1[1]])
            x = [point1[0]]
            y = [point1[1]]
            plt.plot(x, y, 'o', color = '#055583', alpha = 0.65)
    #Line
    if entity.dxftype == 'LINE':
        lineStart1 = entity.start
        lineEnd1 = entity.end
        line1.append([lineStart1, lineEnd1])
        x = [lineStart1[0], lineEnd1[0]]
        y = [lineStart1[1], lineEnd1[1]]
        plt.plot(x, y, color = '#055583', alpha = 0.65)

    #Circle
    if entity.dxftype == 'CIRCLE':
        CenterPoints1 = entity.center
        Radius1 = entity.radius
        for x in range(len(CenterPoints1)):
            print(" CIRCLE \n")
            print("Center Point")
            print(CenterPoints1[x])
            print("Radius")
            print(Radius1[x])
    #Ellipse
    if entity.dxftype == 'ELLIPSE':
        EllipseCenter1 = entity.center
        EllipseMajorAxis1 = entity.major_axis
        for x in range(len(EllipseCenter1)):
            print(" ELLIPSE\n")
            print(EllipseCenter1)
            print(EllipseMajorAxis1)
    #Arc
    if entity.dxftype == 'ARC':
        ArcCenter1 = entity.center
        ArcRadius1 = entity.radius
        ArcStartAngle1 = entity.start_angle
        ArcEndAngle1 = entity.end_angle
        for x in range(len(ArcCenter1)):
            print(" ARC\n")
            print(ArcCenter1)
            print(ArcRadius1)
            print(ArcStartAngle1)
            print(ArcEndAngle1)

    #Polyline
    if entity.dxftype == 'POLYLINE':
        PolylineIsClosed1 = entity.is_closed
        PolylineSplineType1 = entity.spline_type
        PolylinePoints1 = entity.points
        PolylineControlPoints1 = entity.control_points
        PolylineBulge1 = entity.bulge
        PolylineVertexCount1 = entity.__len__()
        PointCnt1 = PolylineVertexCount1
        Beziertmp, CP = CompositeBezier(PolylinePoints1, 1)
        Bezier1.append(Beziertmp)
        CP1.append(CP)

    #LWPolyline
    if entity.dxftype == 'LWPOLYLINE':
        LWPolylinePoints1 = entity.points
        LWPolylineIsClosed1 = entity.is_closed
        LWPolylinePointsCount1 = entity.__len__()
        for x in range(len(LWPolylineIsClosed1)):
            print(" LWPOLYLINE\n")
            print(LWPolylinePoints1)
            print(LWPolylineIsClosed1)
            print(LWPolylinePointsCount1)
    #Spline
    if entity.dxftype == 'SPLINE':
        SplineDegree1 = entity.degree
        SplineStartTangent1 = entity.start_tangent
        SplineEndTangent1 = entity.end_tangent
        SplineControlPoints1 = entity.control_points
        SplineFitPoints1 = entity.fit_points
        SplineKnots1 = entity.knots
        SplineIsClosed1 = entity.is_closed
        for x in range(len(SplineIsClosed1)):
            print(" SPLINE\n")
            print(SplineDegree1)
            print(SplineStartTangent1)
            print(SplineEndTangent1)
            print(SplineControlPoints1)
            print(SplineFitPoints1)
            print(SplineKnots1)
            print(SplineIsClosed1)

#Block
#print("BLOCKS")
#BlockBasepoint= [Block.BlockBasepoint for Block in dxf.entities]
#print(BlockBasepoint)

#Bezier formulas of file #2


dxf2 = dxfgrabber.readfile("Parastie_platgurnu_m3_p2_002.dxf")
#type of objects in file
type2 = [entity.dxftype for entity in dxf2.entities]

output2 = [entity for entity in dxf2.entities]

Color2 = [entity.color for entity in output2]
Linetype2 = [entity.linetype for entity in output2 if entity.dxftype != 'POINT']

CP2 = []
points2 = []
line2 = []
Bezier2 = []
#get parameters of objects
for entity in output2:
    #Point
    if entity.dxftype == 'POINT':
        point2 = entity.point
        if len(points2) == 0 :
            points2.append([point2[0], point2[1]])
            continue
        pos = SortInsertPos(point2, points2, 0, len(points2) - 1)
        if pos != -1 :
            x = [point2[0]]
            y = [point2[1]]
            points2.insert(pos, [point2[0], point2[1]])
            plt.plot(x, y, 'o', color = '#bc0e13', alpha = 0.5)
    #Line
    if entity.dxftype == 'LINE':
        lineStart2 = entity.start
        lineEnd2 = entity.end
        line2.append([lineStart2, lineEnd2])
        x = [lineStart2[0], lineEnd2[0]]
        y = [lineStart2[1], lineEnd2[1]]
        plt.plot(x, y, color = '#bc0e13', alpha = 0.5)

    #Circle
    if entity.dxftype == 'CIRCLE':
        CenterPoints2 = entity.center
        Radius2 = entity.radius
        for x in range(len(CenterPoints2)):
            print(" CIRCLE \n")
            print("Center Point")
            print(CenterPoints2[x])
            print("Radius")
            print(Radius2[x])
    #Ellipse
    if entity.dxftype == 'ELLIPSE':
        EllipseCenter2 = entity.center
        EllipseMajorAxis2 = entity.major_axis
        for x in range(len(EllipseCenter2)):
            print(" ELLIPSE\n")
            print(EllipseCenter2)
            print(EllipseMajorAxis2)
    #Arc
    if entity.dxftype == 'ARC':
        ArcCenter2 = entity.center
        ArcRadius2 = entity.radius
        ArcStartAngle2 = entity.start_angle
        ArcEndAngle2 = entity.end_angle
        for x in range(len(ArcCenter2)):
            print(" ARC\n")
            print(ArcCenter2)
            print(ArcRadius2)
            print(ArcStartAngle2)
            print(ArcEndAngle2)

    #Polyline
    if entity.dxftype == 'POLYLINE':
        PolylineIsClosed2 = entity.is_closed
        PolylineSplineType2 = entity.spline_type
        PolylinePoints2 = entity.points
        PolylineControlPoints2 = entity.control_points
        PolylineBulge2 = entity.bulge
        PolylineVertexCount2 = entity.__len__()
        PointCnt2 = PolylineVertexCount2
        #Bx2, By2 = GetBezier(5, min(8, PointCnt2), PolylinePoints2)
#        print(latex(Bx))
#        print(latex(By))
        Beziertmp, CP = CompositeBezier(PolylinePoints2, 2)
        Bezier2.append(Beziertmp)
        CP2.append(CP)

    #LWPolyline
    if entity.dxftype == 'LWPOLYLINE':
        LWPolylinePoints2 = entity.points
        LWPolylineIsClosed2 = entity.is_closed
        LWPolylinePointsCount2 = entity.__len__()
        for x in range(len(LWPolylineIsClosed2)):
            print(" LWPOLYLINE\n")
            print(LWPolylinePoints2)
            print(LWPolylineIsClosed2)
            print(LWPolylinePointsCount2)
    #Spline
    if entity.dxftype == 'SPLINE':
        SplineDegree2 = entity.degree
        SplineStartTangent2 = entity.start_tangent
        SplineEndTangent2 = entity.end_tangent
        SplineControlPoints2 = entity.control_points
        SplineFitPoints2 = entity.fit_points
        SplineKnots2 = entity.knots
        SplineIsClosed2 = entity.is_closed
        for x in range(len(SplineIsClosed2)):
            print(" SPLINE\n")
            print(SplineDegree2)
            print(SplineStartTangent2)
            print(SplineEndTangent2)
            print(SplineControlPoints2)
            print(SplineFitPoints2)
            print(SplineKnots2)
            print(SplineIsClosed2)


for i in CP1 :
    PlotBezier(i, '#055583', 0.65, 'solid')
for i in CP2 :
    PlotBezier(i, '#bc0e13', 0.5, 'solid')

#We have got data from both given files

#File #1
#for i in points1:
#    print('(', i[0], ',', i[1], ')')
#for i in line1:
#    print('(',latex(i[0]),',', latex(i[1]),')')
#for i in Bezier1 :
#    for j in range(len(i)):
#    #    print('(',latex(i[j][0]),',', latex(i[j][1]),')')
#        a = 1




#File #2
#for i in points2:
#    print('(', i[0], ',', i[1], ')')
#for i in line2:
#    print('(',latex(i[0]),',', latex(i[1]),')')
#for i in Bezier2 :
#    for j in range(len(i)):
#        #    print('(',latex(i[j][0]),',', latex(i[j][1]),')')
#        a = 1



#Lets find which points are represented as the 'same' points in both files

#First we will take 3 points from first file and get the distances between them
#Based on those 3 point we are going to take three point sets from second file
#   and compare if they match based on proportions
#(from those three points we save proportion and angle and in the next steps we
#    are going to use vectors to find the corresponding points)
#If we have found the match then - take one point, pair it with all the other
#   points and find their matching pairs in second file
#Because we do that with vector there is no reason to compare each two points


dist11 = PointDist(points1[0], points1[1])
dist12 = PointDist(points1[0], points1[2])
dist13 = PointDist(points1[1], points1[2])
for i in range(len(points2) - 1) :
    t = True
    transf = False
    dist21 = PointDist(points2[i], points2[i + 1])
    for j in range(len(points2)) :
        if j == i:
            continue
        dist22 = PointDist(points2[i], points2[j])
        if dist11 / dist21 == dist12/ dist22 :
            dist23 = PointDist(points2[i + 1], points2[j])
            if dist11 / dist21 == dist13 / dist23 :
                ratio = dist21 / dist11
                dist = PointDist(points1[0], points2[i])
                unitV = []  #Unit vector of difference between points1 and points2
                tmpVx = points1[0][0] - points2[i][0]
                tmpVy = points2[0][1] - points2[i][1]
                if dist != 0 :
                    tmpVx /= abs(dist)
                    tmpVy /= abs(dist)
                unitV.append(tmpVx)
                unitV.append(tmpVy)
                for k in range(3, len(points1)) :
                    dist1 = PointDist(points1[0], points1[k])
                    tmpPx = points2[i][0] + unitV[0] * dist1
                    tmpPy = points2[i][1] + unitV[1] * dist1
                    tmpP = []
                    tmpP.append(tmpPx)
                    tmpP.append(tmpPy)
                    try:
                        points2.index(tmpP)
                    except ValueError :
                        t = False
                        break
                    if k == len(points2) - 1 :
                        d = PointDist(points1[0], points2[i])
                        transfVx = points2[i][0] - points1[0][0] #Transformation vector
                        transfVy = points2[i][1] - points1[0][0]
                        Ptmpx = points2[j + 1][0] - transfVx
                        Ptmpy = points2[j + 1][1] - transfVy
                        if (((points1[0][0] - points2[i][0]) == 0) and ((points1[0][0] - Ptmpx) == 0)) :
                            alpha = 0
                        elif (points1[0][0] - points2[i][0]) == 0 :
                            a2 = (points1[0][1] - Ptmpy) / (points1[0][0] - Ptmpx)
                            alpha = 90 - np.arctan(a2)
                        elif (points1[0][0] - Ptmpx[i][0]) == 0 :
                            a1 = (points1[0][1] - points2[i][1]) / (points1[0][0] - points2[i][0])
                            alpha = 90 - np.arcatn(a1)
                        else:
                            a1 = (points1[0][1] - points2[i][1]) / (points1[0][0] - points2[i][0])
                            a2 = (points1[0][1] - Ptmpy) / (points1[0][0] - Ptmpx)
                            alpha = np.arctan((a2 - a1) / (1 + a1 * a2))
                        points2 = pointTransform(points2, transfVx, transfVy, points1[0], - alpha)
                        transf = True
                        break
                if t == False or transf == True :
                    break
            if t == False or transf == True :
                break
        if transf == True :
            break
    if transf == True :
        break



#Comparison of objects
#we will go through all objects in both files by picking the first object in first files
#then in second file we will find the most similar one (by type and coordinates)
#we also know that some of the objects are connected which we are going to use to specify two 'equal' objects in both files
#to determine how similar are both files we will calculate max distance between respective objects where:
#   *points - distance between points (but this will not be calculated because for our target it`s not important)
#   lines - distance between end-points
#   Bezier curves - distance for set of points on Bezier (with t parameter change of 0.1)
#       for Bezier we need to iterate both lines because one of them can end before the other one and actual max length can be larger

#for second type of comparison we will compare respective Bezier curves by getting distance between points with equal parameter value


#find closest point in lines array
def ClosestPntLine(P, arr):
    min = -1
    minPind = 0
    for i in range(len(arr) - 1) :
        if i == 0 :
            min = PointDist(P, arr[i][0])
            minPind = i
            lineSt = 0
        tmpD = PointDist(P, arr[i][0])
        if tmpD < min :
            min = tmpD
            minPind = i
            lineSt = 0
        tmpD = PointDist(P, arr[i][1])
        if tmpD < min :
            min = tmpD
            minPind = i
            lineSt = 1
    return arr[minPind][lineSt], min

#find closest point in Bezier curves array
def ClosestPntBezier(P, arr) :
    min = -1
    minPind = 0
    t = Symbol('t')
    for i in range(len(arr) - 1) :
        #as the Bezier curve is made from multiple cubic Bezier curves, we have to compare only the start point of first and end point of last cubic curve that belongs to specific composite Bezier
        if i == 0 :
            min = PointDist(P, (arr[i][0][0].subs(t, 0), arr[i][0][1].subs(t, 0)))
            minPind = i
            lineSt = 0
        tmpD = PointDist(P, (arr[i][0][0].subs(t,0), arr[i][0][1].subs(t,0)))
        if tmpD < min :
            min = tmpD
            minPind = i
            lineSt = 0
        tmpD = PointDist(P, (arr[i][-1][0].subs(t, 1), arr[i][-1][1].subs(t, 1)))
        if tmpD < min :
            min = tmpD
            minPind = i
            lineSt = 1
    return (arr[minPind][-lineSt][0].subs(t, lineSt), arr[minPind][-lineSt][1].subs(t, lineSt)), min

def OneFileObject(P, line, Bezier, CP, visited, fnr) :
    for i in line :
        if (P == i[0] or P == i[1]) :
            if fnr == 1 :
                plt.plot((i[0][0], i[0][1]), (i[1][0], i[1][1]), color = '#055583', linestyle = 'dotted')
            if fnr == 2 :
                plt.plot((i[0][0], i[0][1]), (i[1][0], i[1][1]), color = '#bc0e13', linestyle = 'dotted')

    t = Symbol('t')
    lenB = len(Bezier)
    for i in range(lenB) :
        if (P == (Bezier[0][0].subs(t, 0), Bezier[0][1].subs(t, 0)) or (P == Bezier[-1][0].subs(t, 1), Bezier[-1][1].subs(t, 1))) :
            if fnr == 1 :
                PlotBezier(CP[i], '#055583', 1, 'dotted')
            if fnr == 2 :
                PlotBezier(CP[i], '#bc0e13', 1, 'dotted')

    visited[P] = 2
    return visited

def CubicBezierLen(C):
    t = Symbol('t')
    ax = - 3 * C[0][0] + 9 * C[1][0] - 9 * C[2][0] + 3 * C[3][0]
    bx = 6 * C[0][0] - 12 * C[1][0] + 6 * C[2][0]
    cx = - 3 * C[0][0] + 3 * C[1][0]

    ay = - 3 * C[0][1] + 9 * C[1][1] - 9 * C[2][1] + 3 * C[3][1]
    by = 6 * C[0][1] - 12 * C[1][1] + 6 * C[2][1]
    cy = - 3 * C[0][1] + 3 * C[1][1]

    sqtmp = (ax**2 + ay**2) * t**4 + 2 * (ax * bx + ay * by) * t**3 + (2 * ax * cx + bx**2 + 2 * ay * cy + by**2) * t**2 + 2 * (bx * cx + by * cy) * t + (cx**2 + cy**2)
    # from https://pomax.github.io/bezierinfo/legendre-gauss.html
    w = []
    w.append(0.8888888888888888)
    w.append(0.5555555555555556)
    w.append(0.5555555555555556)

    x = []
    x.append(0.0000000000000000)
    x.append(-0.7745966692414834)
    x.append(0.7745966692414834)

    L = 0
    for i in range(3) :
        fttmp = round(sqtmp.subs(t, 0.5 * x[i] + 0.5), 10)
        ft = np.sqrt(fttmp)
        L += 0.5 * (w[i] * ft)
    return L

def DistantPointOnBezier(dist, nr, B, lineSt, C) :
    t = Symbol('t')
    xtmp = B[nr][0].subs(t, lineSt)
    ytmp = B[nr][1].subs(t, lineSt)
    i = 0
    param = lineSt
    cnt = 0.1
    m = 0.1
    disttmp = 0
    dist1 = 0
    lenB = len(B)
    while (nr >= 0) and (nr <= lenB - 1) :
        while (param >= 0 and param <= 1) :
            param = abs(lineSt - cnt)
            if (param <= 0 or param >= 1) :
                break;
            x2tmp = B[nr][0].subs(t, param)
            y2tmp = B[nr][1].subs(t, param)
            cnt += m
            Pdist = PointDist((xtmp, ytmp), (x2tmp, y2tmp))
            disttmp += Pdist
            if (disttmp < dist + 0.1 and disttmp > dist - 0.1):
                return x2tmp, y2tmp, nr
            elif disttmp > dist + 0.1 :
                cnt -= 0.1
                m *= 0.1
                disttmp -=Pdist
                continue
            xtmp = x2tmp
            ytmp = y2tmp
        dist1 = dist1 + CubicBezierLen(C[nr])
        disttmp = dist1
        if (nr == lenB - 1 and lineSt == 0) :
            break
        if (nr == 0 and lineSt == 1) :
            break
        if lineSt == 0 :
            nr += 1
            param = lineSt
            m = 0.1
        else:
            nr -= 1
            param = lineSt
            m = 0.1
    x2tmp = B[nr][0].subs(t, param)
    y2tmp = B[nr][1].subs(t, param)
    return x2tmp, y2tmp, nr

#calculate difference value using square method
def BezierDiff2(B1, int1, lineSt1, B2, int2, lineSt2, CP1, CP2) :
    dist1 = 0
    dist2 = 0
    value = 0
    cnt1 = 0
    cnt2 = 0
    maxdist = 0
    P1 = []
    P1.append([])
    P1.append([])
    len1 = len(B1)
    P2 = []
    P2.append([])
    P2.append([])
    len2 = len(B2)
    minP1 = []
    minP1.append([])
    minP1.append([])
    minP2 = []
    minP2.append([])
    minP2.append([])
    minPtmp = []
    minPtmp.append([])
    minPtmp.append([])
    t = Symbol('t')
    for i in range(11) :
        if lineSt1 == 1 :
            nr1 = len(B1) - 1
        else :
            nr1 = 0
        print((dist1, nr1, lineSt1, 1))
        P1x, P1y, nr1 = DistantPointOnBezier(dist1, nr1, B1, lineSt1, CP1)

        print(P1x, P1y, nr1)
        P1[0] = P1x
        P1[1] = P1y
        if lineSt2 == 1 :
            nr2 = len(B2) - 1
        else :
            nr2 = 0
        P2x, P2y, nr2 = DistantPointOnBezier(dist2, nr2, B2, lineSt2, CP2)
        P2[0] = P2x
        P2[1] = P2y
        value += (PointDist((P1x, P1y), (P2x, P2y)))**2
        dist1 += int1
        dist2 += int2
        plt.plot((P1[0], P2[0]), (P1[1], P2[1]), color = '#af5ba3')
        minDist = 10000
        if cnt1 < len2 - 1 :
            minP1 = (B2[cnt1][0].subs(t, 0), B2[cnt1][1].subs(t, 0))
            cnt1, minP1, minDisttmp1 = MinDistTanBezier(cnt1, B2, P1, minDist, len2, minPtmp, minP1)
            plt.plot((P1[0], minP1[0]), (P1[1], minP1[1]), color = '#6c9f92', linewidth = 3)
        if cnt2 < len1 - 1 :
            minP2 = (B1[cnt2][0].subs(t, 0), B1[cnt2][1].subs(t, 0))
            cnt2, minP2, minDisttmp2 = MinDistTanBezier(cnt2, B1, P2, minDist, len1, minPtmp, minP2)
            plt.plot((P2[0], minP2[0]), (P2[1], minP2[1]), color = '#6c9f92')
    return value

def LeastSquare(B1, int1, lineSt1, B2, int2, lineSt2, CP1, CP2) :
    dist1 = 0
    dist2 = 0
    value = 0

    if lineSt1 == 1 :
        nr1 = len(B1) - 1
    else:
        nr1 = 0

    len1 = len(B1)
    len2 = len(B2)

    if lineSt2 == 1 :
        nr2 = len(B2) - 1
    else:
        nr2 = 0

    for i in range(11) :
        P1x, P1y, nr1 = DistantPointOnBezier(dist1, nr1, B1, lineSt1, CP1)
        P2x, P2y, nr2 = DistantPointOnBezier(dist2, nr2, B2, lineSt2, CP2)

        value += (PointDist((P1x, P1y), (P2x, P2y)))**2
        dist1 += int1
        dist2 += int2
        plt.plot([P1x, P2x], [P1y, P2y], color = '#af5ba3')
    return value

def DiffAll(P1, visited, line1, line2, Bezier1, Bezier2, CP1, CP2) :

    if P1 in visited :
        if visited[P1] == 2 :
            return visited

    P21, min21 = ClosestPntLine(P1, line2)
    P22, min22 = ClosestPntBezier(P1, Bezier2)
    if (min21 < min22) :
        P2 = P21
    else:
        P2 = P22

    if (PointDist(P1, P2) > 25 or ((min21 == -1 and min22 == -1))) :
        visited = OneFileObject(P1, line1, Bezier1, CP1, visited, 1)
        return visited

    lenLine1 = len(line1)
    lenLine2 = len(line2)
    for i in range(lenLine1) :
        P1Ind = -1
        lineSt1 = -1
        if P1 == line1[i][0] :
            P1Ind = i
            lineSt1 = 0
        elif P1 == line1[i][1] :
            P1Ind = i
            lineSt1 = 1
        else:
            continue
        a1 = LineSlope(line1[P1Ind][0], line1[P1Ind][1])

        lineSt2 = -1
        P2Ind = -1
        minEndDist = -1 #to find which line is the closest one to the line in line1
        for j in range(lenLine2):
            if (P2 == line2[j][0] and not line2[j][0] in visited):
                lineSt2 = 0
            elif (P2 == line2[j][1] and not P2 == line2[j][1] in visited) :
                lineSt2 = 1
            else:
                continue

            a2 = LineSlope(line2[j][0], line2[j][1])
            angle = np.arctan((a2 - a1) / (1 + abs(a1 * a2)))
            #if the angle between these two lines is bigger than 5, then we eill assume that it is not the searched for line
            if abs(angle) > 5 :
                continue

            if (lineSt1 != -1 and lineSt2 != -1) :
                endDist = PointDist(line1[P1Ind][(lineSt1 + 1) % 2], line2[j][(lineSt2 + 1) % 2])
                if (minEndDist == -1 and endDist < 25):
                    minEndDist = endDist
                    P2Ind = j
                if endDist < minEndDist :
                    minEndDist = endDist
                    P2Ind = j

        if P2Ind == -1 :
#            print(line1[i][0], line1[i][1])
            plt.plot([line1[P1Ind][0][0], line1[P1Ind][1][0]], [line1[P1Ind][0][1], line1[P1Ind][1][1]], color = '#055583', linestyle = 'dotted')
            visited[(line1[P1Ind][0], line1[P1Ind][1])] = 2

        if P2Ind != -1 :
            dist1 = PointDist(line1[P1Ind][lineSt1], line2[P2Ind][lineSt2])
            dist2 = PointDist(line1[P1Ind][(lineSt1 + 1) % 2], line2[P2Ind][(lineSt2 + 1) % 2])
            plt.plot([line1[P1Ind][lineSt1][0], line2[P2Ind][lineSt2][0]], [line1[P1Ind][lineSt1][1], line2[P2Ind][lineSt2][1]], color = '#6c9f92', linewidth = 10)
            plt.plot([line1[P1Ind][(lineSt1 + 1) % 2][0],line2[P2Ind][(lineSt2 + 1) % 2][0]], [line1[P1Ind][(lineSt1 + 1) % 2][1], line2[P2Ind][(lineSt2 + 1) % 2][1]], color = '#6c9f92')

    lenB1 = len(Bezier1)
    lenB2 = len(Bezier2)
    for i in range(lenB1) :
        P1Ind = -1
        lineSt1 = -1
        t = Symbol('t')
        if P1 == (Bezier1[i][0][0].subs(t, 0), Bezier1[i][0][1].subs(t, 0)) :
            P1Ind = i
            lineSt1 = 0
        elif P1 == (Bezier1[i][-1][0].subs(t, 1), Bezier1[i][-1][1].subs(t, 1)) :
            P1Ind = i
            lineSt1 = 1
        else:
            continue
        a1 = LineSlope((Bezier1[i][0][0].subs(t, 0), Bezier1[i][0][1].subs(t, 0)), (Bezier1[i][-1][0].subs(t, 1), Bezier1[i][-1][1].subs(t, 1)))

        lineSt2 = -1
        P2Ind = -1
        minEndDist = -1
        for j in range(lenB2):
            if (P2 == (Bezier2[j][0][0].subs(t, 0), Bezier2[j][0][1].subs(t, 0)) and not (Bezier2[j][0][0].subs(t, 0), Bezier2[j][0][1].subs(t, 0)) in visited) :
                lineSt2 = 0
            elif (P2 == (Bezier2[j][-1][0].subs(t, 1), Bezier2[j][-1][1].subs(t, 1)) and not (Bezier2[j][-1][0].subs(t, 1), Bezier2[j][-1][1].subs(t, 1)) in visited) :
                lineSt2 = 1
            else:
                continue

            a2 = LineSlope((Bezier2[j][0][0].subs(t, 0), Bezier2[j][0][1].subs(t, 0)), (Bezier2[j][-1][0].subs(t, 1), Bezier2[j][-1][1].subs(t, 1)))
            angtmp = round((a2 - a1) / (1 + abs(a1 * a2)), 5)
            angle = np.arctan(angtmp)
            if abs(angle) > 5 :
                continue

            if (lineSt1 != -1 and lineSt2 != -1) :
                endDist = PointDist((Bezier1[i][-(lineSt1 + 1) % 2][0].subs(t, (lineSt1 + 1) % 2), Bezier1[i][-(lineSt1 + 1) % 2][1].subs(t, (lineSt1 + 1) % 2)), (Bezier2[j][-(lineSt2 + 1) % 2][0].subs(t, (lineSt2 + 1) % 2),Bezier2[j][-(lineSt2 + 1) % 2][1].subs(t, (lineSt2 + 1) % 2)))
                if (minEndDist == -1 and endDist < 25):
                    minEndDist = endDist
                    P2Ind = j
                if endDist < minEndDist :
                    minEndDist = endDist
                    P2Ind = j

        if (P2Ind == -1) :
            PlotBezier(CP1[i], '#055583', 1, 'dotted')
            t = Symbol('t')
            visited[(Bezier1[i][-lineSt1][0].subs(t, lineSt1), Bezier1[i][-lineSt1][1].subs(t, lineSt1))] = 2

        if P2Ind != -1 :
            len1 = 0
            lenB1k = len(Bezier1[i])
            for k in range(lenB1k) :
                len1 += CubicBezierLen(CP1[i][k])

            len2 = 0
            lenB2k = len(Bezier2[P2Ind])
            for k in range(lenB2k) :
                len2 += CubicBezierLen(CP2[P2Ind][k])

            int1 = len1 / 10
            int2 = len2 / 10

            LeastSquare(Bezier1[i], int1, lineSt1, Bezier2[P2Ind], int2, lineSt2, CP1, CP2)

#            BezierDiff(Bezier1[i], int1, lineSt1, Bezier2[P2Ind], int2, lineSt2, CP1, CP2)


            if not (Bezier1[i][0][0].subs(t, 0), Bezier1[i][0][1].subs(t, 0)) in visited :
                visited = DiffAll((Bezier1[i][0][0].subs(t, 0), Bezier1[i][0][1].subs(t, 0)), visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
            if not (Bezier2[P2Ind][-1][0].subs(t, 1), Bezier2[P2Ind][-1][1].subs(t, 1)) in visited :
                visited = DiffAll((Bezier2[P2Ind][-1][0].subs(t, 1), Bezier2[P2Ind][-1][1].subs(t, 1) ), visited, line1, line2, Bezier1, Bezier2, CP1, CP2)

        for i in line2 :
            if (P1Ind == -1 and P2 == i[0]) :
                plt.plot([i[0][0], i[1][0]], [i[0][1],i[1][1]], color = '#bc0e13', linewidth = 3)
                if not i[1] in visited :
                    visited = DiffAll(i[1], visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
            if (P1Ind == -1 and P2 == i[1]) :
                plt.plot([i[0][0], i[1][0]], [i[0][1], i[1][1]], color = '#bc0e13', linewidth = 3)
                if not i[0] in visited :
                    visited = DiffAll(i[0], visited, line1, line2, Bezier1, Bezier2, CP1, CP2)

        k = 0
        for i in Bezier2 :
            if (P1 == (i[0][0].subs(t, 0), i[0][1].subs(t, 0)) and P1Ind == -1) :
                PlotBezier(CP2[k], '#bc0e13', 1, 'dotted')
                if not (i[-1][0].subs(t, 1), i[-1][1].subs(t, 1)) in visited :
                    visited = DiffAll((i[1][0].subs(t, 1), i[1][1].subs(t, 1)), visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
            if (P1 == (i[-1][0].subs(t, 1), i[-1][1].subs(t, 1)) and P1Ind == -1) :
                PlotBezier(CP2[k], '#bc0e13', 1, 'dotted')
                if not (i[0][0].subs(t, 0), i[0][1].subs(t, 0)) in visited :
                    visited = DiffAll((i[0][0].subs(t, 0), i[0][1].subs(t, 0)), visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
            k += 1

    visited[P1] = 2
    return visited

def MaxBezierDist(B1, B2) :
    for i in range(len(B2)) :
        cnt = 0
        param = 0
        P = []
        P.append([])
        P.append([])
        while param <= 1 :
            param += cnt
            cnt += 0.1
            x2 = B2[0].subs(t, param)
            y2 = B2[1].subs(t, param)
            dist, P = BezierMinDist(B1[0], B1[1], (x2, y2))
            if dist > max :
                P1[0] = P[0]
                P1[1] = P[1]
                P2[0] = x2
                P2[1] = y2
                max = dist
    if max <= 1.5 :
        return
    plt.plot([P1[0], P2[0]], [P1[1], P2[1]], color = '#667281')
    return

#Search for distance between point and Bezier curve
def BezierMinDist(Bx, By, P) :
    print(P)
    print(Bx)
    t = Symbol('t')
    param = 0
    #get distance values from points on bezier that differs from each other by 0.1 parameter value
    for i in range (10) :
        param = param + 0.1
        x = Bx.subs(t, param)
        y = By.subs(t, param)
        dist = PointDist((x, y), P)
        if i == 0 :
            min = PointDist((x, y), P)
            minParam = round(param, 5)
        elif dist < min :
            min = dist
            minParam = round(param, 5)

    #check if in segments between previously determined points there is not point with smaller distance to the given point
    param = minParam
    dist1 = 0
    dist2 = 0
    x = Bx.subs(t, param)
    y = By.subs(t, param)
    minP = []
    minP.append(x)
    minP.append(y)
    while (param >= minParam - 0.1 and param <= minParam + 0.1 and param <= 1 and param >= 0) :
        if (param <= 0.99 and param >= 0.01):
            x1 = Bx.subs(t, param + 0.01)
            y1 = By.subs(t, param + 0.01)
            dist1 = PointDist((x, y), P)
            x2 = Bx.subs(t, param - 0.01)
            y2 = By.subs(t, param - 0.01)
            dist2 = PointDist((x, y), P)
        if (param < 0.01 and param > 0.99) :
            return min, minP
        if (dist1 < min and dist1 < dist2) :
            param += 0.01
            min = dist1
            minP[0] = x1
            minP[1] = y1
            continue
        elif (dist2 < min and dist2 < dist1) :
            param -= 0.01
            min = dist2
            minP[0] = x2
            minP[1] = y2
            continue
        else:
            return min, minP
    return min, minP

#Min distance from given point to Bezier curve
#given: cnt - nr of Bezier curve segment
#       B - Bezier curve on which we need to find point on
#       P - given point from which the closest line will be calculated
#       minDist - that has been found already
#       len - the length of B array
#       minPtmp - just some point, that we dont need to create every time the function is called
#       minP - coordinates of the closest point on Bezier curve to P
def MinDistTanBezier(cnt, B, P, minDist, len, minPtmp, minP) :
    if cnt > len - 1 :
        return len - 1, minP, minDist

    dist, minPtmp = BezierMinDist(B[cnt][0], B[cnt][1], P)
    if minDist >= dist :
        minDist = dist
        minP = minPtmp
        cnt += 1
        return MinDistTanBezier(cnt, B, P, minDist, len)
    elif dist > minDist :
        cnt -= 1
        return cnt, minP, minDist
    return cnt, minP, minDist

visited = {}

for i in line1 :
    if not i[0] in visited :
        visited = DiffAll(i[0], visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
print("aaa")
for i in line2 :
    if not i[0] in visited :
        visited = OneFileObject(i[0], line2, Bezier2, CP2, visited, 2)
t = Symbol('t')
for i in Bezier1 :
    if not (i[0][0].subs(t, 0), i[0][1].subs(t, 0)) in visited :
        visited = DiffAll((i[0][0].subs(t, 0), i[0][1].subs(t, 0)), visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
for i in Bezier2 :
    if not (i[0][0].subs(t, 0), i[0][1].subs(t, 0)) in visited :
        visited = OneFileObject((i[0][0].subs(t, 0), i[0][1].subs(t, 0)), line2, Bezier2, CP2, visited, 2)

plt.show()
