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


#find closest point from specified objects
def ClosestPntLine(P, arr):
    min = -1
    minP = 0
    for i in range(len(arr) - 1) :
        if i == 0 :
            min = PointDist(P, arr[i][0])
            minP = i
            lineSt = 0
        tmpD = PointDist(P, arr[i][0])
        if tmpD < min :
            min = tmpD
            minP = i
            lineSt = 0
        tmpD = PointDist(P, arr[i][1])
        if tmpD < min :
            min = tmpD
            minP = i
            lineSt = 1
    return arr[minP][lineSt], min

def ClosestPntBezier(P, arr):
    min = -1
    minP = 0
    t = Symbol('t')
    for i in range(len(arr) - 1) : #as the Bezier curve is made from multiple cubic Bezier curves, we have to compare only the start point of first and end point of last cubic curve that belongs to specific composite Bezier
        if i == 0 :
            min = PointDist(P, (arr[i][0][0].subs(t, 0), arr[i][0][1].subs(t, 0)))
            minP = i
            lineSt = 0
        tmpD = PointDist(P, (arr[i][0][0].subs(t,0), arr[i][0][1].subs(t,0)))
        if tmpD < min :
            min = tmpD
            minP = i
            lineSt = 0
        tmpD = PointDist(P, (arr[i][-1][0].subs(t, 1), arr[i][-1][1].subs(t, 1)))
        if tmpD < min :
            min = tmpD
            minP = i
            lineSt = 1
    return arr[minP][0], min


#color = '#6c9f92' (color for largest distane lines)

def DiffAll(P1, visited, line1, line2, Bezier1, Bezier2, CP1, CP2) :
    if P1 in visited :
        if visited[P1] == 2 :
            return visited
    visited[P1] = 1
    P21, min21 = ClosestPntLine(P1, line2)
    P22, min22 = ClosestPntBezier(P1, Bezier2)
    if (min21 < min22) :
        P2 = P21
    else:
        P2 = P22

    if (min21 == -1 and min22 == -1) :
        visited[P1] = 2
        return

    for i in range(len(line1)) :
        P1ind = -1
        lineSt1 = -1
        if P1 == line1[i][0] :
            P1ind = i
            lineSt1 = 0
        elif P1 == line1[i][1] :
            P1ind = i
            lineSt1 = 1
        else:
            continue

        a1 = LineSlope(line1[i][0], line1[i][1])

        lineSt2 = -1
        P2ind = -1
        minEndDist = -1 #to find which line is the closest one to the line in line1
        for j in range(len(line2)):
            if P2 == line2[j][0] :
                lineSt2 = 0
            elif P2 == line2[j][1] :
                lineSt2 = 1
            else:
                continue
            a2 = LineSlope(line2[j][0], line2[j][1])
            angle = np.arctan((a2 - a1) / (1 + abs(a1 * a2)))
            if abs(angle) > 5 :
                continue
            if (lineSt1 != -1 and lineSt2 != -1) :
                endDist = PointDist(line1[i][(lineSt1 + 1) % 2], line2[j][(lineSt2 + 1) % 2])
                if minEndDist == -1 :
                    minEndDist = endDist
                    P2Ind = j
                if endDist < minEndDist :
                    minEndDist = endDist
                    P2ind = j

        if (P2ind == -1) :
            print(line1[i][0], line1[i][1])
            plt.plot([line1[i][0][0], line1[i][1][0]], [line1[i][0][1], line1[i][1][1]], color = '#055583', linestyle = 'dotted')
            visited[(line1[i][0], line1[i][1])] = 2
            if not line1[i][(lineSt1 + 1) % 2] in visited :
                visited = DiffAll(line1[i][(lineSt1 + 1) % 2], visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
        if (P1ind != -1 and P2ind != -1) :
            dist1 = PointDist(line1[i][lineSt1], line2[P2Ind][lineSt2])
            dist2 = PointDist(line1[i][(lineSt1 + 1) % 2], line2[P2Ind][(lineSt2 + 1) % 2])
            plt.plot([line1[i][lineSt1][0], line2[P2Ind][lineSt2][0]], [line1[i][lineSt1][1], line2[P2Ind][lineSt2][1]], color = '#6c9f92')
            plt.plot([line1[i][(lineSt1 + 1) % 2][0],line2[P2Ind][(lineSt2 + 1) % 2][0]], [line1[i][(lineSt1 + 1) % 2][1], line2[P2Ind][(lineSt2 + 1) % 2][1]], color = '#6c9f92')
            if  not line1[i][(lineSt1 + 1) % 2] in visited :
                visited = DiffAll(line1[i][(lineSt1 + 1) % 2], visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
            if not line2[P2Ind][(lineSt2 + 1) % 2] in visited :
                visited = DiffAll(line2[P2Ind][(lineSt2 + 1) % 2], visited, line1, line2, Bezier1, Bezier2, CP1, CP2)


    for i in range(len(Bezier1)) :
        P1ind = -1
        BezierSt1 = -1
        t = Symbol('t')
        if P1 == (Bezier1[i][0][0].subs(t, 0), Bezier1[i][0][1].subs(t, 0)) :
            P1ind = i
            BezierSt1 = 0
        elif P1 == (Bezier1[i][-1][0].subs(t, 1), Bezier1[i][-1][1].subs(t, 1)) :
            P1ind = i
            BezierSt1 = 1
        else:
            continue
        a1 = LineSlope((Bezier1[i][0][0].subs(t, 0), Bezier1[i][0][1].subs(t, 0)), (Bezier1[i][-1][0].subs(t, 1), Bezier1[i][-1][1].subs(t, 1)))

        BezierSt2 = -1
        P2ind = -1
        for j in range(len(Bezier2)):
            if P2 == (Bezier2[j][0][0].subs(t, 0), Bezier2[j][0][1].subs(t, 0)) :
                BezierSt2 = 0
            elif P2 == (Bezier2[j][-1][0].subs(t, 1), Bezier2[j][-1][1].subs(t, 1)) :
                BezierSt2 = 1
            else:
                continue
            a2 = LineSlope((Bezier1[j][0][0].subs(t, 0), Bezier1[j][0][1].subs(t, 0)), (Bezier1[j][-1][0].subs(t, 1), Bezier1[j][-1][1].subs(t, 1)))
            angle = np.arctan((a2 - a1) / (1 + abs(a1 * a2)))
            if abs(angle) > 5 :
                continue
            if (BezierSt1 != -1 and BezierSt2 != -1) :
                endDist = PointDist((Bezier1[i][-(BezierSt1 + 1) % 2][0].subs(t, (BezierSt1 + 1) % 2), Bezier1[i][-(BezierSt1 + 1) % 2][1].subs(t, (BezierSt1 + 1) % 2)), (Bezier2[j][-(BezierSt2 + 1) % 2][0].subs(t, (BezierSt2 + 1) % 2),Bezier2[j][-(BezierSt2 + 1) % 2][1].subs(t, (BezierSt2 + 1) % 2)))
                if minEndDist == -1 :
                    minEndDist = endDist
                    P2Ind = j
                if endDist < minEndDist :
                    minEndDist = endDist
                    P2ind = j

        if (P2ind == -1) :
            PlotBezier(CP1[i], '#055583', 1, 'dotted')
            t = Symbol('t')
            if not (Bezier1[i][-(BezierSt1 + 1) % 2][0].subs(t, (BezierSt1 + 1) % 2), Bezier1[i][-(BezierSt1 + 1) % 2][1].subs(t, (BezierSt1 + 1) % 2)) in visited :
                visited = DiffAll(((Bezier1[i][-(BezierSt1 + 1) % 2][0].subs(t, (BezierSt1 + 1) % 2)), Bezier1[i][-(BezierSt1 + 1) % 2][1].subs(t, (BezierSt1 + 1) % 2)), visited, line1, line2, Bezier1, Bezier2, CP1, CP2)

        if (P1ind != -1 and P2ind != -1) :
            len1 = 0
            for k in range(len(Bezier1[i])) :
                len1 += CubicBezierLen(CP1[i])
            len2 = 0
            for k in range(len(Bezier2[P2ind])) :
                len2 += CubicBezierLen(CP2[P2ind])

            int1 = len1 / 100
            int2 = len2 / 100

            BezierDiff(Bezier1[i], int1, Bezier2[P2ind], int2)

#####Šis also jāsalabo!!!!!
            MaxBezierDist(Bezier1[i], Bezier2[P2ind])
            if not ((Bezier1[i][0].subs(t, (BezierSt1 + 1) % 2)), Bezier1[i][1].subs(t, (BezierSt1 + 1) % 2)) in visited :
                visited = DiffAll((Bezier1[i][0].subs(t, (BezierSt1 + 1) % 2)), visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
            if not ((Bezier2[P2ind][0].subs(t, (BezierSt2 + 1) % 2)), Bezier2[P2Ind][1].subs(t, (BezierSt2 + 1) % 2)) in visited :
                visited = DiffAll((Bezier2[P2Ind][0].subs(t, (BezierSt2 + 1) % 2), Bezier2[P2Ind][1].subs(t, (BezierSt2 + 1) % 2) ), visited, line1, line2, Bezier1, Bezier2, CP1, CP2)

        for i in line2 :
            if P2 == i[0] :
                plt.plot([i[0][0], i[1][0]], [i[0][1],i[1][1]], color = '#bc0e13')
                if not i[1] in visited :
                    visited = DiffAll(i[1], visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
            if P2 == i[1] :
                plt.plot([i[0][0], i[1][0]], [i[0][1], i[1][1]], color = '#bc0e13')
                if not i[0] in visited :
                    visited = DiffAll(i[0], visited, line1, line2, Bezier1, Bezier2, CP1, CP2)

        k = 0
        for i in Bezier2 :
            if P1 == (i[0][0].subs(t, 0), i[0][1].subs(t, 0)) :
                PlotBezier(CP2[k], '#bc0e13', 1, 'dotted')
                if not (i[-1][0].subs(t, 1), i[-1][1].subs(t, 1)) in visited :
                    visited = DiffAll((i[1][0].subs(t, 1), i[1][1].subs(t, 1)), visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
            if P1 == (i[-1][0].subs(t, 1), i[-1][1].subs(t, 1)) :
                PlotBezier(CP2[k], '#bc0e13', 1, 'dotted')
                if not (i[0][0].subs(t, 0), i[0][1].subs(t, 0)) in visited :
                    visited = DiffAll((i[0][0].subs(t, 0), i[0][1].subs(t, 0)), visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
            ++k

    visited[P1] = 2
    return visited


def CubicBezierLen(C):
    ax = - 3 * C[0][0] + 9 * C[1][0] - 9 * C[2][0] + 3 * C[3][0]
    bx = 6 * C[0][0] - 12 * C[1][0] + 6 * C[2][0]
    cx = - 3 * C[0][0] + 3 * C[1][0]

    ay = - 3 * C[0][1] + 9 * C[1][1] - 9 * C[2][1] + 3 * C[3][1]
    by = 6 * C[0][1] - 12 * C[1][1] + 6 * C[2][1]
    cy = - 3 * C[0][1] + 3 * C[1][1]

    ft = np.sqrt((ax**2 + ay**2) * t**4 + 2 * (ax * bx + ay * by) * t**3 + (2 * ax * cx + bx**2 + 2 * ay * cy + by**2) * t**2 + 2 * (bx * cx + by * cy) * t + (bx**2 + by**2))

    # from https://pomax.github.io/bezierinfo/legendre-gauss.html
    w[1] = 0.8888888888888888
    w[2] = 0.5555555555555556
    w[3] = 0.5555555555555556

    x[1] = 0.0000000000000000
    x[2] = -0.7745966692414834
    x[3] = 0.7745966692414834

    L = 0
    for i in range(1, 4) :
        L += 0.5 * (w[i] * ft.subs(t, 0.5 * x[i] + 0.5))

    return L

def MaxBezierDist(B1, B2):
    for i in range(len(B1)) :
        cnt = 0
        param = 0
        max = 0
        while param <= 1 :
            param += cnt
            cnt += 0.1
            x2 = Bx2.subs(t, param)
            y2 = By2.subs(t, param)
            dist, P = BezierMinDist(Bx1, By1, (x2, y2))
            if cnt == 0 :
                P1 = []
                P1.append([Bx1, By1])
                P2 = []
                P2. append([x2, y2])
            if dist > max :
                P1[0] = P[0]
                P1[1] = P[1]
                P2[0] = x2
                P2[1] = y2
                max = dist

    for i in range(len(B2)) :
        cnt = 0
        param = 0
        while param <= 1 :
            param += cnt
            cnt += 0.1
            x1 = Bx1.subs(t, param)
            y1 = By1.subs(t, param)
            dist = BezierMinDist(Bx2, By2, (x1, y1))
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
    param = 0
    #get distance values from points on bezier that differs from each other by 0.1 parameter value
    for i in range (11) :
        param = param + i * 0.1
        x = Bx.subs(t, param)
        y = By.subs(t, param)
        dist = PointDist((x, y), P)
        if i == 0 :
            min = PointDist((x, y), P)
            minParam = param
        elif dist < min :
            min = dist
            minParam = param
    #check if in segments between previously determined points there is not point with smaller distance to the given point
    while (param >= minParam - 0.1 and param <= minParam + 0.1) :
        x = Bx.subs(t, param + 0.1)
        y = By.subs(t, param + 0.1)
        dist1 = PointDist((x, y), P)
        x = Bx.subs(t, param - 0.1)
        y = By.subs(t, param - 0.1)
        dist2 = PointDist((x, y), P)
        if (dist1 < min and dist1 < dist2) :
            param += 0.1
            min = dist1
            continue
        elif (dist2 < min and dist2 < dist1) :
            param -= 0.1
            min = dist2
            continue
        else:
            x = Bx.subs(t, param)
            y = By.subs(t, param)
            plt.plot((x, y), P, color = '#6c9f92')
            return min

#calculate difference value using square method
def BezierDiff(B1, int1, lineSt1, B2, int2, lineSt2, CP1, CP2) :
    dist1 = 0
    dist2 = 0
    value = 0
    nr1 = - lineSt1
    nr2 = - lineSt2
    for i in range(101) :
        P1, nr1 = DistantPointOnBezier(dist1, nr1, B1, lineSt1, CP1, 1)
        P2, nr2 = DistantPointOnBezier(dist2, nr2, B2, lineSt2, CP2, 2)
        value += DistantPoint
        dist1 += int1
        dist2 += int2
    return value


def DistantPointOnBezier(dist, nr, B, lineSt, C, fnr) :
    xtmp = B[nr][0].subs(t, lineSt)
    ytmp = B[nr][1].subs(t, lineSt)
    i = 0
    param = lineSt
    cnt = 0.1
    m = 0.1
    disttmp = 0
    while (nr >= 0) and (nr <= len(B) - 1) :
        while (param != 0 and param != 1) :
            param = abs(lineSt - cnt)
            x2tmp = B[nr][0].subs(t, param)
            y2tmp = B[nr][1].subs(t, param)
            cnt += m
            disttmp += PointDist((xtmp, ytmp), (x2tmp, y2tmp))
            if (disttmp < dist + 0.1 and disttmp > dist - 0.1):
                return (x2tmp, y2tmp)
            elif disttmp > dist + 0.1 :
                cnt -= 0.1
                m *= 0.1
            xtmp = x2tmp
            ytmp = y2tmp

        dist1 = dist1 + CubicBezierLen(C[fnr][nr])
        disttmp = dist1
        if lineSt == 0 :
            ++nr
        else:
            --nr


visited = {}
for i in line1 :
    if not i[0] in visited :
        visited = DiffAll(i[0], visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
for i in line2 :
    if not i[0] in visited :
        visited = DiffAll(i[0], visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
t = Symbol('t')
for i in Bezier1 :
    if not (i[0][0].subs(t, 0), i[0][1].subs(t, 0)) in visited :
        visited = DiffAll((i[0][0].subs(t, 0), i[0][1].subs(t, 0)), visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
for i in Bezier2 :
    if not (i[0][0].subs(t, 0), i[0][1].subs(t, 0)) in visited :
        visited = DiffAll((i[0][0].subs(t, 0), i[0][1].subs(t, 0)), visited, line1, line2, Bezier1, Bezier2, CP1, CP2)


plt.show()
