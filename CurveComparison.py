#PF problem - (something with const y line)
import dxfgrabber
import numpy as np
from numpy.linalg import inv
import array as arr
import sympy
from sympy import init_printing, Symbol, UnevaluatedExpr, expand, pretty_print as pprint, latex
from sympy.solvers import solve
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.axes as axes


#Filenames should be written in the following two lines, as shown in example:
#   dxf1(2) = dxfgrabber.readfile("filename")
dxf1 = dxfgrabber.readfile("ieva_likne2_1_002.dxf")
dxf2 = dxfgrabber.readfile("ieva_likne2_002.dxf")


firstFileCol = '#055583'
secondFileCol = '#bc0e13'

#Graph formatting
def MyForm(x, p):
    return x / 10

plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
major_ticks = np.arange(-1000, 1000, 100)
minor_ticks = np.arange(-1000, 1000, 10)
ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor = True)
ax.set_yticks(major_ticks)
ax.set_yticks(minor_ticks, minor = True)
ax.grid(which = 'both')
ax.grid(which = 'major', alpha = 0.5)
ax.grid(which = 'minor', alpha = 0.2)
ax.get_xaxis().set_major_formatter(
    mpl.ticker.FuncFormatter(MyForm))
ax.get_yaxis().set_major_formatter(
    mpl.ticker.FuncFormatter(MyForm))

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

#type of objects in file
type1 = [entity.dxftype for entity in dxf1.entities]
output1 = [entity for entity in dxf1.entities]
Color1 = [entity.color for entity in output1]
Linetype1 = [entity.linetype for entity in output1 if entity.dxftype != 'POINT']
CP1 = []
points1 = []
line1 = []
Bezier1 = []
circle1 = []
ellipse1 = []
arc1 = []
LWPolyline1 = []
splineCP1 = []
object1nr = 1.000
object2nr = 2.000
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
        centerPoints1 = entity.center
        radius1 = entity.radius
        lenCirclePnts = len(centerPoints1)
        for i in lenCirclePnts :
            circle1.append(centerPoints1[i], radius1[i])

    #Ellipse
    if entity.dxftype == 'ELLIPSE':
        ellipseCenter1 = entity.center
        ellipseMajorAxis1 = entity.major_axis
        ellipse.append([ellipseCenter1, ellipseMajorAxis1])

    #Arc
    if entity.dxftype == 'ARC':
        arcCenter1 = entity.center
        arcRadius1 = entity.radius
        arcStartAngle1 = entity.start_angle
        arcEndAngle1 = entity.end_angle
        arc1.append([arcCenter1, arcRadius1, arcStartAngle1, arcEndAngle1])
    #Polyline
    if entity.dxftype == 'POLYLINE':
        PolylinePoints1 = entity.points
        Beziertmp, CP = CompositeBezier(PolylinePoints1, 1)
        Bezier1.append(Beziertmp)
        CP1.append(CP)

    #LWPolyline
    if entity.dxftype == 'LWPOLYLINE':
        LWPolylinePoints1 = entity.points
        Beziertmp, CP = CompositeBezier(LWPolylinePoints1, 1)
        Bezier1.append(Beziertmp)
        CP1.append(CP)

    #Spline
    if entity.dxftype == 'SPLINE':
        splineControlPoints1 = entity.control_points
        splineCP1.append(splineControlPoints1)
#Block
#print("BLOCKS")
#BlockBasepoint= [Block.BlockBasepoint for Block in dxf.entities]
#print(BlockBasepoint)
#Bezier formulas of file #2

#type of objects in file
type2 = [entity.dxftype for entity in dxf2.entities]
output2 = [entity for entity in dxf2.entities]
Color2 = [entity.color for entity in output2]
Linetype2 = [entity.linetype for entity in output2 if entity.dxftype != 'POINT']
CP2 = []
points2 = []
line2 = []
Bezier2 = []
circle2 = []
ellipse2 = []
arc2 = []
LWPolyline2 = []
splineCP2 = []
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
        centerPoints2 = entity.center
        radius2 = entity.radius
        lenCirclePnts = len(centerPoints2)
        for i in lenCirclePnts :
            circle2.append(centerPoints2[i], radius2[i])
    #Ellipse
    if entity.dxftype == 'ELLIPSE':
        EllipseCenter2 = entity.center
        EllipseMajorAxis2 = entity.major_axis
    #Arc
    if entity.dxftype == 'ARC':
        arcCenter2 = entity.center
        arcRadius2 = entity.radius
        arcStartAngle2 = entity.start_angle
        arcEndAngle2 = entity.end_angle
        arc2.append([arcCenter2, arcRadius2, arcStartAngle2, arcEndAngle2])
    #Polyline
    if entity.dxftype == 'POLYLINE':
        PolylinePoints2 = entity.points
        Beziertmp, CP = CompositeBezier(PolylinePoints2, 2)
        Bezier2.append(Beziertmp)
        CP2.append(CP)
    #LWPolyline
    if entity.dxftype == 'LWPOLYLINE':
        LWPolylinePoints2 = entity.points
        Beziertmp, CP = CompositeBezier(LWPolylinePoints2, 2)
        Bezier2.append(Beziertmp)
        CP2.append(CP)
    #Spline
    if entity.dxftype == 'SPLINE':
        splineControlPoints2 = entity.control_points
        splineCP2.append(splineControlPoints2)

for i in CP1 :
    PlotBezier(i, '#055583', 0.2, 'solid')
for i in CP2 :
    PlotBezier(i, '#bc0e13', 0.2, 'solid')
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
def rotate(points1, points2, line1, line2, Bezier1, Bezier2, arc1, arc2, circle1, circle2, ellipse1, ellipse2) :
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
uniqueObjectsLine1 = {}
uniqueObjectsLine2 = {}
ObjectsLine1 = {}
ObjectsLine2 = {}
uniqueObjectsBezier1 = {}
uniqueObjectsBezier2 = {}
ObjectsBezier1 = {}
ObjectsBezier2 = {}
uniqueObjectsArc1 = {}
uniqueObjectsArc2 = {}
ObjectsArc1 = {}
ObjectsArc2 = {}
uniqueObjectsCircle1 = {}
uniqueObjectsCircle2 = {}
ObjectCircle = {}
ObjectCircle1 = {}
ObjectCircle2 = {}
uniqueObjectsEllipse1 = {}
uniqueObjectsEllipse2 = {}
ObjectEllipse1 = {}
ObjectEllipse2 = {}

equalObjectsLine = {}
equalObjectsBezier = {}


#find closest point in lines array
def ClosestPntLine(P, arr, visited):
    min = -1
    minPind = 0
    lineSt = 0
    for i in range(len(arr)) :
        if min == -1 :
            if not arr[i][0] in visited :
                min = PointDist(P, arr[i][0])
                minPind = i
                lineSt = 0
            else :
                if visited[arr[i][0]] == 1 :
                    min = PointDist(P, arr[i][0])
                    minPind = i
                    lineSt = 0
        if not arr[i][0] in visited :
            tmpD = PointDist(P, arr[i][0])
            if tmpD < min :
                min = tmpD
                minPind = i
                lineSt = 0
        else :
            if visited[arr[i][0]] == 1 :
                tmpD = PointDist(P, arr[i][0])
                if tmpD < min :
                    min = tmpD
                    minPind = i
                    lineSt = 0
        if not arr[i][1] in visited :
            tmpD = PointDist(P, arr[i][1])
            if tmpD < min :
                min = tmpD
                minPind = i
                lineSt = 1
        else :
            if visited[arr[i][1]] == 1 :
                tmpD = PointDist(P, arr[i][1])
                if tmpD < min :
                    min = tmpD
                    minPind = i
                    lineSt = 1
    return arr[minPind][lineSt], min


#find closest point in Bezier curves array
def ClosestPntBezier(P, arr, visited) :
    min = -1
    minPind = 0
    t = Symbol('t')
    lineSt = 0
    for i in range(len(arr)) :
        #as the Bezier curve is made from multiple cubic Bezier curves, we have to compare only the start point of first and end point of last cubic curve that belongs to specific composite Bezier
        if min == -1 :
            if not (arr[i][0][0].subs(t, 0), arr[i][0][1].subs(t, 0)) in visited:
                min = PointDist(P, (arr[i][0][0].subs(t, 0), arr[i][0][1].subs(t, 0)))
                minPind = i
                lineSt = 0
            elif visited[(arr[i][0][0].subs(t, 0), arr[i][0][1].subs(t, 0))] == 1 :
                min = PointDist(P, (arr[i][0][0].subs(t, 0), arr[i][0][1].subs(t, 0)))
                minPind = i
                lineSt = 0
        if not (arr[i][0][0].subs(t,0), arr[i][0][1].subs(t,0)) in visited :
            tmpD = PointDist(P, (arr[i][0][0].subs(t,0), arr[i][0][1].subs(t,0)))
            if tmpD < min :
                min = tmpD
                minPind = i
                lineSt = 0
        elif visited[(arr[i][0][0].subs(t,0), arr[i][0][1].subs(t,0))] == 1 :
            tmpD = PointDist(P, (arr[i][0][0].subs(t,0), arr[i][0][1].subs(t,0)))
            if tmpD < min :
                min = tmpD
                minPind = i
                lineSt = 0

        if not (arr[i][-1][0].subs(t, 1), arr[i][-1][1].subs(t, 1)) in visited :
            tmpD = PointDist(P, (arr[i][-1][0].subs(t, 1), arr[i][-1][1].subs(t, 1)))
            if tmpD < min :
                min = tmpD
                minPind = i
                lineSt = 1
        elif visited[(arr[i][-1][0].subs(t, 1), arr[i][-1][1].subs(t, 1))] == 1 :
            tmpD = PointDist(P, (arr[i][-1][0].subs(t, 1), arr[i][-1][1].subs(t, 1)))
            if tmpD < min :
                min = tmpD
                minPind = i
                lineSt = 1
    return (arr[minPind][-lineSt][0].subs(t, lineSt), arr[minPind][-lineSt][1].subs(t, lineSt)), min

def OneFileObject(P, line, Bezier, CP, visited, fnr) :
    lenLine = len(line)
    for i in range(lenLine) :
        if (P == line[i][0] or P == line[i][1]) :
            if fnr == 1 :
                plt.plot((line[i][0][0], line[i][1][0]), (line[i][0][1], line[i][1][1]), color = '#055583', linestyle = 'dotted')
                if uniqueObjectsLine1.get(i) == None :
                    length = PointDist(line[i][0], line[i][1])
                    uniqueObjectsLine1[i] = length
            if fnr == 2 :
                plt.plot((line[i][0][0], line[i][1][0]), (line[i][0][1], line[i][1][1]), color = '#bc0e13', linestyle = 'dotted')
                if uniqueObjectsLine2.get(i) == None :
                    length = PointDist(line[i][0], line[i][1])
                    uniqueObjectsLine2[i] = length
            if not line[i][0] in visited :
                visited[line[i][0]] = 1
            if not line[i][1] in visited :
                visited[line[i][1]] = 1

    t = Symbol('t')
    lenB = len(Bezier)
    for i in range(lenB) :
        if (P == (Bezier[i][0][0].subs(t, 0), Bezier[i][0][1].subs(t, 0)) or (P == (Bezier[i][-1][0].subs(t, 1), Bezier[i][-1][1].subs(t, 1)))) :
            if fnr == 1 :
                PlotBezier(CP[i], '#055583', 1, 'dotted')
                if uniqueObjectsBezier1.get(i) == None :
                    length = CompositeBezierLen(CP[i])
                    uniqueObjectsBezier1[i] = length
            if fnr == 2 :
                PlotBezier(CP[i], '#bc0e13', 1, 'dotted')
                if uniqueObjectsBezier2.get(i) == None :
                    length = CompositeBezierLen(CP[i])
                    uniqueObjectsBezier2[i] = length
            if not (Bezier[i][0][0].subs(t, 0), Bezier[i][0][1].subs(t, 0)) in visited :
                visited[(Bezier[i][0][0].subs(t, 0), Bezier[i][0][1].subs(t, 0))] = 1
            if not (Bezier[i][-1][0].subs(t, 1), Bezier[i][-1][1].subs(t, 1)) in visited :
                visited[Bezier[i][-1][0].subs(t, 1), Bezier[i][-1][1].subs(t, 1)] = 1
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
    param = lineSt
    cnt = 1
    m = 0.5
    disttmp = 0
    dist1 = 0
    lenB = len(B)
    while (nr >= 0) and (nr <= lenB - 1) :
        while (param >= 0 and param <= 1) :
            param = abs(round((lineSt - cnt), 5))
            if (param < 0 or param > 1) :
                if param < 0 :
                    param = 0
                else :
                    param = 1
                break
            x2tmp = B[nr][0].subs(t, param)
            y2tmp = B[nr][1].subs(t, param)
            Pdist = PointDist((xtmp, ytmp), (x2tmp, y2tmp))
            disttmp += Pdist
            if (disttmp <= dist + 0.05 and disttmp >= dist - 0.05) :
                return x2tmp, y2tmp, nr
            elif m == 0 :
                break
            elif disttmp > dist + 0.05 :
                cnt = round((cnt - m), 5)
                m = round(m * 0.5, 5)
                disttmp -= Pdist
                continue
            xtmp = x2tmp
            ytmp = y2tmp
            if (param == 0 or param == 1) :
                break
            cnt = round((cnt + m), 5)
        dist1 = dist1 + CubicBezierLen(C[nr])
        disttmp = dist1
        if (nr == lenB - 1 and lineSt == 0) :
            break
        if (nr == 0 and lineSt == 1) :
            break
        if lineSt == 0 :
            nr += 1
            param = lineSt
            m = 0.5
            cnt = 1
        else:
            nr -= 1
            param = lineSt
            m = 0.5
            cnt = 1
    x2tmp = B[nr][0].subs(t, param)
    y2tmp = B[nr][1].subs(t, param)
    return x2tmp, y2tmp, nr

#Calculate LeastSquare value between two given lines
def LeastSquare(B1, int1, lineSt1, B2, int2, lineSt2, CP1, CP2) :
    dist1 = 0
    dist2 = 0
    value = 0
    lenB1 = len(B1)
    lenB2 = len(B2)
    for i in range(11) :
        if lineSt1 == 1 :
            nr1 = lenB1 - 1
        else :
            nr1 = 0
        if lineSt2 == 1 :
            nr2 = lenB2 - 1
        else :
            nr2 = 0
        P1x, P1y, nr1 = DistantPointOnBezier(dist1, nr1, B1, lineSt1, CP1)
        P2x, P2y, nr2 = DistantPointOnBezier(dist2, nr2, B2, lineSt2, CP2)
        value += (PointDist((P1x, P1y), (P2x, P2y)))**2
        dist1 += int1
        dist2 += int2
#        plt.plot(P1x, P1y, 'ro')
#        plt.plot(P2x, P2y, 'bo')
#        plt.plot([P1x, P2x], [P1y, P2y], color = '#af5ba3')
    t = Symbol('t')
    a = B1[0][0].subs(t, 0)
    b = B1[0][1].subs(t, 0)
    ax.annotate(value, xy = (a, b), xytext = (a, b))
    return value


def MaxBezierDist(B1, B2) :
    t = Symbol('t')
    lenB1 = len(B1)
    lenB2 = len(B2)
    max = 0
    P = []
    P.append([])
    P.append([])
    P1 = []
    P1.append([])
    P1.append([])
    P2 = []
    P2.append([])
    P2.append([])
    P11 = []
    P11.append([])
    P11.append([])
    P22 = []
    P22.append([])
    P22.append([])
    for i in range(lenB2) :
        cnt = 0
        param = 0
        while param <= 1 :
            cnt += 0.1
            x2 = B2[i][0].subs(t, param)
            y2 = B2[i][1].subs(t, param)
            min = 10000
            for j in range(lenB1) :
                dist, P = BezierMinDist(B1[j][0], B1[j][1], (x2, y2))
                if dist < min :
                    P1[0] = P[0]
                    P1[1] = P[1]
                    P2[0] = x2
                    P2[1] = y2
                    min = dist
            if min >= max :
                max = min
                P11[0] = P1[0]
                P11[1] = P1[1]
                P22[0] = P2[0]
                P22[1] = P2[1]

            param += cnt
        #print(P1, P2)
        #plt.plot([P1[0], P2[0]], [P1[1], P2[1]])
    #plt.plot([P1[0], P2[0]], [P1[1], P2[1]], color = '#667281')
    a = ((P11[0] + P22[0]) / 2)
    b = ((P11[1] + P22[1]) / 2)
    ax.annotate(round(max / 10, 2), xy = (a, b), xytext = (a, b), alpha = 0.5, size = 'small')
    plt.plot([P11[0], P22[0]], [P11[1], P22[1]], color = 'black', alpha = 0.4)
    return max

#Search for distance between point and Bezier curve
def BezierMinDist(Bx, By, P) :
    t = Symbol('t')
    param = 0
    #get distance values from points on bezier that differs from each other by 0.1 parameter value
    for i in range(11) :
        param = 0.1 * i
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
        if param <= 0.99:
            x1 = Bx.subs(t, param + 0.01)
            y1 = By.subs(t, param + 0.01)
            dist1 = PointDist((x1, y1), P)
        else :
            dist1 = 10000
        if param >= 0.01 :
            x2 = Bx.subs(t, param - 0.01)
            y2 = By.subs(t, param - 0.01)
            dist2 = PointDist((x2, y2), P)
        else :
            dist2 = 10000
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
        if (param < 0.01 and param > 0.99) :
            return min, minP
        else:
            return min, minP

    #This point technically should never been met
    return min, minP

def DiffAll(P1, visited, line1, line2, Bezier1, Bezier2, CP1, CP2) :
    DrawnLines2Ind = {}
    #Check if given point isn`t already visited

    if P1 in visited :
        if visited[P1] == 2 :
            return visited

    #Get coordinates of the closest point in second file
    if len(line2) > 0 :
        P21, min21 = ClosestPntLine(P1, line2, visited)
    else :
        min21 = -1
    if len(Bezier2) > 0 :
        P22, min22 = ClosestPntBezier(P1, Bezier2, visited)
    else :
        min22 = -1
    if  (min21 == -1 and min22 == -1) :
        visited = OneFileObject(P1, line1, Bezier1, CP1, visited, 1)
        return visited
    #Select closest point
    if (min21 < min22 and min21 != -1) :
        P2 = P21
    elif (min21 > min22 and min22 != -1) :
        P2 = P22
    elif min21 != -1 :
        P2 = P21
    elif min22 != -1 :
        P2 = P22

    #if the point is too far, then probably that is not the respective point
    if PointDist(P1, P2) > 15 :
        visited = OneFileObject(P1, line1, Bezier1, CP1, visited, 1)
        return visited

    #Find all lines with one of the enpoints matching the given point
    lenLine1 = len(line1)
    lenLine2 = len(line2)
    for i in range(lenLine1) :
        P1Ind = -1
        lineSt1 = -1
        #if this one doesn`t correspond, then ignore
        if P1 == line1[i][0] :
            if not line1[i][0] in visited :
                if not line1[i][1] in visited :
                    P1Ind = i
                    lineSt1 = 0
                    visited[line1[i][1]] = 1
                elif visited[line1[i][1]] == 1 :
                    P1Ind = i
                    lineSt1 = 0
            elif visited[line1[i][0]] == 1 :
                if not line1[i][1] in visited :
                    P1Ind = i
                    lineSt1 = 0
                elif visited[line1[i][1]] == 1 :
                    P1Ind = i
                    lineSt1 = 0
        elif P1 == line1[i][1] :
            if not line1[i][1] in visited :
                if not line1[i][0] in visited :
                    P1Ind = i
                    lineSt1 = 1
                elif visited[line1[i][0]] == 1 :
                    P1Ind = i
                    lineSt1 = 1
            elif visited[line1[i][0]] == 1 :
                if not line1[i][0] in visited :
                    P1Ind = i
                    lineSt1 = 1
                elif visited[line1[i][0]] == 1 :
                    P1Ind = i
                    lineSt1 = 1
        else:
            continue
        print(line1[P1Ind])
        a1 = LineSlope(line1[P1Ind][0], line1[P1Ind][1])
        lineSt2 = -1
        P2Ind = -1
        minEndDist = -1 #to find which line is the closest one to the line in line1
        #find the respective line in second file
        for j in range(lenLine2):
            lineSt2tmp = -1
            #if endpoints doesn`t match P2 then ignore and iterate forward
            if P2 == line2[j][0] :
                if not line2[j][0] in visited :
                    if not line2[j][1] in visited :
                        lineSt2tmp = 0
                    elif visited[line2[j][1]] == 1 :
                        lineSt2tmp = 0
                elif visited[line2[j][0]] == 1 :
                    if not line2[j][1] in visited :
                        lineSt2tmp = 0
                    elif visited[line2[j][1]] == 1 :
                        lineSt2tmp = 0
            elif P2 == line2[j][1] :
                if not line2[j][1] in visited :
                    if not line2[j][0] in visited :
                        lineSt2tmp = 1
                    elif visited[line2[j][0]] == 1 :
                        lineSt2tmp = 1
                elif visited[line2[j][1]] == 1 :
                    if not line2[j][0] in visited :
                        lineSt2tmp = 1
                    elif visited[line2[j][0]] == 1 :
                        lineSt2tmp = 1
            else:
                continue
            print("((((((((((()))))))))))")
            print(line2[j])

            a2 = LineSlope(line2[j][0], line2[j][1])
            angle = np.arctan((a2 - a1) / (1 + abs(a1 * a2)))
            #if the angle between these two lines is bigger than 5, then we eill assume that it is not the searched for line
            if abs(angle) > 5 :
                if (i == lenLine1 - 1 and not j in DrawnLines2Ind) :
#                    plt.plot([line2[j][0][0], line2[j][1][0]], [line2[j][0][1], line2[j][1][1]], color = '#6c9f92', linestyle = 'dotted')
                    if uniqueObjectsLine2.get(j) == None :
                        length = PointDist(line2[j][0], line2[j][1])
                        uniqueObjectsLine2[j] = length
                    visited[line2[j][lineSt2tmp]] = 2
                    if not line2[j][lineSt2tmp - 1] in visited :
                        visited[line2[j][lineSt2tmp - 1]] = 1
                continue

            #if lines do exist, then choose the closest one from second file
            if (lineSt1 != -1 and lineSt2tmp != -1) :
                endDist = PointDist(line1[P1Ind][(lineSt1 + 1) % 2], line2[j][(lineSt2tmp + 1) % 2])
                #if the endpoint distance is larger than 25 then its probably not the respective line
                if (endDist >= 25 and not j in DrawnLines2Ind and i == lenLine1 - 1) :
#                    plt.plot([line2[j][0][0], line2[j][1][0]], [line2[j][0][1], line2[j][1][1]], color = '#6c9f92', linestyle = 'dotted')
                    if uniqueObjectsLine2.get(j) == None :
                        length = PointDist(line2[j][0], line2[j][1])
                        uniqueObjectsLine2[j] = length
                    visited[line2[j][lineSt2tmp]] = 2
                    if not line2[j][lineSt2tmp - 1] in visited :
                        visited[line2[j][lineSt2tmp - 1]] = 1
                    continue

                if (minEndDist == -1 and endDist < 25):
                    minEndDist = endDist
                    P2Ind = j
                    lineSt2 = lineSt2tmp
                elif endDist < minEndDist :
                    if (i == lenLine1 - 1 and not P2Ind in DrawnLines2Ind and P2Ind != -1) :
#                        plt.plot([line2[P2Ind][0][0], line2[P2Ind][1][0]], [line2[P2Ind][0][1], line2[P2Ind][1][1]], color = '#6c9f92', linestyle = 'dotted')
                        if uniqueObjectsLine2.get(P2Ind) == None :
                            length = PointDist(line2[P2Ind][0], line2[P2Ind][1])
                            uniqueObjectsLine2[P2Ind] = length
                        visited[line2[P2Ind][lineSt2]] = 2
                        if not line2[P2Ind][lineSt2 - 1] in visited :
                            visited[line2[P2Ind][lineSt2 - 1]] = 1

                    minEndDist = endDist
                    P2Ind = j
                    lineSt2 = lineSt2tmp
        print("1. var : ", visited)
        #if there is no respective line in second file then 1. file line is found only in one file -> lets mark it differently
        if (P2Ind == -1 and P2Ind != -1):
#            plt.plot([line1[P1Ind][0][0], line1[P1Ind][1][0]], [line1[P1Ind][0][1], line1[P1Ind][1][1]], color = '#055583', linestyle = 'dotted')
            if uniqueObjectsLine1.get(P1Ind) == None :
                length = PointDist(line1[P1Ind][0], line1[P2Ind][1])
                uniqueObjectsLine1[line1[P1Ind][0], line1[P2Ind][1]] = length
            if line1[P1Ind][0] == P2 :
                if not line1[P1Ind][1] in visited :
                    visited[line1[P1Ind][1]] = 1
            if line1[P1Ind][1] == P2 :
                if not line1[P1Ind][0] in visited :
                    visited[line1[P1Ind][0]] = 1
        print("2. var : ", visited)
        #if the corresponding line in second file is found, then connect respective endpoints and get the max distance between the lines, which also is the respective endpoint distance
        if (P2Ind != -1 and P1Ind != -1):
            if not line1[P1Ind][(lineSt1 + 1) % 2] in visited :
                visited[line1[P1Ind][(lineSt1 + 1) % 2]] = 1
            if not line2[P2Ind][(lineSt2 + 1) % 2] in visited :
                visited[line2[P2Ind][(lineSt2 + 1) % 2]] = 1

            DrawnLines2Ind[P2Ind] = 1
            dist1 = PointDist(line1[P1Ind][lineSt1], line2[P2Ind][lineSt2])
            dist2 = PointDist(line1[P1Ind][(lineSt1 + 1) % 2], line2[P2Ind][(lineSt2 + 1) % 2])
            if dist1 >= dist2 :
                max = dist1
            else :
                max = dist2
#            plt.plot([line1[P1Ind][lineSt1][0], line2[P2Ind][lineSt2][0]], [line1[P1Ind][lineSt1][1], line2[P2Ind][lineSt2][1]], color = '#6c9f92')
#            plt.plot([line1[P1Ind][(lineSt1 + 1) % 2][0], line2[P2Ind][(lineSt2 + 1) % 2][0]], [line1[P1Ind][(lineSt1 + 1) % 2][1], line2[P2Ind][(lineSt2 + 1) % 2][1]], color = '#6c9f92')

            #equal objects will be drawn from file 1
            print("3. var : ", visited)
            if max <= 0.1 :
                print(line1[P1Ind])
                if equalObjectsLine.get(P1Ind) == None :
                    length = PointDist(line1[P1Ind][lineSt1], line1[P1Ind][(lineSt1 + 1) % 2])
                    equalObjectsLine[P1Ind] = length
                continue

            if ObjectsLine1.get(P1Ind) == None :
                length = PointDist(line1[P1Ind][lineSt1], line1[P1Ind][(lineSt1 + 1) % 2])
                ObjectsLine1[P1Ind] = length
            if ObjectsLine2.get(P2Ind) == None :
                length = PointDist(line2[P2Ind][lineSt2], line2[P2Ind][(lineSt2 + 1) % 2])
                ObjectsLine2[P2Ind] = length
        print("4. var : ", visited)

    DrawnLines2Ind.clear()

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
                lineSt2tmp = 0
            elif (P2 == (Bezier2[j][-1][0].subs(t, 1), Bezier2[j][-1][1].subs(t, 1)) and not (Bezier2[j][-1][0].subs(t, 1), Bezier2[j][-1][1].subs(t, 1)) in visited) :
                lineSt2tmp = 1
            else:
                continue
            a2 = LineSlope((Bezier2[j][0][0].subs(t, 0), Bezier2[j][0][1].subs(t, 0)), (Bezier2[j][-1][0].subs(t, 1), Bezier2[j][-1][1].subs(t, 1)))
            angtmp = round((a2 - a1) / (1 + abs(a1 * a2)), 5)
            angle = np.arctan(angtmp)
            if abs(angle) > 5 :
                if (i == lenLine1 - 1 and not j in DrawnLines2Ind) :
#                    PlotBezier(CP2[j], '#6c9f92', 1, 'dotted')
                    visited[(Bezier2[j][-lineSt2tmp][0].subs(t, lineSt2tmp), Bezier2[j][-lineSt2tmp][1].subs(t, lineSt2tmp))] = 2
                    if not (Bezier2[j][- abs(-lineSt2tmp + 1)][0].subs(t, abs(lineSt2tmp - 1)), Bezier2[j][- abs(-lineSt2tmp + 1)][1].subs(t, abs(lineSt2tmp - 1))) :
                        visited[Bezier2[j][- abs(-lineSt2tmp + 1)][0].subs(t, abs(lineSt2tmp - 1)), Bezier2[j][- abs(-lineSt2tmp + 1)][1].subs(t, abs(lineSt2tmp))] = 1
                    if uniqueObjectsBezier2.get(j) == None :
                        length = CompositeBezierLen(CP2[j])
                        uniqueObjectsBezier2[j] = length
                continue
            if (lineSt1 != -1 and lineSt2tmp != -1) :
                endDist = PointDist((Bezier1[i][-(lineSt1 + 1) % 2][0].subs(t, (lineSt1 + 1) % 2), Bezier1[i][-(lineSt1 + 1) % 2][1].subs(t, (lineSt1 + 1) % 2)), (Bezier2[j][-(lineSt2tmp + 1) % 2][0].subs(t, (lineSt2tmp + 1) % 2),Bezier2[j][-(lineSt2tmp + 1) % 2][1].subs(t, (lineSt2tmp + 1) % 2)))
                if (i == lenLine1 - 1 and not j in DrawnLines2Ind and endDist >= 25) :
#                    PlotBezier(CP2[j], '#6c9f92', 1, 'dotted')
                    visited[(Bezier2[j][-lineSt2tmp][0].subs(t, lineSt2tmp), Bezier2[j][-lineSt2tmp][1].subs(t, lineSt2tmp))] = 2

                    if not (Bezier1[i][-(lineSt1 + 1) % 2][0].subs(t, (lineSt1 + 1) % 2), Bezier1[i][-(lineSt1 + 1) % 2][1].subs(t, (lineSt1 + 1) % 2)) in visited :
                        visited[(Bezier1[i][-(lineSt1 + 1) % 2][0].subs(t, (lineSt1 + 1) % 2), Bezier1[i][-(lineSt1 + 1) % 2][1].subs(t, (lineSt1 + 1) % 2))] = 1
                    if uniqueObjectsBezier2.get(i) == None :
                        length = CompositeBezierLen(CP2[i])
                        uniqueObjectsBezier2[i] = length
                if (minEndDist == -1 and endDist < 25):
                    minEndDist = endDist
                    P2Ind = j
                    lineSt2 = lineSt2tmp
                if endDist < minEndDist :
                    if (i == lenLine1 - 1 and not P2Ind in DrawnLines2Ind) :
#                        PlotBezier(CP2[P2Ind], '#6c9f92', 1, 'dotted', linewidth = 4)
                        visited[(Bezier2[P2Ind][-lineSt2][0].subs(t, lineSt2), Bezier2[P2Ind][-lineSt2][1].subs(t, lineSt2))] = 2
                        if not (Bezier2[P2Ind][-abs(-lineSt2 + 1)][0].subs(t, abs(lineSt2 - 1)), Bezier2[P2Ind][-abs(-lineSt2 + 1)][1].subs(t, abs(lineSt2 - 1))) :
                            visited[(Bezier2[P2Ind][-abs(-lineSt2 + 1)][0].subs(t, abs(lineSt2 - 1)), Bezier2[P2Ind][-abs(-lineSt2 + 1)][1].subs(t, abs(lineSt2 - 1)))] = 1
                        if uniqueObjectsBezier2.get(P2Ind) == None :
                            length = CompositeBezierLen(CP2[P2Ind])
                            uniqueObjectsBezier2[P2Ind] = length
                    minEndDist = endDist
                    P2Ind = j
                    lineSt2 = lineSt2tmp

        DrawnLines2Ind.clear()

        if (P2Ind == -1) :
#            PlotBezier(CP1[i], '#055583', 1, 'dotted')
            t = Symbol('t')
            visited[(Bezier1[i][-lineSt1][0].subs(t, lineSt1), Bezier1[i][-lineSt1][1].subs(t, lineSt1))] = 2
            if not (Bezier1[i][-abs(-lineSt1 + 1)][0].subs(t, abs(lineSt1 - 1)), Bezier1[i][-abs(-lineSt1 + 1)][1].subs(t, abs(lineSt1 - 1))) :
                visited[(Bezier1[i][-abs(-lineSt1 + 1)][0].subs(t, abs(lineSt1 - 1)), Bezier1[i][-abs(-lineSt1 + 1)][1].subs(t, abs(lineSt1 - 1)))] = 1
            if uniqueObjectsBezier1.get(i) == None :
                length = CompositeBezierLen(CP1[i])
                uniqueObjectsLine1[i] = length
        if (P2Ind != -1 and P1Ind != -1):
            if not (Bezier1[i][-abs(-lineSt1 + 1)][0].subs(t, abs(lineSt1 - 1)), Bezier1[i][-abs(-lineSt1 + 1)][1].subs(t, abs(lineSt1 - 1))) in visited :
                visited[Bezier1[i][-abs(-lineSt1 + 1)][0].subs(t, abs(lineSt1 - 1)), Bezier1[i][-abs(-lineSt1 + 1)][1].subs(t, abs(lineSt1 - 1))] = 1
            if not (Bezier2[P2Ind][-abs(-lineSt2 + 1)][0].subs(t, abs(lineSt2 - 1)), Bezier2[P2Ind][-abs(-lineSt2 + 1)][1].subs(t, abs(lineSt2 - 1))) :
                visited[Bezier2[P2Ind][-abs(-lineSt2 + 1)][0].subs(t, abs(lineSt2 - 1)), Bezier2[P2Ind][-abs(-lineSt2 + 1)][1].subs(t, abs(lineSt2 - 1))] = 1

            max = MaxBezierDist(Bezier1[i], Bezier2[i])
            #equal objects will be drawn from file 1
            if max <= 1 :
                if equalObjectsBezier.get(i) == None :
                    length = CompositeBezierLen(CP1[i])
                    equalObjectsBezier[i] = length
                continue

            if ObjectsBezier1.get(i) == None :
                length = CompositeBezierLen(CP1[i])
                ObjectsBezier1[i] = length
            if ObjectsBezier2.get(P2Ind) == None :
                length = CompositeBezierLen(CP2[P2Ind])
                ObjectsBezier2[P2Ind] = length

            len1 = CompositeBezierLen(CP1[i])
            len2 = CompositeBezierLen(CP2[P2Ind])

            int1 = len1 / 10
            int2 = len2 / 10
            LeastSquare(Bezier1[i], int1, lineSt1, Bezier2[P2Ind], int2, lineSt2, CP1[i], CP2[P2Ind])

            MaxBezierDist(Bezier1[i], Bezier2[P2Ind])
#            BezierDiff(Bezier1[i], int1, lineSt1, Bezier2[P2Ind], int2, lineSt2, CP1, CP2)
    visited[P2] = 1
    visited[P1] = 2
    print(P1)
    print("5. var : ", visited)
    return visited

def CompositeBezierLen(CP) :
    length = 0
    lenCP = len(CP)
    for i in range(lenCP) :
        length += CubicBezierLen(CP[i])
    return length



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
print("line1")
for i in line1 :
    if not i[0] in visited :
        visited = DiffAll(i[0], visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
    elif not visited[i[0]] == 1 :
        visited = DiffAll(i[0], visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
print(visited)
print("line2")
for i in line2 :
    if not i[0] in visited and not i[1] in visited :
        visited = OneFileObject(i[0], line2, Bezier2, CP2, visited, 2)

t = Symbol('t')
for i in Bezier1 :
    if not (i[0][0].subs(t, 0), i[0][1].subs(t, 0)) in visited :
        visited = DiffAll((i[0][0].subs(t, 0), i[0][1].subs(t, 0)), visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
for i in Bezier2 :
    if not (i[0][0].subs(t, 0), i[0][1].subs(t, 0)) in visited and not (i[-1][0].subs(t, 1), i[-1][1].subs(t, 1)) in visited :
        visited = OneFileObject((i[0][0].subs(t, 0), i[0][1].subs(t, 0)), line2, Bezier2, CP2, visited, 2)

def CheckWhereVisited1(P, visited, line1, line2, Bezier1, Bezier2, CP1, CP2) :
    return visited

for i in line1 :
    if i[0] in visited :
        if visited[i[0]] == 1 :
            visited = CheckWhereVisited1(i[0], visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
    if i[1] in visited :
        if visited[i[1]] == 1 :
            visited = CheckWhereVisited1(i[0], visited, line1, line2, Bezier1, Bezier2, CP1, CP2)
print(visited)
visited.clear()
def ClosestCircle(CP, r, circle2) :
    for i in circle :
        dist = PointDist(CP, i)
        if dist < 0.1:
            if (r + 0.05 > circle2[i] and circle2[i] > r - 0.05) :
                closest = i
#find differences between circles in bothe files
def CircleDiff(CP, r, circle1, circle2, visited) :
    P2 = ClosestCircle(CP, r, circle2)

for i in circle1 :
    if not centerPoints1 in visited :
        visited = CircleDiff(i, circle1, circle2, visited)

print("equal Line")
print(equalObjectsLine)

print("equal Bezier")
print(equalObjectsBezier)

print("line1")
print(ObjectsLine1)

print("line2")
print(ObjectsLine2)

print("UniqueLine1")
print(uniqueObjectsLine1)

print("Unique Line2")
print(uniqueObjectsLine2)


print("Bezier1")
print(ObjectsBezier1)

print("Bezier2")
print(ObjectsBezier2)

print("UniqueBezier1")
print(uniqueObjectsBezier1)

print("Unique Bezier2")
print(uniqueObjectsBezier2)




nreqf = 0.0
nr1f = 1.0
nr2f = 2.0
nrf = 3.0

lenLine1 = len(line1)
lenLine2 = len(line2)
lenBezier1 = len(Bezier1)
lenBezier2 = len(Bezier2)
#ax.annotate(nr1f, xy = (200, 200), xytext = (200,200) )


for i in range(lenLine1) :
    if i in equalObjectsLine :
        nreqf += 0.01
        plt.plot([line1[i][0][0], line1[i][1][0]], [line1[i][0][1], line1[i][1][1]],
            color = '#717577', alpha = 0.65)
        a = ((line1[i][0][0] + line1[i][1][0]) / 2)
        b = ((line1[i][0][1] + line1[i][1][1]) / 2)
        ax.annotate(round(nreqf, 2), xy = (a, b), xytext = (a, b), family = 'monospace',
            color = '#3d4042', ha = 'center', va = 'center', weight = 'bold')
        print(round(nreqf, 2), ' - ', equalObjectsLine[i] / 10)
    if i in ObjectsLine1 :
        nr1f += 0.01
        plt.plot([line1[i][0][0], line1[i][1][0]], [line1[i][0][1], line1[i][1][1]],
            color = firstFileCol, alpha = 0.65)
        a = ((line1[i][0][0] + line1[i][1][0]) / 2)
        b = ((line1[i][0][1] + line1[i][1][1]) / 2)
        ax.annotate(round(nr1f, 2), xy = (a, b), xytext = (a, b), family = 'monospace',
            color = '#00304c', ha = 'center', va = 'center', weight = 'bold')
        print(round(nr1f, 2), ' - ', ObjectsLine1[i] / 10)
    if i in uniqueObjectsLine1 :
        nr1f += 0.01
        plt.plot([line1[i][0][0], line1[i][1][0]], [line1[i][0][1], line1[i][1][1]],
            color = firstFileCol, alpha = 0.65, linestyle = 'dotted')
        a = ((line1[i][0][0] + line1[i][1][0]) / 2)
        b = ((line1[i][0][1] + line1[i][1][1]) / 2)
        ax.annotate(round(nr1f, 2), xy = (a, b), xytext = (a, b), family = 'monospace',
            color = '#00304c', ha = 'center', va = 'center', weight = 'bold')
        print(round(nr1f, 2), ' - ', uniqueObjectsLine1[i] / 10)

print()
for i in range(lenLine2) :
    if i in ObjectsLine2 :
        nreqf += 0.01
        plt.plot([line2[i][0][0], line2[i][1][0]], [line2[i][0][1], line2[i][1][1]],
            color = secondFileCol, alpha = 0.5)
        a = ((line2[i][0][0] + line2[i][1][0]) / 2)
        b = ((line2[i][0][1] + line2[i][1][1]) / 2)
        ax.annotate(round(nreqf, 2), xy = (a, b), xytext = (a, b), family = 'monospace',
            color = '#59080a', ha = 'center', va = 'center', weight = 'bold')
        print(round(nreqf, 2), ' - ', ObjectsLine2[i] / 10 )
    if i in uniqueObjectsLine2 :
        nr2f += 0.01
        plt.plot([line2[i][0][0], line2[i][1][0]], [line2[i][0][1], line2[i][1][1]],
            color = secondFileCol, alpha = 0.5, linestyle = 'dotted')
        a = ((line2[i][0][0] + line2[i][1][0]) / 2)
        b = ((line2[i][0][1] + line2[i][1][1]) / 2)
        ax.annotate(round(nr2f, 2), xy = (a, b), xytext = (a, b), family = 'monospace',
            color = '#59080a', ha = 'center', va = 'center', weight = 'bold')
        print(round(nr2f, 2), ' - ', uniqueObjectsLine2[i] / 10 )

print()
lenBezier1 = len(Bezier1)
lenBezier2 = len(Bezier2)

for i in range(lenBezier1) :
    if i in equalObjectsBezier :
        PlotBezier(CP1[i], '#717577', 0.65, 'solid')
        nr1f += 0.01
        pos = int(round(len(Bezier1[i]) / 3))
        a = Bezier1[i][pos][0].subs(t, 0)
        b = Bezier1[i][pos][1].subs(t, 0)
        ax.annotate(round(nreqf, 2), xy = (a, b), xytext = (a, b), family = 'monospace',
            color = '#3d4042', ha = 'center', va = 'center', weight = 'bold')
        print(round(nreqf, 2), ' - ', equalObjectsBezier[i] / 10)
    if i in ObjectsBezier1 :
        PlotBezier(CP1[i], firstFileCol, 0.65, 'solid')
        nr1f += 0.01
        pos = int(round(len(Bezier1[i]) / 3))
        a = Bezier1[i][pos][0].subs(t, 0)
        b = Bezier1[i][pos][1].subs(t, 0)
        ax.annotate(round(nr1f, 2), xy = (a, b), xytext = (a, b), family = 'monospace',
            color = '#00304c', ha = 'center', va = 'center', weight = 'bold')
        print(round(nr1f, 2), ' - ', ObjectsBezier1[i] / 10)
    if i in uniqueObjectsBezier1 :
        PlotBezier(CP1[i], firstFileCol, 0.65, 'dotted')
        nr1f += 0.01
        pos = int(round(len(Bezier1[i]) / 3))
        a = Bezier1[i][pos][0].subs(t, 0)
        b = Bezier1[i][pos][1].subs(t, 0)
        ax.annotate(round(nr1f, 2), xy = (a, b), xytext = (a, b), family = 'monospace',
            color = '#00304c', ha = 'center', va = 'center', weight = 'bold')
        print(round(nr1f, 2), ' - ', uniqueObjectsBezier1[i] / 10)

print()
for i in range(lenBezier2) :
    if i in ObjectsBezier2 :
        PlotBezier(CP2[i], secondFileCol, 0.5, 'solid')
        nr2f += 0.01
        pos = int(round(len(Bezier2[i]) / 4))
        a = Bezier2[i][pos][0].subs(t, 0)
        b = Bezier2[i][pos][1].subs(t, 0)
        ax.annotate(round(nr2f, 2), xy = (a, b), xytext = (a, b), family = 'monospace',
            color = '#59080a', ha = 'center', va = 'center', weight = 'bold')
        print(round(nr2f, 2), ' - ', ObjectsBezier2[i] / 10)
    if i in uniqueObjectsBezier2 :
        PlotBezier(CP2[i], secondFileCol, 0.5, 'dotted')
        nr2f += 0.01
        pos = int(round(len(Bezier2[i]) / 4))
        a = Bezier2[i][pos][0].subs(t, 0)
        b = Bezier2[i][pos][1].subs(t, 0)
        ax.annotate(round(nr2f, 2), xy = (a, b), xytext = (a, b), family = 'monospace',
            color = '#59080a', ha = 'center', va = 'center', weight = 'bold')
        print(round(nr2f, 2), ' - ', uniqueObjectsBezier2[i] / 10)

plt.show()
