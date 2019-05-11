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
    Bx, By = BezierFormulaComp(C)
    Bezier.append([Bx, By])
    for i in range(1, len(P)): #Calculate controlpoints for P[i], P[i + 1] segment
        dtmp = PointDist(P[i], P[i + 1]) / 3
        d = dtmp
        C[2] = DistantPoint(P[i + 1], P[i + 2], P[i], d)
        C[0] = P[i]
        C[3] = P[i + 1]
        C[1] = DistantPoint(P[i], P[i - 1], P[i + 1], d)
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
            Bx, By = BezierFormulaComp(C)
            Bezier.append([Bx, By])
            break
    return Bezier

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
        Bxtmp = (UnevaluatedExpr(BinC))*(1 - t)**(4 - i - 1)*t**(i)*(C[i][0])    #Create Bx(t) formula
        Bx = UnevaluatedExpr(Bx) + Bxtmp
        Bytmp = (UnevaluatedExpr(BinC))*(1 - t)**(4 - i - 1)*t**(i)*(C[i][1])  #Create Bx(t) formula
        By = UnevaluatedExpr(By) + Bytmp
#        Bx = f"{Bx} + ({BinC})(1 - t)^{n - i}*t^{i - 1}({P[i - 1][0]}) " #Create Bx(t) formula as a string
#        By = f"{By} + ({BinC})(1 - t)^{n - i}*t^{i - 1}({P[i - 1][1]}) " #Create Bx(t) formula as a string
    return (Bx, By)

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
    return np.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)

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






dxf1 = dxfgrabber.readfile("svarki_002.dxf")

#type of objects in file
print("type of objects in file")
type1 = [entity.dxftype for entity in dxf1.entities]
print(type1)

output1 = [entity for entity in dxf1.entities]

Color1 = [entity.color for entity in output1]
print(Color1)
Linetype1 = [entity.linetype for entity in output1 if entity.dxftype != 'POINT']
print(Linetype1)

point1 = []
line1 = []
Bezier1 = []

#get parameters of objects
for entity in output1:
    #Point
    if entity.dxftype == 'POINT':
        pointCoord1 = entity.point
        point1.append([pointCoord1[0], pointCoord1[1]])
    #Line
    if entity.dxftype == 'LINE':
        lineStart1 = entity.start
        lineEnd1 = entity.end
        LFx, LFy = ParamLineFormula(lineStart1, lineEnd1)
        line1.append([LFx, LFy])

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

#        print(" POLYLINE\n")
#        print("is closed")
#        print(PolylineIsClosed1)
#        print("Spline type")
#        print(PolylineSplineType1)
#        print("line points")
#        print(PolylinePoints1)
#        print("Control points")
#        print(PolylineControlPoints)
#        print("bulge")
#        print(PolylineBulge)
        print("vertex count:")
        print(PolylineVertexCount1)
        PointCnt1 = PolylineVertexCount1
        Bezier1.append(CompositeBezier(min(8, PointCnt1), PolylinePoints1))


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
print("=============")

#Bezier formulas of file #1


dxf2 = dxfgrabber.readfile("aplis_002.dxf")
#type of objects in file
print("type of objects in file")
type2 = [entity.dxftype for entity in dxf2.entities]
print(type2)

output2 = [entity for entity in dxf2.entities]

Color2 = [entity.color for entity in output2]
print(Color2)
Linetype2 = [entity.linetype for entity in output2 if entity.dxftype != 'POINT']
print(Linetype2)

point2 = []
line2 = []
Bezier2 = []
#get parameters of objects
for entity in output2:
    #Point
    if entity.dxftype == 'POINT':
        pointCoord2 = entity.point
        point2.append([pointCoord2[0], pointCoord2[1]])
    #Line
    if entity.dxftype == 'LINE':
        lineStart2 = entity.start
        lineEnd2 = entity.end
        LFx, LFy = ParamLineFormula(lineStart2, lineEnd2)
        line2.append([LFx, LFy])

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

        print(" POLYLINE\n")
        print("is closed")
        print(PolylineIsClosed2)
        print("Spline type")
        print(PolylineSplineType2)
        print("line points")
        print(PolylinePoints2)
#        print("Control points")
#        print(PolylineControlPoints)
#        print("bulge")
#        print(PolylineBulge)
        print("vertex count:")
        print(PolylineVertexCount2)
        PointCnt2 = PolylineVertexCount2
        Bx2, By2 = GetBezier(5, min(8, PointCnt2), PolylinePoints2)
        print("formula")
#        print(latex(Bx))
#        print(latex(By))
        Bezier2.append(CompositeBezier(min(8, PointCnt2), PolylinePoints2))

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



#We have got data from both given files

#File #1
#for i in point1:
#    print('(', i[0], ',', i[1], ')')
#for i in line1:
#    print('(',latex(i[0]),',', latex(i[1]),')')
#for i in Bezier1 :
#    for j in range(len(i)):
#    #    print('(',latex(i[j][0]),',', latex(i[j][1]),')')
#        a = 1




#File #2
#for i in point2:
#    print('(', i[0], ',', i[1], ')')
#for i in line2:
#    print('(',latex(i[0]),',', latex(i[1]),')')
#for i in Bezier2 :
#    for j in range(len(i)):
#        #    print('(',latex(i[j][0]),',', latex(i[j][1]),')')
#        a = 1



#Lets find which points are represented as the 'same' points in both files

#First we will take 3 points from first file and get the distances between them
#Based on those 3 point we are going to take three point sets from second file and compare if they match based on proportions
#(from those three points we save proportion and angle and in the next steps we are going to use vectors to find the corresponding points)
#If we have found the match then - take one point, pair it with all the other points and find their matching pairs in second file
#Because we do that with vector there is no reason to compare each two points

for i in points1 :    
