import numpy
from sympy import init_printing, Symbol, UnevaluatedExpr, expand, pretty_print as pprint, latex





def __init__(calc) :
    pass

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
    a = numpy.array([[S1, y1],[S2, y2]], dtype = 'float')
    b = numpy.array([-b1, -b2], dtype = 'float')
    Ptmp = numpy.linalg.solve(a,b)
    P = []
    P.append(Ptmp[0])
    P.append(Ptmp[1])
    return P
#Returns slope of line trough two given points
def LineSlope(P1, P2):
    if P1[0] == P2[0]:
        return 1000
    return (P1[1] - P2[1]) / (P1[0] - P2[0])
#get vector of given line segment (P1 - start point, P2 - end point)
def GetVector(P1, P2) :
    vx = P2[0] - P1[0]
    vy = P2[1] - P1[1]
    return [vx, vy]
#returns angle between two given vectors
def GetAngle(a, b) :
    alen = PointDist((0, 0), a)
    blen = PointDist((0, 0), b)
    alpha = numpy.arccos(round((a[0] * b[0] + a[1] * b[1]) / (alen * blen), 5))
    return alpha
#Rotate point P around point rP by given angle alpha
def Rotate(P, rP, alpha):
    x = numpy.cos(alpha) * (P[0] - rP[0]) - numpy.sin(alpha) * (P[1] - rP[1]) + rP[0]
    y = numpy.sin(alpha) * (P[0] - rP[0]) + numpy.cos(alpha) * (P[1] - rP[1]) + rP[1]
    return (x, y)
#returns point that lies on bisector of angle given by three points of which the vertex is the middle one
def Bisector(P1, P2, P3) :
    a = GetVector(P2, P1)
    b = GetVector(P2, P3)
    alpha = GetAngle(a, b) / 2
    nP = Rotate(P1, P2, -alpha)
    return nP
#Returns distance between two points
def PointDist(a, b):
    dist = numpy.sqrt(round((a[0] - b[0])**2 + (a[1] - b[1])**2, 10))
    return dist
#Calculate distant point along a line a certain distance away given point
#includes transportation using paralel line, if the given point doesn`t lie on line
def DistantPoint(P, P1, P2, d):
    Ptmp = []
    dist = 0
    v = []
    Pn = []
    P1n = []
    P2n = []
    a = LineSlope(P1, P2)
    alpha = numpy.arctan(round(a / (1), 5))
    Pn.append(numpy.cos(alpha) * ((P[0] + d) - P[0]) - numpy.sin(alpha) * (P[1] - P[1])  + P[0])
    Pn.append(numpy.sin(alpha) * ((P[0] + d) - P[0]) + numpy.cos(alpha) * (P[1] - P[1])  + P[1])
    return(Pn)
#Get Perpendicular formula
def PerpFormula(P1, P2, P):
    x = Symbol('x')
    if P1[1] == P2[1]:
        PF = P[1] + 0 * x
        return PF
    PF = - (P1[0] - P2[0]) / (P1[1] - P2[1]) * x + P[0] * ((P1[0] - P2[0]) / (P1[1] - P2[1])) + P[1]
    return(PF)
#get factorial
def Fact(n):
    f = 1
    for i in range(n):
        f = f * (i + 1)
    return f
#get binominal coefficients
def BinCoef(n, k):
    return int(Fact(n)/(Fact(k) * Fact(n - k)))
#transform points
def pointTransform(P, Vx, Vy, rP, alpha) :
    for i in range(len(P) - 1) :
        tmpx = (P[i][0] - Vx)
        tmpy = (P[i][1] - Vy)
        tmpx = numpy.cos(alpha) * (tmpx - rP[0]) - numpy.sin(alpha) * (tmpy - rP[1]) + rP[0]
        tmpy = numpy.sin(alpha) * (tmpx - rP[0]) + numpy.cos(alpha) * (tmpy - rP[1]) + rP[1]
        P[i][0] = tmpx
        P[i][1] = tmpy
    return P
