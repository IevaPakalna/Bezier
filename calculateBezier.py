from sympy import init_printing, Symbol, UnevaluatedExpr, expand, pretty_print as pprint, latex
import numpy
import matplotlib.pyplot as plt
from numpy.linalg import inv

import calculations

def __init__(Bezier) :
    pass
#Crete Composite Bezier from given points
def CompositeBezier(Points, nr):
    #composite bezier segments controlpoins are placed respectively:
    #1. calculate bisector between two sequent lines crosspoint
    #2. calculate perpendicular for this bisector at its "start point"
    #3. set controlpoint in distance 1/3 of respective sequent point distance
    P = PMtrx(len(Points), Points)
    C = [[0],[1],[2],[3]]
    CP = []
    Bezier = []
    t = Symbol('t')
    x = Symbol('x')
    dtmp = calculations.PointDist(P[0], P[1])
    d = dtmp
    if P[0][0] < P[2][0] :
        d = - dtmp
    bP = calculations.Bisector(P[0], P[1], P[2])
    perpF = calculations.PerpFormula(bP, P[1], P[1])
    C[2] = calculations.DistantPoint(P[1], P[1], (-1, perpF.subs(x, -1)), d / 3) #Point on (P[1], P[2]) line in 1/3 distance of [P[0],P[1]] (C2 controlpoint)
    if P[0][0] > P[1][0] :
        d = - dtmp
    else :
        d = dtmp
    middlePnt = calculations.DistantPoint(P[0], P[0], P[1], d / 2) #Middle point of [P[0],P[1]]
    PF1 = calculations.PerpFormula(P[0], P[1], middlePnt)   #Middle perpendicular of [P[0],P[1]]
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
    PF2 = calculations.PerpFormula(middlePnt, PF1tmp, C[2]) #Perpendicular of PF1 that goes trough C[2]
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
    PF1xPF2 = calculations.LineIntersect(middlePnt, PF1tmp, C[2], PF2tmp)   #Intersection point of PF1 and PF2
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
    lenP = len(P)
    for i in range(1, lenP + 1): #Calculate controlpoints for P[i], P[i + 1] segment
        if lenP > 3 and i <= lenP - 3 :
            dtmp = calculations.PointDist(P[i], P[i + 1]) / 3
            d = dtmp
            if P[i][0] < P[i + 2][0] :
                d = - dtmp
            bP = calculations.Bisector(P[i], P[i + 1], P[i + 2])
            perpF = calculations.PerpFormula(bP, P[i + 1], P[i + 1])
            C[2] = calculations.DistantPoint(P[i + 1], P[i + 1], (-1, perpF.subs(x, -1)), d)
            C[0] = P[i]
            C[3] = P[i + 1]
            if P[i - 1] < P[i + 1] :
                d = dtmp
            else :
                d = - dtmp
            bP = calculations.Bisector(P[i - 1], P[i], P[i + 1])
            perpF = calculations.PerpFormula(bP, P[i], P[i])
            C[1] = calculations.DistantPoint(P[i], P[i], (-1, perpF.subs(x, -1)), d)
            CP.append([C[0], C[1], C[2], C[3]])
            Bx, By = BezierFormulaComp(C)
            Bezier.append([Bx, By])
        elif i == len(P) - 2:
            C[0] = P[i]

            dtmp = calculations.PointDist(P[i + 1], P[i]) / 3
            d = dtmp
            if P[i][0] > P[i + 1][0] :
                d = - dtmp
            bP = calculations.Bisector(P[i - 1], P[i], P[i + 1])
            perpF = calculations.PerpFormula(bP, P[i], P[i])
            C[1] = calculations.DistantPoint(P[i], P[i], (-1, perpF.subs(x, -1)), d)
            Ptmp1 = calculations.DistantPoint(P[i], P[i], P[i + 1], d)
            perpF = calculations.PerpFormula(P[i + 1], P[i], Ptmp1)
            dtmp = calculations.PointDist(P[i + 1], P[i]) / 2
            if P[i][0] > P[i + 1][0] :
                d = - dtmp
            else :
                d = dtmp
            Ptmp = calculations.DistantPoint(P[i], P[i], P[i + 1], d)
            perptmpF = calculations.PerpFormula(C[1], Ptmp1, C[1])
            Ptmp2 = calculations.LineIntersect(C[1], (-1, perptmpF.subs(x, -1)), Ptmp, (-1, perpF.subs(x, -1)))
            C[2] = [2 * Ptmp2[0] - C[1][0], 2 * Ptmp2[1] - C[1][1]]
            C[3] = P[i + 1]

            CP.append([C[0], C[1], C[2], C[3]])
            Bx, By = BezierFormulaComp(C)
            Bezier.append([Bx, By])
            break
    return Bezier, CP

#returns Bezier curves formula (for given controlpoints)
def BezierFormulaComp(C):
    t = Symbol('t')
    Bx = 0
    By = 0
    for i in range(4):
        BinC = calculations.BinCoef(3, i)
        Bxtmp = BinC*(1 - t)**(4 - i - 1)*t**(i)*(C[i][0])    #Create Bx(t) formula
        Bx = Bxtmp + Bx
        Bytmp = BinC*(1 - t)**(4 - i - 1)*t**(i)*(C[i][1])  #Create Bx(t) formula
        By = Bytmp + By
#        Bx = f"{Bx} + ({BinC})(1 - t)^{n - i}*t^{i - 1}({P[i - 1][0]}) " #Create Bx(t) formula as a string
#        By = f"{By} + ({BinC})(1 - t)^{n - i}*t^{i - 1}({P[i - 1][1]}) " #Create Bx(t) formula as a string
    if Bx == 0 :
        Bx == 0 * t
    if By == 0 :
        By == 0 * t
    return Bx, By

#Plots Bezier curve
def PlotBezier(C, clr, transp, lstyle) :
    for i in range(len(C)) :
        t1 = numpy.linspace(0, 1, 50)
        Bx1 = 0
        By1 = 0
        for j in range(4) :
            BinC = calculations.BinCoef(3, j)
            Bxtmp = BinC * (1 - t1)**(j) * t1**(4 - j - 1) * (C[i][3 - j][0])    #Create Bx(t) formula
            Bx1 = (Bx1) + Bxtmp
            Bytmp = BinC * (1 - t1)**(j) * t1**(4 - j - 1) * (C[i][3 - j][1])  #Create Bx(t) formula
            By1 = (By1) + Bytmp
        plt.plot(Bx1, By1, color = clr, alpha = transp, linestyle = lstyle)
    return

#Returns parametric line Formula
def ParamLineFormula(P1, P2):
    t = Symbol('t')
    PLFx = ((1 - t) * P1[0] + t * P2[0])
    PLFy = ((1 - t) * P1[1] + t * P2[1])
    return PLFx, PLFy

#Returns line formula
def LineFormula(P1, P2):
    x = Symbol('x')
    if P1[0] == P2[0]:
        LF = (P1[0] + 0 * x)
        return LF
    LF = (UnevaluatedExpr((P1[1] - P2[1]) / (P1[0] - P2[0])) * x + (P1[0] * P2[1] - P2[0] * P1[1]) / (P1[0] - P2[0]))
    return LF

#Returns Bezier formula through given points
def GetBezier(n, Pn, Points):
    P = PMtrx(Pn, Points)
    Pn = len(P)
    T = TMtrx(n, Pn, P)
    M = MMtrx(n)
    C = CMtrx(Pn, M, P, T)
    return BezierFormula(n, C)

#Bezier formula for n-th degree curve
def BezierFormula(n, C):
    Bx = 0
    By = 0
    t = Symbol('t')
    for i in range(1, n + 1):
        BinC = calculations.BinCoef(n - 1, i - 1)
        Bxtmp = (UnevaluatedExpr(BinC))*(1 - t)**(n - i)*t**(i - 1)*(C[i - 1][0])    #Create Bx(t) formula
        Bx = Bx + Bxtmp
        Bytmp = (UnevaluatedExpr(BinC))*(1 - t)**(n - i)*t**(i - 1)*(C[i - 1][1])  #Create Bx(t) formula
        By = By + Bytmp
#        Bx = f"{Bx} + ({BinC})(1 - t)^{n - i}*t^{i - 1}({P[i - 1][0]}) " #Create Bx(t) formula as a string
#        By = f"{By} + ({BinC})(1 - t)^{n - i}*t^{i - 1}({P[i - 1][1]}) " #Create Bx(t) formula as a string
    return (Bx, By)

#https://pomax.github.io/bezierinfo/#curvefitting
def MMtrx(n):
    t = Symbol('t')
    Mtrx = []
    for i in range(n):
        Mtrx.append([])
        BinC = calculations.BinCoef(n - 1 , i)
        Mtmp = expand((BinC)*(1 - t)**(n - i - 1)*t**(i))
        for j in range(n):
            Mtrx[i].append(float(Mtmp.coeff(t, n - j - 1)))
    return(Mtrx)
#Get fit-point matrix
def PMtrx(n, Points) :
    P = []
    for i in range(n):
        P.append(Points[i])
    return P
#get T matrix with parameter values
def TMtrx(n, Pn, P):
    d = []
    d.append(0)
    for i in range(1, Pn):   #create point (t) distance array
        dist = d[i - 1] + calculations.PointDist(P[i-1], P[i])
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
    M = numpy.flip(M, 0)
    Tt = numpy.transpose(T)
    Mi = inv(M)
    C = numpy.matmul(Tt, T)
    C = inv(C)
    C = numpy.matmul(Mi, C)
    C = numpy.matmul(C, Tt)
    C = numpy.matmul(C, P)
    return C
