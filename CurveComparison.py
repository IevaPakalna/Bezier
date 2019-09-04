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
import subprocess
import os
import sys
import logging
import fileinput


maxOffset = 0.4 #cm maximal acceptable similar object offset



filename1 = sys.argv[1]
filename2 = sys.argv[2]
showPlot = sys.argv[3].replace('showPlot=', '')
if showPlot == "True" :
    showPlot = True
else:
    showPlot = False

statFile = sys.argv[4]
if not os.path.exists(statFile):
    stat = open(statFile, 'w')
    stat.close()

handler = logging.FileHandler(statFile.replace('.log', '_INFO.log'))
handler.setFormatter(logging.Formatter('%(message)s'))
statInfolog = logging.getLogger("statInfolog")
statInfolog.setLevel(logging.INFO)
statInfolog.addHandler(handler)

maxOffset = float(sys.argv[5].replace('maxOffset(cm)=', ''))



#logging.basicConfig(filename = "info.log", format = '%(message)s', level = logging.INFO)







##DXF
#fileName1 = sys.argv[1]
##fileName1 = str(input("Enter DXF file path:   "))
##fileName1 = 'D:/prog/svgFiles/torss/SIEVIETES_TORSS__AR__KOEF.svgSheet.ii.svg'
##SVG
#fileName2 = sys.argv[2]
##fileName2 = str(input("Enter SVG file path:   "))
##fileName2 = 'D:/prog/svgFiles/torss/SIEVIETES_TORSS__AR__KOEF.svgSheet.iii.svg'
#
##subprocess.run(["c:\Program Files\Inkscape\inkscape.com", "-f", fileName2, "-E", "D:/prog/file2.eps","--export-ignore-filters", "--export-ps-level=3"])
##subprocess.run(["c:\Program Files\pstoedit\pstoedit.exe", "-f", "dxf:-polyaslines", "D:/prog/file2.eps", "D:/prog/file2.dxf"])
#
##fileName2 = 'D:/prog/file2.dxf'
#
##Filenames should be written in the following two lines, as shown in example:
##   dxf1(2) = dxfgrabber.readfile("filename")
##dxf1 = dxfgrabber.readfile(fileName1)
##dxf2 = dxfgrabber.readfile(fileName2)


firstFileCol = '#000ec7'
secondFileCol = '#bc0e13'



class Calculations :
    def __init__(calc) :
        pass

    #Coordinates of projected point on a line
    def PointProjecOnLine(P, P1, P2):
        slope = Calculations.LineSlope(P1, P2)
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
            S1 = Calculations.LineSlope(P11, P12)
        if P21[1] == P22[1] :
            S2 = 0
        else:
            S2 = Calculations.LineSlope(P21, P22)
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
            return 1000
        return (P1[1] - P2[1]) / (P1[0] - P2[0])
    #get vector of given line segment (P1 - start point, P2 - end point)
    def GetVector(P1, P2) :
        vx = P2[0] - P1[0]
        vy = P2[1] - P1[1]
        return [vx, vy]
    #returns angle between two given vectors
    def GetAngle(a, b) :
        alen = Calculations.PointDist((0, 0), a)
        blen = Calculations.PointDist((0, 0), b)
        alpha = np.arccos((a[0] * b[0] + a[1] * b[1]) / (alen * blen))
        return alpha
    #Rotate point P around point rP by given angle alpha
    def Rotate(P, rP, alpha):
        x = np.cos(alpha) * (P[0] - rP[0]) - np.sin(alpha) * (P[1] - rP[1]) + rP[0]
        y = np.sin(alpha) * (P[0] - rP[0]) + np.cos(alpha) * (P[1] - rP[1]) + rP[1]
        return (x, y)
    #returns point that lies on bisector of angle given by three points of which the vertex is the middle one
    def Bisector(P1, P2, P3) :
        a = Calculations.GetVector(P2, P1)
        b = Calculations.GetVector(P2, P3)
        alpha = Calculations.GetAngle(a, b) / 2
        nP = Calculations.Rotate(P1, P2, -alpha)
        return nP
    #Returns distance between two points
    def PointDist(a, b):
        dist = np.sqrt(round((a[0] - b[0])**2 + (a[1] - b[1])**2, 10))
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
        a = Calculations.LineSlope(P1, P2)
        alpha = np.arctan(round(a / (1), 5))
        Pn.append(np.cos(alpha) * ((P[0] + d) - P[0]) - np.sin(alpha) * (P[1] - P[1])  + P[0])
        Pn.append(np.sin(alpha) * ((P[0] + d) - P[0]) + np.cos(alpha) * (P[1] - P[1])  + P[1])
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
        return int(Calculations.Fact(n)/(Calculations.Fact(k) * Calculations.Fact(n - k)))
    #transform points
    def pointTransform(P, Vx, Vy, rP, alpha) :
        for i in range(len(P) - 1) :
            tmpx = (P[i][0] - Vx)
            tmpy = (P[i][1] - Vy)
            tmpx = np.cos(alpha) * (tmpx - rP[0]) - np.sin(alpha) * (tmpy - rP[1]) + rP[0]
            tmpy = np.sin(alpha) * (tmpx - rP[0]) + np.cos(alpha) * (tmpy - rP[1]) + rP[1]
            P[i][0] = tmpx
            P[i][1] = tmpy
        return P

class CalculateBezier(Calculations) :
    def __init__(Bezier) :
        pass
    #Crete Composite Bezier from given points
    def CompositeBezier(Points, nr):
        P = CalculateBezier.PMtrx(len(Points), Points)
        C = [[0],[1],[2],[3]]
        CP = []
        Bezier = []
        t = Symbol('t')
        x = Symbol('x')
        dtmp = Calculations.PointDist(P[0], P[1])
        d = dtmp
        if P[0][0] < P[2][0] :
            d = - dtmp
        bP = Calculations.Bisector(P[0], P[1], P[2])
        perpF = Calculations.PerpFormula(bP, P[1], P[1])
        C[2] = Calculations.DistantPoint(P[1], P[1], (-1, perpF.subs(x, -1)), d / 3) #Point on (P[1], P[2]) line in 1/3 distance of [P[0],P[1]] (C2 controlpoint)
        if P[0][0] > P[1][0] :
            d = - dtmp
        else :
            d = dtmp
        middlePnt = Calculations.DistantPoint(P[0], P[0], P[1], d / 2) #Middle point of [P[0],P[1]]
        PF1 = Calculations.PerpFormula(P[0], P[1], middlePnt)   #Middle perpendicular of [P[0],P[1]]
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
        PF2 = Calculations.PerpFormula(middlePnt, PF1tmp, C[2]) #Perpendicular of PF1 that goes trough C[2]
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
        PF1xPF2 = Calculations.LineIntersect(middlePnt, PF1tmp, C[2], PF2tmp)   #Intersection point of PF1 and PF2
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
        Bx, By = CalculateBezier.BezierFormulaComp(C)
        Bezier.append([Bx, By])
        lenP = len(P)
        for i in range(1, lenP + 1): #Calculate controlpoints for P[i], P[i + 1] segment
            if lenP > 3 and i <= lenP - 3 :
                dtmp = Calculations.PointDist(P[i], P[i + 1]) / 3
                d = dtmp
                if P[i][0] < P[i + 2][0] :
                    d = - dtmp
                bP = Calculations.Bisector(P[i], P[i + 1], P[i + 2])
                perpF = Calculations.PerpFormula(bP, P[i + 1], P[i + 1])
                C[2] = Calculations.DistantPoint(P[i + 1], P[i + 1], (-1, perpF.subs(x, -1)), d)
                C[0] = P[i]
                C[3] = P[i + 1]
                if P[i - 1] < P[i + 1] :
                    d = dtmp
                else :
                    d = - dtmp
                bP = Calculations.Bisector(P[i - 1], P[i], P[i + 1])
                perpF = Calculations.PerpFormula(bP, P[i], P[i])
                C[1] = Calculations.DistantPoint(P[i], P[i], (-1, perpF.subs(x, -1)), d)
                CP.append([C[0], C[1], C[2], C[3]])
                Bx, By = CalculateBezier.BezierFormulaComp(C)
                Bezier.append([Bx, By])
            elif i == len(P) - 2:
                C[0] = P[i]

                dtmp = Calculations.PointDist(P[i + 1], P[i]) / 3
                d = dtmp
                if P[i][0] > P[i + 1][0] :
                    d = - dtmp
                bP = Calculations.Bisector(P[i - 1], P[i], P[i + 1])
                perpF = Calculations.PerpFormula(bP, P[i], P[i])
                C[1] = Calculations.DistantPoint(P[i], P[i], (-1, perpF.subs(x, -1)), d)
                Ptmp1 = Calculations.DistantPoint(P[i], P[i], P[i + 1], d)
                perpF = Calculations.PerpFormula(P[i + 1], P[i], Ptmp1)
                dtmp = Calculations.PointDist(P[i + 1], P[i]) / 2
                if P[i][0] > P[i + 1][0] :
                    d = - dtmp
                else :
                    d = dtmp
                Ptmp = Calculations.DistantPoint(P[i], P[i], P[i + 1], d)
                perptmpF = Calculations.PerpFormula(C[1], Ptmp1, C[1])
                Ptmp2 = Calculations.LineIntersect(C[1], (-1, perptmpF.subs(x, -1)), Ptmp, (-1, perpF.subs(x, -1)))
                C[2] = [2 * Ptmp2[0] - C[1][0], 2 * Ptmp2[1] - C[1][1]]
                C[3] = P[i + 1]

                CP.append([C[0], C[1], C[2], C[3]])
                Bx, By = CalculateBezier.BezierFormulaComp(C)
                Bezier.append([Bx, By])
                break
        return Bezier, CP

    def BezierFormulaComp(C):
        t = Symbol('t')
        Bx = 0
        By = 0
        for i in range(4):
            BinC = Calculations.BinCoef(3, i)
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

    def PlotBezier(C, clr, transp, lstyle) :
        for i in range(len(C)) :
            t1 = np.linspace(0, 1, 50)
            Bx1 = 0
            By1 = 0
            for j in range(4) :
                BinC = Calculations.BinCoef(3, j)
                Bxtmp = BinC * (1 - t1)**(j) * t1**(4 - j - 1) * (C[i][3 - j][0])    #Create Bx(t) formula
                Bx1 = (Bx1) + Bxtmp
                Bytmp = BinC * (1 - t1)**(j) * t1**(4 - j - 1) * (C[i][3 - j][1])  #Create Bx(t) formula
                By1 = (By1) + Bytmp
            plt.plot(Bx1, By1, color = clr, alpha = transp, linestyle = lstyle)
        return


    #Get parametric line Formula
    def ParamLineFormula(P1, P2):
        t = Symbol('t')
        PLFx = ((1 - t) * P1[0] + t * P2[0])
        PLFy = ((1 - t) * P1[1] + t * P2[1])
        return PLFx, PLFy

    #Get line formula
    def LineFormula(P1, P2):
        x = Symbol('x')
        if P1[0] == P2[0]:
            LF = (P1[0] + 0 * x)
            return LF
        LF = (UnevaluatedExpr((P1[1] - P2[1]) / (P1[0] - P2[0])) * x + (P1[0] * P2[1] - P2[0] * P1[1]) / (P1[0] - P2[0]))
        return LF

    #Getting Bezier (currently from whole polyline)
    def GetBezier(n, Pn, Points):
        P = CalculateBezier.PMtrx(Pn, Points)
        Pn = len(P)
        T = CalculateBezier.TMtrx(n, Pn, P)
        M = CalculateBezier.MMtrx(n)
        C = CalculateBezier.CMtrx(Pn, M, P, T)
        return CalculateBezier.BezierFormula(n, C)

    def BezierFormula(n, C):
        Bx = 0
        By = 0
        t = Symbol('t')
        for i in range(1, n + 1):
            BinC = Calculations.BinCoef(n - 1, i - 1)
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
            BinC = Calculations.BinCoef(n - 1 , i)
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
            dist = d[i - 1] + Calculations.PointDist(P[i-1], P[i])
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

class Rotate(Calculations) :
    def __init__(rotate) :
        pass

    def FindFrame(CP) :
        P11tmp = []
        P12tmp = []
        P21tmp = []
        P22tmp = []
        P31tmp = []
        P32tmp = []
        max1 = 0
        max2 = 0
        P11 = (0, 0)
        P12 = (0, 0)
        P21 = (0, 0)
        P22 = (0, 0)
        lenCP = len(CP)
        for i in range(0, lenCP) :
            for j in range(i, lenCP) :
                maxDist1 = 0
                maxDist2 = 0
                maxDist3 = 0
                dist1 = Calculations.PointDist(CP[i][0][0], CP[j][0][0])
                dist2 = Calculations.PointDist(CP[i][0][0], CP[j][-1][3])
                dist3 = Calculations.PointDist(CP[i][-1][3], CP[j][-1][3])

                a1 = Calculations.LineSlope(P11, P12)
                a2 = Calculations.LineSlope(P21, P22)

                #find the largest distance
                if dist1 >= dist2 and dist1 >= dist3 and dist1 > maxDist1 :
                    maxDist1 = dist1
                    slope = Calculations.LineSlope(CP[i][0][0], CP[j][0][0])
                    alpha1 = np.degrees(np.arctan((slope - a1) / (1 + a1 * slope)))
                    alpha2 = np.degrees(np.arctan((slope - a2) / (1 + a2 * slope)))

                    angle = np.degrees(np.arctan((slope) / (1)))
                if dist2 >= dist1 and dist2 >= dist3 and dist2 > maxDist2:
                    maxDist2 = dist2
                    slope = Calculations.LineSlope(CP[i][0][0], CP[j][-1][3])
                    alpha1 = np.degrees(np.arctan((slope - a1) / (1 + a1 * slope)))
                    alpha2 = np.degrees(np.arctan((slope - a2) / (1 + a2 * slope)))

                    angle = np.degrees(np.arctan((slope) / (1)))
                if dist3 >= dist1 and dist3 >= dist2 and dist3 > maxDist3:
                    maxDist = dist3
                    slope = Calculations.LineSlope(CP[i][-1][3], CP[j][-1][3])
                    alpha1 = np.degrees(np.arctan((slope - a1) / (1 + a1 * slope)))
                    alpha2 = np.degrees(np.arctan((slope - a2) / (1 + a2 * slope)))

                    angle = np.degrees(np.arctan((slope) / (1)))

                angle1 = np.degrees(np.arctan((a1) / (1)))
                angle2 = np.degrees(np.arctan((a2) / (1)))

                #update points placed in second furthest distance
                if (max1 <= max2) or abs(alpha1) < 5 :
                    if maxDist1 > max1 :
                        if abs(angle - angle2) < 5 :
                            continue
                        max1 = maxDist1
                        P11 = CP[i][0][0]
                        P12 = CP[j][0][0]
                    if maxDist2 > max1 :
                        if abs(angle - angle2) < 5 :
                            continue
                        max1 = maxDist2
                        P11 = CP[i][0][0]
                        P12 = CP[j][-1][3]
                    if maxDist3 > max1 :
                        if abs(angle - angle2) < 5 :
                            continue
                        max1 = maxDist3
                        P11 = CP[i][-1][3]
                        P12 = CP[j][-1][3]
                elif (max2 < max1) or abs(alpha2) < 5 :
                    if maxDist1 > max2 :
                        if abs(angle - angle1) < 5 :
                            continue
                        max2 = maxDist1
                        P21 = CP[i][0][0]
                        P22 = CP[j][0][0]
                    if maxDist2 > max2 :
                        if abs(angle - angle1) < 5 :
                            continue
                        max2 = maxDist2
                        P21 = CP[i][0][0]
                        P22 = CP[j][-1][3]
                    if maxDist3 > max2 :
                        if abs(angle - angle1) < 5 :
                            continue
                        max2 = maxDist3
                        P21 = CP[i][-1][3]
                        P22 = CP[j][-1][3]

        return P11, P12, max1, P21, P22, max2


    def Transformation(points2, line2, CP2) :
#        P111, P112, dist11, P121, P122, dist12 = Rotate.FindFrame(CP1)
#        P211, P212, dist21, P221, P222, dist22 = Rotate.FindFrame(CP2)
#
#        #find left vertical edge of frame in file 1
#        if P111[0] < P112[0] and P111[1] > P112[1] :
#            P11 = P111
#            if P121[1] < P122[1] :
#                P12 = P121
#            else :
#                P12 = P122
#        elif P112[0] < P111[0] and P112[1] > P111[1] :
#            P11 = P112
#            if P121[1] < P122[1] :
#                P12 = P121
#            else :
#                P12 = P122
#        elif P121[0] < P122[0] and P121[1] > P122[1] :
#            P11 = P121
#            if P111[1] < P112[1] :
#                P12 = P111
#            else :
#                P12 = P112
#        elif P122[0] < P121[0] and P122[1] > P121[1] :
#            P11 = P122
#            if P121[1] < P122[1] :
#                P12 = P111
#            else :
#                P12 = P112
#
#        #find respective left vertical edge of frame in file 2
#        if P211[0] < P212[0] and P211[1] < P212[1] :
#            P21 = P211
#            if P221[1] > P222[1] :
#                P22 = P221
#            else :
#                P22 = P222
#        elif P212[0] < P211[0] and P212[1] < P211[1] :
#            P21 = P212
#            if P221[1] > P222[1] :
#                P22 = P221
#            else :
#                P22 = P222
#        elif P221[0] < P222[0] and P221[1] < P222[1] :
#            P21 = P221
#            if P211[1] > P212[1] :
#                P22 = P211
#            else :
#                P22 = P212
#        elif P222[0] < P221[0] and P222[1] < P221[1] :
#            P21 = P222
#            if P211[1] > P212[1] :
#                P22 = P211
#            else :
#                P22 = P212
#
#        #a1 = Calculations.LineSlope(P11, P12)
#        #a2 = Calculations.LineSlope(P21, P22)
#        #alpha = np.degrees(np.arctan((a2 - a1) / (1 + a1 * a2)))
#        #dist1tmp = Calculations.PointDist(P11, P12)
#        #dist2tmp = Calculations.PointDist(P21, P22)
#        #unit = dist1tmp / dist2tmp
        unit = 10
        points = []
        lenPoints2 = len(points2)
        for i in range(lenPoints2):
            #print(points2[i][0], " - ", P11[0], " * ", unit)
            #move
            points.append([])
            points[i].append((points2[i][0]) * unit)
            if points2[i][1] == 0.0 :
                points[i].append(( points2[i][1]) * unit)
            else :
                points[i].append((- points2[i][1]) * unit)
            #rotate
#            points[i][0] = np.cos(alpha) * (points[i][0] - P11[0]) - np.sin(alpha) * (points[i][1] - P11[1]) + P11[0]
#            points[i][1] = np.sin(alpha) * (points[i][0] - P11[0]) + np.cos(alpha) * (points[i][1] - P11[1]) + P11[1]

            points[i] = tuple(points[i])
        points2.clear()

        CP = []
        lenCP2 = len(CP2)
        for i in range(lenCP2) :
            CP.append([])
            lenCP2i = len(CP2[i])
            for j in range(lenCP2i) :
                CP[i].append([])
                for k in range(4) :
                    #move
                    CP[i][j].append([])
                    CP[i][j][k].append((CP2[i][j][k][0]) * unit)
                    if CP2[i][j][k][1] == 0.0 :
                        CP[i][j][k].append((CP2[i][j][k][1]) * unit)
                    else :
                        CP[i][j][k].append((- CP2[i][j][k][1]) * unit)

                    #rotate
#                    CP[i][j][k][0] = np.cos(alpha) * (CP[i][j][k][0] - P11[0]) - np.sin(alpha) * (CP[i][j][k][1] - P11[1]) + P11[0]
#                    CP[i][j][k][1] = np.sin(alpha) * (CP[i][j][k][0] - P11[0]) + np.cos(alpha) * (CP[i][j][k][1] - P11[1]) + P11[1]

                    CP[i][j][k] = tuple(CP[i][j][k])
        CP2.clear()
        return points, line2, CP


class File1(CalculateBezier, Calculations) :
    def __init__(file1, filename) :
        file1.file = dxfgrabber.readfile(filename)

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
            return File1.SortInsertPos(P, points, l, med)
        if points[med][0] < P[0] :
            return File1.SortInsertPos(P, points, med + 1, r)
        if points[med][0] == P[0] :
            if points[med][1] == P[1] :
                return -1 #there already exists identical point, therefore we will not save it
            if points[med][1] > P[1] :
                return File1.SortInsertPos(P, points, l, med)
            if points[med][1] < P[1] :
                return File1.SortInsertPos(P, points, med + 1, r)

    #Binary sort (for abscissa)
    def SortInsertPosLines(P, CP, l, r): #l - left side, r - right side of segment in array
        if r == -1:
            return 0
        if l == r :
            if CP[l][0][0][0] == P[0] and CP[l][0][0][0] == CP[l][0][1][0] and CP[l][0][0][1] == CP[l][0][1][1] :
                if CP[l][0][0][1] == P[1] :
                    return -1
                if CP[l][0][0][1] < P[1] :
                    return l + 1
                return l
            if CP[l][0][0][0] < P[0] :
                return l + 1
            return l
        med = l + (r - l) // 2
        if CP[med][0][0][0] > P[0] :
            return File1.SortInsertPosLines(P, CP, l, med)
        if CP[med][0][0][0] < P[0] :
            return File1.SortInsertPosLines(P, CP, med + 1, r)
        if CP[med][0][0][0] == P[0] :
            if CP[med][0][0][1] == P[1] :
                return -1 #there already exists identical point, therefore we will not save it
            if CP[med][0][0][1] > P[1] :
                return File1.SortInsertPosLines(P, CP, l, med)
            if CP[med][0][0][1] < P[1] :
                return File1.SortInsertPosLines(P, CP, med + 1, r)

    def Lines(P1, P2, Bezier1, CP1, linesSum) :
        if P1[0] < P2[0] :
            lineStart1 = [P1[0], P1[1]]
            lineEnd1 = [P2[0], P2[1]]
        elif P1[0] > P2[0] :
            lineStart1 = [P2[0], P2[1]]
            lineEnd1 = [P1[0], P1[1]]
        elif P1[1] < P2[1] :
            lineStart1 = [P1[0], P1[1]]
            lineEnd1 = [P2[0], P2[1]]
        else :
            lineStart1 = [P2[0], P2[1]]
            lineEnd1 = [P1[0], P1[1]]

        linetmp = []
#                try:
#                    CP1.index([[lineStart1, lineStart1, lineEnd1, lineEnd1]])
#                except:
        pos = File1.SortInsertPosLines(lineStart1, CP1, 0, linesSum)
        if pos == -1 :
            pos = linesSum

        Bezier1.insert(pos, [CalculateBezier.BezierFormulaComp([lineStart1, lineStart1, lineEnd1, lineEnd1])])
        CP1.insert(pos, [[lineStart1, lineStart1, lineEnd1, lineEnd1]])

        linesSum += 1

        x = [lineStart1[0], lineEnd1[0]]
        y = [lineStart1[1], lineEnd1[1]]

        plt.plot(x, y, color = firstFileCol, alpha = 0.2)
        plt.plot(lineStart1[0], lineStart1[1], 'o', ms = 1, color = firstFileCol, alpha = 0.5)
        plt.plot(lineEnd1[0], lineEnd1[1], 'o', ms = 1, color = firstFileCol, alpha = 0.5)

        return Bezier1, CP1, linesSum

    def CalculateObjects(file1) :
        type1 = [entity.dxftype for entity in file1.file.entities]
        output1 = [entity for entity in file1.file.entities]
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
        linesSum = -1
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
                r = len(points1) - 1
                pos = File1.SortInsertPos(point1, points1, 0, r)
                if pos != -1 :
                    points1.insert(pos, [point1[0], point1[1]])
                    x = [point1[0]]
                    y = [point1[1]]
                    plt.plot(x, y, 'o', ms = 1,  color = firstFileCol, alpha = 0.65)
            #Line
            if entity.dxftype == 'LINE':
                P1 = entity.start
                P2 = entity.end
                Bezier1, CP1, linesSum = File1.Lines(P1, P2, Bezier1, CP1, linesSum)
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
                if len(PolylinePoints1) == 2 :
                    Bezier1, CP1, linesSum = File1.Lines(PolylinePoints1[0], PolylinePoints1[1], Bezier1, CP1, linesSum)
                    continue
                if len(PolylinePoints1) == 3 :
                    PolylinePoints1.append(PolylinePoints[-1])
                Beziertmp, CP = CalculateBezier.CompositeBezier(PolylinePoints1, 1)
                Bezier1.append(Beziertmp)
                CP1.append(CP)
            #LWPolyline
            if entity.dxftype == 'LWPOLYLINE':
                LWPolylinePoints1 = entity.points
                Beziertmp, CP = CalculateBezier.CompositeBezier(LWPolylinePoints1, 1)
                Bezier1.append(Beziertmp)
                CP1.append(CP)
            #Spline
            if entity.dxftype == 'SPLINE':
                splineControlPoints1 = entity.control_points
                splineCP1.append(splineControlPoints1)

        linetmp = []
        linetmp.append([])
        #compress lines, so there are not seperate line segments
        lenCP1 = len(CP1)
        j = 0
        while j < len(CP1) :
            if (CP1[j][0][0][0] == CP1[j][0][1][0] and CP1[j][0][0][1] == CP1[j][0][1][1] and CP1[j][0][2][0] == CP1[j][0][3][0] and CP1[j][0][2][1] == CP1[j][0][3][1]) :
                #length of CP1 is changing as lines are being compressed
                i = j + 1
                while i < len(CP1) :
                    if len(CP1[i]) == 1 and (CP1[i][0][0][0] == CP1[i][0][1][0] and CP1[i][0][0][1] == CP1[i][0][1][1] and CP1[i][0][2][0] == CP1[i][0][3][0] and CP1[i][0][2][1] == CP1[i][0][3][1]) :
                        #this is yet not in use, if statement - kreiss to not work
                        a1 = Calculations.LineSlope(CP1[j][-1][0], CP1[j][-1][3])
                        a3 = Calculations.LineSlope(CP1[i][0][0], CP1[i][0][3])
                        if a1 == a3 :
                            t = Symbol('t')
                            a1Fx, a1Fy = CalculateBezier.ParamLineFormula(CP1[j][-1][0], CP1[j][-1][3])
                            a3Fx, a3Fy = CalculateBezier.ParamLineFormula(CP1[i][0][0], CP1[i][0][3])
                            y1 = a1Fx.subs(t, 0)
                            y2 = a1Fx.subs(t, 0)
                            if y1 == y2:
                                if CP1[i][0][0][0] >= CP1[j][-1][0][0] and CP1[i][0][3][0] >= CP1[j][-1][3][0] and CP1[i][0][0][0] <= CP1[j][-1][3][0] and CP1[j][0][0][1] == CP1[i][0][3][1]:
                                    if CP1[i][0][0][0] == CP1[j][-1][0][0] and CP1[i][0][3][0] == CP1[j][-1][3][0] and CP1[i][0][0][1] == CP1[j][-1][0][1] and CP1[i][0][3][1] == CP1[j][-1][3][1]:
                                        i += 1
                                        continue
                                    CP1[j][-1][3] = CP1[i][0][3]
                                    CP1[j][-1][2] = CP1[i][0][3]
                                    Bezier1[j][-1] = CalculateBezier.BezierFormulaComp(CP1[j][-1])
                                    Bezier1.pop(i)
                                    CP1.pop(i)
                                    i -= 1
                                elif CP1[i][0][0][1] >= CP1[j][-1][0][1] and CP1[i][0][3][1] >= CP1[j][-1][3][1] and CP1[i][0][0][1] <= CP1[j][-1][3][1] and CP1[j][0][0][0] == CP1[i][0][3][0]:
                                    if CP1[i][0][0][0] == CP1[j][-1][0][0] and CP1[i][0][3][0] == CP1[j][-1][3][0] and CP1[i][0][0][1] == CP1[j][-1][0][1] and CP1[i][0][3][1] == CP1[j][-1][3][1]:
                                        i += 1
                                        continue
                                    CP1[j][-1][3] = CP1[i][0][3]
                                    CP1[j][-1][2] = CP1[i][0][3]
                                    Bezier1[j][-1] = CalculateBezier.BezierFormulaComp(CP1[j][-1])
                                    Bezier1.pop(i)
                                    CP1.pop(i)
                                    i -= 1
                    i += 1

            j += 1
        return points1, line1, Bezier1, CP1


class File2(CalculateBezier, Rotate) :
    def __init__(file2, filename) :
        file2.file = open(filename)
    #Binary sort (for abscissa)
    def SortInsertPosLines(P, CP, l, r): #l - left side, r - right side of segment in array
        if r == -1 :
            return 0
        if l == r :
            if CP[l][0][0][0] == P[0] and CP[l][0][0][0] == CP[l][0][1][0] and CP[l][0][0][1] == CP[l][0][1][1] :
                if CP[l][0][0][1] == P[1] :
                    return -1
                if CP[l][0][0][1] < P[1] :
                    return l + 1
                return l
            if CP[l][0][0][0] < P[0] :
                return l + 1
            return l
        med = l + (r - l) // 2
        if CP[med][0][0][0] > P[0] :
            return File2.SortInsertPosLines(P, CP, l, med)
        if CP[med][0][0][0] < P[0] :
            return File2.SortInsertPosLines(P, CP, med + 1, r)
        if CP[med][0][0][0] == P[0] :
            if CP[med][0][0][1] == P[1] :
                return -1 #there already exists identical point, therefore we will not save it
            if CP[med][0][0][1] > P[1] :
                return File2.SortInsertPosLines(P, CP, l, med)
            if CP[med][0][0][1] < P[1] :
                return File2.SortInsertPosLines(P, CP, med + 1, r)

    def GetValue(text, pos, length) :
        value = ''
        for i in range(pos, length + 1) :
            if text[i] == '\"' :
                i += 1
                while text[i] != '\"' :
                    value += text[i]
                    i += 1
                value = float(value)
                return value, i
    def GetStrValue(text, pos, length) :
        value = ''
        for i in range(pos, length + 1) :
            if text[i] == '\"' :
                i += 1
                while text[i] != '\"' :
                    value += text[i]
                    i += 1
                return value, i
#returns control points or polylinePoints based on which way the line is defined in given file
    def GetControlPoints(text, pos, length) :
        CP = []
        polyPoints = []
        value = ''
        while pos < length :
            if text[pos] == '\"' :
                return CP, polyPoints, pos
            #get first point in path
            if text[pos] == 'M' :
                pos += 1
                while text[pos] != ' ' and text[pos] != '\"' :
                    value += text[pos]
                    pos += 1
                pos += 1
                value = float(value)
                P1x = value
                value = ''

                while text[pos] != ' ' and text[pos] != '\"' :
                    value += text[pos]
                    pos += 1
                value = float(value)
                P1y = value
                value = ''
            #get controlpoints of given Bezier
            if text[pos] == 'C' :
                CP.append([])
                CP[0].append(P1x)
                CP[0].append(P1y)
                pos += 1
                while text[pos] != ' ' and text[pos] != '\"' :
                    value += text[pos]
                    pos += 1
                pos += 1
                value = float(value)
                CP.append([])
                CP[1].append(value)
                value = ''

                while text[pos] != ' ' and text[pos] != '\"' :
                    value += text[pos]
                    pos += 1
                pos += 1
                value = float(value)
                CP[1].append(value)
                value = ''

                while text[pos] != ' ' and text[pos] != '\"' :
                    value += text[pos]
                    pos += 1
                pos += 1
                value = float(value)
                CP.append([])
                CP[2].append(value)
                value = ''

                while text[pos] != ' ' and text[pos] != '\"' :
                    value += text[pos]
                    pos += 1
                pos += 1
                value = float(value)
                CP[2].append(value)
                value = ''

                while text[pos] != ' ' and text[pos] != '\"' :
                    value += text[pos]
                    pos += 1
                pos += 1
                value = float(value)
                CP.append([])
                CP[3].append(value)
                value = ''

                while text[pos] != ' ' and text[pos] != '\"' :
                    value += text[pos]
                    pos += 1
                value = float(value)
                CP[3].append(value)
                value = ''
                pos -= 1

            #get parameters of an ellipse (first parameter is under 'M')
            if text[pos] == 'A' :
                pos += 1
                #get horizontal axis
                while text[pos] != ' ' and text[pos] != '\"' :
                    value += text[pos]
                    pos += 1
                pos += 1
                value = float(value)
                rx = value
                value = ''
                #get vertical axis
                while text[pos] != ' ' and text[pos] != '\"' :
                    value += text[pos]
                    pos += 1
                pos += 1
                value = float(value)
                ry = value
                value = ''
                #get ellipse rotation angle
                while text[pos] != ' ' and text[pos] != '\"' :
                    value += text[pos]
                    pos += 1
                pos += 1
                value = float(value)
                angle = value
                value = ''
                #get the properties if the large or small arc is used
                while text[pos] != ' ' and text[pos] != '\"' :
                    value += text[pos]
                    pos += 1
                pos += 1
                value = float(value)
                largeArcFlag = value
                value = ''
                #get which of two possible ellipses is drawn
                while text[pos] != ' ' and text[pos] != '\"' :
                    value += text[pos]
                    pos += 1
                pos += 1
                value = float(value)
                sweepFlag = value
                value = ''
                #get the second point on ellipse between which the arc will be drawn
                while text[pos] != ' ' and text[pos] != '\"' and text[pos] != '.':
                    value += text[pos]
                    pos += 1
                if text[pos] == '.' :
                    pos += 1
                    valuetmp = '.'
                    while text[pos] != ' ' and text[pos] != '\"' and text[pos] != '.':
                        valuetmp += text[pos]
                        pos += 1
                value = float(value)
                valuetmp = float(valuetmp)
                if value < 0 :
                    P2x = value - valuetmp
                else :
                    P2x = value + valuetmp
                pos += 1
                value = ''
                while text[pos] != ' ' and text[pos] != '\"' and text[pos] != '.' :
                    value += text[pos]
                    pos += 1
                if text[pos] == '.' :
                    pos += 1
                    valuetmp = '.'
                    while text[pos] != ' ' and text[pos] != '\"' and text[pos] != '.':
                        valuetmp += text[pos]
                        pos += 1
                pos -= 1
                value = float(value)
                valuetmp = float(valuetmp)
                if value < 0 :
                    P2y = value - valuetmp
                else :
                    P2y = value + valuetmp
                value = ''
                #for now its given that every ellipse will be a circle

                #let 'draw' radical axis, then calculate two possible central Points
                #based on radius, and half of radical axis, and we will use Pithagorean theorem
                midx = (P1x + P2x) / 2
                midy = (P1y + P2y) / 2
                dist = Calculations.PointDist((P1x, P1y), (midx, midy))

                edge3 = np.sqrt(round((rx)**2 - (dist)**2, 5))
                #get perpendicular formula which also will coincide with 3rd edge
                PF = Calculations.PerpFormula((P1x, P1y), (midx, midy), (midx, midy))

                x = Symbol('x')
                tmpx = -1
                tmpy = PF.subs(x, tmpx)

                #get two possible circle centers
                c1 = []
                c2 = []
                c1 = Calculations.DistantPoint((midx, midy), (midx, midy), (tmpx, tmpy), edge3)
                c2 = Calculations.DistantPoint((midx, midy), (midx, midy), (tmpx, tmpy), - edge3)
                c = []
                #get the correct arc based on largeArcFlag and sweepFlag from 4 possible options
                radToDeg = 57.2957795

                a1 = Calculations.LineSlope(c1, (P1x, P1y))
                a2 = Calculations.LineSlope(c1, (P2x, P2y))
                t1 = np.arctan(a1 / 1)
                if P1x < c1[0] :
                    t1 = ((t1 + 180 / radToDeg) % (2 * np.pi))
                a2 = Calculations.LineSlope(c1, (P2x, P2y))
                t2 = np.arctan(a2 / 1)
                if P2x < c1[0] :
                    t2 = ((t2 + 180 / radToDeg) % (2 * np.pi))
                if t1 > t2 :
                    arcang = 360 - t1 + t2
                else :
                    arcang = t2 - t1

                if largeArcFlag == 0 and sweepFlag == 0 and 360 - arcang <= 180 :
                    c = c1
                elif largeArcFlag == 0 and sweepFlag == 0 :
                    c = c2
                elif largeArcFlag == 0 and sweepFlag == 1 and arcang <= 180 :
                    c = c1
                elif largeArcFlag == 0 and sweepFlag == 1 :
                    c = c2
                elif largeArcFlag == 1 and sweepFlag == 0 and 360 - arcang > 180 :
                    c = c1
                elif largeArcFlag == 1 and sweepFlag == 0 :
                    c = c2
                elif largeArcFlag == 1 and sweepFlag == 1 and arcang > 180 :
                    c = c1
                else :
                    c = c2

                a1 = Calculations.LineSlope(c, (P1x, P1y))
                a2 = Calculations.LineSlope(c, (P2x, P2y))
                t1 = np.arctan(a1 / 1)
                if P1x < c[0] :
                    t1 = ((t1 + 180 / radToDeg) % (2 * np.pi))
                a2 = Calculations.LineSlope(c, (P2x, P2y))
                t2 = np.arctan(a2 / 1)
                if P2x < c[0] :
                    t2 = ((t2 + 180 / radToDeg) % (2 * np.pi))
                t = t1

                if sweepFlag == 0 :
                    cnt = -2 * np.pi / 18
                else :
                    cnt = 2 * np.pi / 18

                if (t1 < t2 and sweepFlag == 1) or (t1 > t2 and sweepFlag == 0):
                    while t < t2 :
                        xprim = rx * np.cos(t) + c[0]
                        yprim = rx * np.sin(t) + c[1]
                        polyPoints.append((xprim, yprim))
                        t += cnt
                    if t - cnt == t1 :
                        t += cnt / 2
                        xprim = rx * np.cos(t) + c[0]
                        yprim = rx * np.sin(t) + c[1]
                        polyPoints.append((xprim, yprim))
                    xprim = rx * np.cos(t2) + c[0]
                    yprim = rx * np.sin(t2) + c[1]
                    polyPoints.append((xprim, yprim))
                else :
                    t = t2
                    while t < 360 / radToDeg :
                        xprim = rx * np.cos(t) + c[0]
                        yprim = rx * np.sin(t) + c[1]
                        polyPoints.append((xprim, yprim))
                        t += cnt
                    xprim = rx * np.cos(0) + c[0]
                    yprim = rx * np.sin(0) + c[1]
                    polyPoints.append((xprim, yprim))
                    t = 0
                    while t < t1 :
                        xprim = rx * np.cos(t) + c[0]
                        yprim = rx * np.sin(t) + c[1]
                        polyPoints.append((xprim, yprim))
                        t += cnt
                    xprim = rx * np.cos(t1) + c[0]
                    yprim = rx * np.sin(t1) + c[1]
                    polyPoints.append((xprim, yprim))
            pos += 1
        return CP, polyPoints, pos


    def WriteCircle(text, pos, length, points, polylinePoints) :
        cx = -1
        cy = +1
        isShape = False
        while pos < length :
            attrib = ''
            if text[pos] == '/' :
                break
            if text[pos] == ' ' :
                pos += 1
                while text[pos] != '=' :
                    if text[pos] == '>' :
                        break
                    attrib += text[pos]
                    pos += 1
                pos -= 1
                if attrib == 'cx' :
                    cx, pos = File2.GetValue(text, pos, length)
                    continue
                if attrib == 'cy' :
                    cy, pos = File2.GetValue(text, pos, length)
                    continue
                if attrib == 'r' :
                    r, pos = File2.GetValue(text, pos, length)
                    continue
                if attrib == 'pathLength' :

                    continue
                if attrib == 'stroke' :
                    color, pos = File2.GetStrValue(text, pos, length)
                    if color == '#000000' :
                        isShape = True
                    else :
                        return points, polylinePoints, isShape
            pos += 1
        if r <= 0.1 :
            if cx == -1 :
                cx = 0
            if cy == 1 :
                cy = 0
            points.append((cx, cy))
        else :
            t = 0
            while t <= 360 :
                x = r * np.cos(t)
                y = r * np.sin(t)
                polylinePoints.append((x, y))
                t += 36

        return points, polylinePoints, isShape

                #these we will add later...or not
                #clip-path, clip-rule, color, color-interpolation, color-rendering, cursor, display, fill, fill-opacity, fill-rule, filter, mask, opacity, pointer-events, shape-rendering, stroke, stroke-dasharray, stroke-dashoffset, stroke-linecap, stroke-linejoin, stroke-miterlimit, stroke-opacity, stroke-width, transform, vector-effect, visibility
                #class, style, also xml prikoli

    #for text ellement its different


    def WriteText(text, pos, length, pointTags, group) :
        tag = ''
        while pos < length :
            attrib = ''
            if text[pos] == '/' :
                break
            if text[pos] == ' ' :
                pos += 1
                while text[pos] != '=' and pos < length :
                    attrib += text[pos]
                    pos += 1
            if text[pos] == '>' :
                pos += 1
                while text[pos] != '<' and pos < length :
                    tag += text[pos]
                    pos += 1
                if attrib == 'x' :
                    x, pos = File2.GetValue(text, pos, length)
                if attrib == 'y' :
                    y, pos = File2.GetValue(text, pos, length)
                if attrib == 'dx' :
                    dx, pos = File2.GetValue(text, pos, length)
                if attrib == 'dy' :
                    dy, pos = File2.GetValue(text, pos, length)
                if attrib == 'rotate' :
                    alpha, pos = File2.GetValue(text, pos, length)
    #            if attrib == 'lengthAdjust' :

    #            if attrib == 'textLength' :

            pos += 1
        pointTags.append(tag)
        return pointTags

    def WriteEllipse(text, pos, length) :
        isShape = False
        while pos < length :
            attrib = ''
            if text[pos] == '/' :
                break
            if text[pos] != ' ' :
                pos = pos + 1
                while text[pos] != '=' and pos < length :
                    if text[pos] == '>' :
                        break
                    attrib += text[pos]
                    pos += 1

                if attrib == 'cx' :

                    continue
                if attrib == 'cy' :

                    continue
                if attrib == 'rx' :

                    continue
                if attrib == 'ry' :

                    continue
                if attrib == 'pathLength' :

                    continue
                if attrib == 'stroke' :
                    isShape = True
                else :
                    return
            pos += 1
        return

    def WriteLine(text, pos, length) :
        isShape = False
        while pos < length :
            attrib = ''
            if text[pos] == '/' :
                break
            if text[pos] == ' ' :
                pos += 1
                while text[pos] != '=' and pos < length :
                    if text[pos] == '>' :
                        break
                    attrib += text[pos]
                    pos += 1
                pos -= 1

                if attrib == 'x1' :
                    x1, pos = File2.GetValue(text, pos + 1, length)
                if attrib == 'y1' :
                    y1, pos = File2.GetValue(text, pos + 1, length)
                if attrib == 'x2' :
                    x2, pos = File2.GetValue(text, pos + 1, length)
                if attrib == 'y2' :
                    y2, pos = File2.GetValue(text, pos + 1, length)
    #            if attrib == 'pathLength' :
                if attrib == 'stroke' :
                    color, pos = File2.GetStrValue(text, pos + 1, length)
                    if color == '#000000' :
                        isShape = True
                    else :
                        return [], isShape
            pos += 1
        if x1 < x2 :
            CP = ([[x1, y1], [x1, y1], [x2, y2], [x2, y2]])
        elif x1 > x2 :
            CP = ([[x2, y2], [x2, y2], [x1, y1], [x1, y1]])
        elif y1 < y2 :
            CP = ([[x1, y1], [x1, y1], [x2, y2], [x2, y2]])
        else :
            CP = ([[x2, y2], [x2, y2], [x1, y1], [x1, y1]])
        return CP, isShape

    def WritePath(text, pos, length, polylinePoints) :
        isShape = False
        while pos < length :
            attrib = ''
            if text[pos] == '/' :
                break
            if text[pos] == ' ' :
                pos += 1
                while text[pos] != '=' and pos < length :
                    if text[pos] == '>' :
                        break

                    attrib += text[pos]
                    pos += 1
                if attrib == 'd' :
                    CPoints, polyPoints, pos = File2.GetControlPoints(text, pos + 2, length)
                    continue
                if attrib == 'pathLength' :

                    continue
                if attrib == 'stroke' :
                    color, pos = File2.GetStrValue(text, pos + 1, length)
                    if color == '#000000' :
                        isShape = True
                    else :
                        return CPoints, polylinePoints, isShape
            pos += 1
        if len(polyPoints) > 0:
            polylinePoints.append(polyPoints)
        return CPoints, polylinePoints, isShape


    def WritePolygon(text, pos, length) :
        while pos < length :
            attrib = ''
            if text[pos] == '/' :
                break
            if text[pos] != ' ' :
                pos += 1
                while text[pos] != '=' and pos < length :
                    if text[pos] == '>' :
                        break
                    attrib += text[pos]
                    pos += 1

                if attrib == 'points' :

                    continue
                if attrib == 'pathLength' :

                    continue
            pos += 1

    def WritePolyline(text, pos, length) :
        while pos < length :
            attrib = ''
            if text[pos] == '/' :
                break
            if text[pos] != ' ' :
                pos += 1
                while text[pos] != '=' and pos < length:
                    if text[pos] == '>' :
                        break
                    attrib += text[pos]
                    pos += 1

                if attrib == 'points' :

                    continue
                if attrib == 'pathLength' :

                    continue
            pos += 1

    def WriteRect(text, pos, length) :
        while pos < length :
            attrib = ''
            if text[pos] == '/' :
                break
            if text[pos] != ' ' :
                pos += 1
                while text[pos] != '=' and pos < length:
                    if text[pos] == '>' :
                        break
                    attrib += text[pos]
                    pos += 1

                if attrib == 'x' :

                    continue
                if attrib == 'y' :

                    continue
                if attrib == 'width' :

                    continue
                if attrib == 'height' :

                    continue
                if attrib == 'rx' :

                    continue
                if attrib == 'ry' :

                    continue
                if attrib == 'pathLength' :

                    continue
            pos += 1

    def CalculateObjects(file2) :
        textLines = file2.file.readlines()

        group = False
        points = []
        polylinePoints = []
        pointTags = []
        line = []
        CP2 = []
        linesSum = -1

        for textLine in textLines :
            lenText = len(textLine)
            for i in range(lenText) :
                element = ''
                if textLine[i] == '<' :
                    j = i + 1
                    while (textLine[j] != '>' and textLine[j] != ' ') :
                        element += textLine[j]
                        i = j
                        j += 1
                        #<svg>

                    if element == 'g':
                        group = True
                    if element == '/text' :
                        group = False
                    if element == 'circle' :
                        points, polylinePoints, isShape = File2.WriteCircle(textLine, i + 1, lenText, points, polylinePoints)
                    if element == 'text' :
                        pointTags = File2.WriteText(textLine, i + 1, lenText, pointTags, group)
                    if element == 'ellipse' :
                        File2.WriteEllipse(textLine, i + 1, lenText)
                    if element == 'line' :
                        CPtmp, isShape = File2.WriteLine(textLine, i + 1, lenText)
                        if isShape == False :
                            continue
#                        try :
#                            CP2.index([CPtmp])
#                        except :
                        pos = File2.SortInsertPosLines(CPtmp[0], CP2, 0, linesSum)
                        CP2.insert(pos, [CPtmp])
                        linesSum += 1
                    if element == 'path' :
                        CPtmp, polylinePoints, isShape = File2.WritePath(textLine, i + 1, lenText, polylinePoints)
                        if isShape == False :
                            continue
                        if len(CPtmp) != 0 :
                            #for sorting we need to make sure Bezier is formed in the right way, meaning that the segments are corresponding
                            if CPtmp[0][0] > CPtmp[3][0] :
                                CPtmp.reverse()
                            elif CPtmp[0][1] > CPtmp[3][1] :
                                CPtmp.reverse()
                            pos = File2.SortInsertPosLines(CPtmp[0], CP2, 0, linesSum)
                            CP2.append([CPtmp])
                            linesSum += 1
                    if element == 'polygon' :
                        File2.WritePolygon(textLine, i + 1, lenText)
                    if element == 'polyline' :
                        File2.WritePolyline(textLine, i + 1, lenText)
                    if element == 'rect' :
                        File2.WriteRect(textLine, i + 1, lenText)
                    break
        #compress lines, so there are not seperate line segments
        lenCP2 = len(CP2)
        j = 0
        while j < len(CP2) :
            #length of CP1 is changing as lines are being compressed
            i = j + 1
            while i < len(CP2) :
                a1 = Calculations.PointDist(CP2[i][0][0], CP2[i][0][1])
                a2 = Calculations.PointDist(CP2[i][0][1], CP2[i][0][2])
                a3 = Calculations.PointDist(CP2[i][0][2], CP2[i][0][3])

                a4 = Calculations.PointDist(CP2[j][-1][0], CP2[j][-1][1])
                a5 = Calculations.PointDist(CP2[j][-1][1], CP2[j][-1][2])
                a6 = Calculations.PointDist(CP2[j][-1][2], CP2[j][-1][3])
                if (((a1 == 0 and a2 == 0) or (a2 == 0 and a3 == 0) or (a1 == 0 and a3 == 0)) and ((a4 == 0 and a5 == 0) or (a5 == 0 and a6 == 0) or (a4 == 0 and a6 == 0))) :
                    a1 = Calculations.LineSlope(CP2[j][-1][0], CP2[j][-1][3])
                    a3 = Calculations.LineSlope(CP2[i][0][0], CP2[i][0][3])
                    if a1 == a3 :
                        t = Symbol('t')
                        a1Fx, a1Fy = CalculateBezier.ParamLineFormula(CP2[j][-1][0], CP2[j][-1][3])
                        a3Fx, a3Fy = CalculateBezier.ParamLineFormula(CP2[i][0][0], CP2[i][0][3])
                        y1 = a1Fx.subs(t, 0)
                        y2 = a1Fx.subs(t, 0)
                        if y1 == y2:
                            if CP2[i][0][0][0] > CP2[j][-1][0][0] and CP2[i][0][3][0] > CP2[j][-1][3][0] and CP2[i][0][0][0] <= CP2[j][-1][3][0] and CP2[j][0][0][1] == CP2[i][0][3][1] :
                                if CP2[i][0][0][0] == CP2[j][-1][0][0] and CP2[i][0][3][0] == CP2[j][-1][3][0] and CP2[i][0][0][1] == CP2[j][-1][0][1] and CP2[i][0][3][1] == CP2[j][-1][3][1] :
                                    i += 1
                                    continue
                                CP2[j][-1][3] = CP2[i][0][3]
                                CP2[j][-1][2] = CP2[i][0][3]
                                CP2.pop(i)
                                i -= 1
                            elif CP2[i][0][0][1] > CP2[j][-1][0][1] and CP2[i][0][3][1] > CP2[j][-1][3][1] and CP2[i][0][0][1] <= CP2[j][-1][3][1] and CP2[j][-1][0][0] == CP2[i][0][3][0] :
                                if CP2[i][0][0][0] == CP2[j][-1][0][0] and CP2[i][0][3][0] == CP2[j][-1][3][0] and CP2[i][0][0][1] == CP2[j][-1][0][1] and CP2[i][0][3][1] == CP2[j][0][3][1] :
                                    i += 1
                                    continue
                                CP2[j][-1][3] = CP2[i][0][3]
                                CP2[j][-1][2] = CP2[i][0][3]
                                CP2.pop(i)
                                i -= 1
                else :
                    t = Symbol('t')
                    dist1 = Calculations.PointDist(CP2[j][-1][0], CP2[j][-1][1])
                    dist2 = Calculations.PointDist(CP2[j][-1][1], CP2[j][-1][2])
                    dist3 = Calculations.PointDist(CP2[j][-1][2], CP2[j][-1][3])
                    if dist3 != 0 :
                        P11 = CP2[j][-1][2]
                        P12 = CP2[j][-1][3]
                        a1 = round(Calculations.LineSlope(CP2[j][-1][2], CP2[j][-1][3]), 3)
                        a1Fx, a1Fy = CalculateBezier.ParamLineFormula(CP2[j][-1][2], CP2[j][-1][3])
                    elif dist2 != 0 :
                        P11 = CP2[j][-1][1]
                        P12 = CP2[j][-1][3]
                        a1 = round(Calculations.LineSlope(CP2[j][-1][1], CP2[j][-1][3]), 3)
                        a1Fx, a1Fy = CalculateBezier.ParamLineFormula(CP2[j][-1][1], CP2[j][-1][3])
                    else :
                        P11 = CP2[j][-1][0]
                        P12 = CP2[j][-1][3]
                        a1 = round(Calculations.LineSlope(CP2[j][-1][0], CP2[j][-1][3]), 3)
                        a1Fx, a1Fy = CalculateBezier.ParamLineFormula(CP2[j][-1][0], CP2[j][-1][3])

                    dist1 = Calculations.PointDist(CP2[j][0][0], CP2[j][0][1])
                    dist2 = Calculations.PointDist(CP2[j][0][1], CP2[j][0][2])
                    dist3 = Calculations.PointDist(CP2[j][0][2], CP2[j][0][3])

                    if dist1 != 0 :
                        P21 = CP2[j][0][0]
                        P22 = CP2[j][0][1]
                        a2 = round(Calculations.LineSlope(CP2[j][0][0], CP2[j][0][1]), 3)
                        a2Fx, a2Fy = CalculateBezier.ParamLineFormula(CP2[j][0][0], CP2[j][0][1])
                    elif dist2 != 0 :
                        P21 = CP2[j][0][0]
                        P22 = CP2[j][0][2]
                        a2 = round(Calculations.LineSlope(CP2[j][0][0], CP2[j][0][2]), 3)
                        a2Fx, a2Fy = CalculateBezier.ParamLineFormula(CP2[j][0][0], CP2[j][0][2])
                    else :
                        P21 = CP2[j][0][0]
                        P22 = CP2[j][0][3]
                        a2 = round(Calculations.LineSlope(CP2[j][0][0], CP2[j][0][3]), 3)
                        a2Fx, a2Fy = CalculateBezier.ParamLineFormula(CP2[j][0][0], CP2[j][0][3])

                    dist1 = Calculations.PointDist(CP2[i][-1][0], CP2[i][-1][1])
                    dist2 = Calculations.PointDist(CP2[i][-1][1], CP2[i][-1][2])
                    dist3 = Calculations.PointDist(CP2[i][-1][2], CP2[i][-1][3])
                    if dist3 != 0 :
                        P31 = CP2[i][-1][2]
                        P32 = CP2[i][-1][3]
                        a3 = round(Calculations.LineSlope(CP2[i][-1][2], CP2[i][-1][3]), 3)
                        a3Fx, a3Fy = CalculateBezier.ParamLineFormula(CP2[i][-1][3], CP2[i][-1][2])
                    elif dist2 != 0 :
                        P31 = CP2[i][-1][1]
                        P32 = CP2[i][-1][3]
                        a3 = round(Calculations.LineSlope(CP2[i][-1][1], CP2[i][-1][3]), 3)
                        a3Fx, a3Fy = CalculateBezier.ParamLineFormula(CP2[i][-1][3], CP2[i][-1][1])
                    else :
                        P31 = CP2[i][-1][0]
                        P32 = CP2[i][-1][3]
                        a3 = round(Calculations.LineSlope(CP2[i][-1][0], CP2[i][-1][3]), 3)
                        a3Fx, a3Fy = CalculateBezier.ParamLineFormula(CP2[i][-1][3], CP2[i][-1][0])

                    dist1 = Calculations.PointDist(CP2[i][0][0], CP2[i][0][1])
                    dist2 = Calculations.PointDist(CP2[i][0][1], CP2[i][0][2])
                    dist3 = Calculations.PointDist(CP2[i][0][2], CP2[i][0][3])
                    if dist1 != 0 :
                        P41 = CP2[i][0][0]
                        P42 = CP2[i][0][1]
                        a4 = round(Calculations.LineSlope(CP2[i][0][0], CP2[i][0][1]), 3)
                        a4Fx, a4Fy = CalculateBezier.ParamLineFormula(CP2[i][0][0], CP2[i][0][1])
                    elif dist2 != 0 :
                        P41 = CP2[i][0][0]
                        P42 = CP2[i][0][2]
                        a4 = round(Calculations.LineSlope(CP2[i][0][0], CP2[i][0][2]), 3)
                        a4Fx, a4Fy = CalculateBezier.ParamLineFormula(CP2[i][0][0], CP2[i][0][2])
                    else :
                        P41 = CP2[i][0][0]
                        P42 = CP2[i][0][3]
                        a4 = round(Calculations.LineSlope(CP2[i][0][0], CP2[i][0][3]), 3)
                        a4Fx, a4Fy = CalculateBezier.ParamLineFormula(CP2[i][0][0], CP2[i][0][3])

                    if a1 == a3 or a1 == a4 or a2 == a3 or a2 == a4 :
                        y1 = a1Fx.subs(t, 0)
                        y2 = a2Fx.subs(t, 0)
                        y3 = a3Fx.subs(t, 0)
                        y4 = a4Fx.subs(t, 0)
                        if y1 == y3 or y1 == y4 or y2 == y3 or y2 == y4 :
                            if round(P41[0], 1) == round(P21[0], 1) and round(P41[1], 1) == round(P21[1], 1) and ((P42[0] >= P41[0] and P42[1] >= P41[1] and P21[1] >= P22[1] and P21[1] >= P22[1]) or (P42[0] <= P41[0] and P42[1] <= P41[1] and P21[1] <= P22[1] and P21[1] <= P22[1])) :
                                if CP2[j][-1][3][0] < CP2[i][-1][3][0] :
                                    CP2[j].reverse()
                                    for k in CP2[j]:
                                        k.reverse()
                                    CP2[j].extend(CP2[i])
                                    CP2.pop(i)
                                    i -= 1
                                elif CP2[i][-1][3][0] < CP2[j][-1][3][0] :
                                    CP2[i].reverse()
                                    for k in CP2[i]:
                                        k.reverse()
                                    CP2[j].insert(0, CP2[i][0])
                                    CP2.pop(i)
                                    i -=1
                                elif CP2[i][-1][3][1] < CP2[j][-1][3][1] :
                                    CP2[i].reverse()
                                    for k in CP2[i]:
                                        k.reverse()
                                    CP2[j].insert(0, CP2[i][0])
                                    CP2.pop(i)
                                    i -= 1
                                else :
                                    CP2[j].reverse()
                                    for k in CP2[j]:
                                        k.reverse()
                                    CP2[j].extend(CP2[i])
                                    CP2.pop(i)
                                    i -= 1
                            elif round(P41[0], 1) == round(P12[0], 1) and round(P41[1], 1) == round(P12[1], 1) and ((P42[0] >= P41[0] and P42[1] >= P41[1] and P12[1] >= P11[1] and P12[1] >= P11[1]) or (P42[0] <= P41[0] and P42[1] <= P41[1] and P12[1] <= P11[1] and P12[1] <= P11[1])) :
                                if CP2[j][0][0][0] < CP2[i][-1][3][0] :
                                    CP2[j].extend(CP2[i])
                                    CP2.pop(i)
                                    i -= 1
                                elif CP2[i][-1][3][0] < CP2[j][0][0][0] :
                                    CP2[j].insert(0, CP2[i][0])
                                    CP2.pop(i)
                                    i -= 1
                                elif CP2[i][-1][3][1] < CP2[j][0][0][1] :
                                    CP2[j].insert(0, CP2[i][0])
                                    CP2.pop(i)
                                    i -= 1
                                else :
                                    CP2[j].extend(CP2[i])
                                    CP2.pop(i)
                                    i -= 1
                            elif round(P32[0], 1) == round(P21[0], 1) and round(P32[1], 1) == round(P21[1], 1) and ((P31[0] >= P32[0] and P31[1] >= P32[1] and P21[1] >= P22[1] and P21[1] >= P22[1]) or (P31[0] <= P32[0] and P31[1] <= P32[1] and P21[1] <= P22[1] and P21[1] <= P22[1])) :
                                if CP2[j][-1][3][0] < CP2[i][0][0][0] :
                                    CP2[j].extend(CP2[i])
                                    CP2.pop(i)
                                    i -= 1
                                elif CP2[i][0][0][0] < CP2[j][-1][3][0] :
                                    CP2[j].insert(0, CP2[i][0])
                                    CP2.pop(i)
                                    i -= 1
                                elif CP2[i][0][0][1] < CP2[j][-1][3][1] :
                                    CP2[j].insert(0, CP2[i][0])
                                    CP2.pop(i)
                                    i -= 1
                                else :
                                    CP2[j].extend(CP2[i])
                                    CP2.pop(i)
                                    i -= 1
                            elif round(P32[0], 1) == round(P12[0], 1) and round(P32[1], 1) == round(P12[1], 1) and ((P31[0] >= P32[0] and P31[1] >= P32[1] and P12[1] >= P11[1] and P12[1] >= P22[1]) or (P31[0] <= P32[0] and P31[1] <= P32[1] and P12[1] <= P11[1] and P12[1] <= P11[1])) :
                                if CP2[j][0][0][0] < CP2[i][0][0][0] :
                                    CP2[j].reverse()
                                    for k in CP2[j]:
                                        k.reverse()
                                    CP2[j].extend(CP2[i])
                                    CP2.pop(i)
                                    i -= 1
                                elif CP2[i][0][0][0] < CP2[j][0][0][0] :
                                    CP2[i].reverse()
                                    for k in CP2[i]:
                                        k.reverse()
                                    CP2[j].insert(0, CP2[i][0])
                                    CP2.pop(i)
                                    i -= 1
                                elif CP2[i][0][0][1] < CP2[j][0][0][1] :
                                    CP2[i].reverse()
                                    for k in CP2[i]:
                                        k.reverse()
                                    CP2[j].insert(0, CP2[i][0])
                                    CP2.pop(i)
                                    i -= 1
                                else :
                                    CP2[j].reverse()
                                    for k in CP2[j]:
                                        k.reverse()
                                    CP2[j].extend(CP2[i])
                                    CP2.pop(i)
                                    i -= 1
                i += 1
            j += 1
        Bezier = []
        for i in polylinePoints :
            if len(polylinePoints) == 2 :
                Beziertmp, CPtmp = File2.Lines(polylinePoints[0], PolylinePoints[1], Bezier, CP2, linesSum)
                continue
            Beziertmp, CPtmp = CalculateBezier.CompositeBezier(i, 2)
            CP2.append(CPtmp)
        points, line, CP2 = Rotate.Transformation(points, line, CP2)

        Beziertmp = [[]]
        for i in CP2 :
            for j in i:
                Bx, By = CalculateBezier.BezierFormulaComp(j)
                Beziertmp[0].append([Bx, By])
            Bezier.extend(Beziertmp)

            Beziertmp.clear()
            Beziertmp.append([])

#        for i in polylinePoints :
#            Beziertmp, CPtmp = CalculateBezier.CompositeBezier(i, 2)
#            CP2.append(CPtmp)
#            Bezier.append(Beziertmp)
        return points, line, Bezier, CP2

class Plot(Calculations) :
    def __init__(plot) :
        pass

    def PlotBezier(C, clr, transp, lstyle) :
        for i in range(len(C)) :
            t1 = np.linspace(0, 1, 50)
            Bx1 = 0
            By1 = 0
            for j in range(4) :
                BinC = Calculations.BinCoef(3, j)
                Bxtmp = BinC * (1 - t1)**(j) * t1**(4 - j - 1) * (C[i][3 - j][0])    #Create Bx(t) formula
                Bx1 = (Bx1) + Bxtmp
                Bytmp = BinC * (1 - t1)**(j) * t1**(4 - j - 1) * (C[i][3 - j][1])  #Create Bx(t) formula
                By1 = (By1) + Bytmp
            plt.plot(Bx1, By1, color = clr, alpha = transp, linestyle = lstyle)
        return






#Lets find which points are represented as the 'same' points in both files
#First we will take 3 points from first file and get the distances between them
#Based on those 3 point we are going to take three point sets from second file
#   and compare if they match based on proportions
#(from those three points we save proportion and angle and in the next steps we
#    are going to use vectors to find the corresponding points)
#If we have found the match then - take one point, pair it with all the other
#   points and find their matching pairs in second file
#Because we do that with vector there is no reason to compare each two points





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

class BezierCalculations(Calculations) :
    def __init__(BezierCalc) :
        pass
    #find closest point in Bezier curves array
    def ClosestPntBezier(P, arr, visited, minDist) :
        min = -1
        minPind = 0
        t = Symbol('t')
        lineSt = 0
        for i in range(len(arr)) :
            #as the Bezier curve is made from multiple cubic Bezier curves, we have to compare only the start point of first and end point of last cubic curve that belongs to specific composite Bezier
            if min == -1 :
                if not (arr[i][0][0][0], arr[i][0][0][1]) in visited:
                    mintmp = Calculations.PointDist(P, (arr[i][0][0][0], arr[i][0][0][1]))
                    if mintmp > minDist :
                        min = mintmp
                        minPind = i
                        lineSt = 0
                elif visited[(arr[i][0][0][0], arr[i][0][0][1])] == 1 :
                    mintmp = Calculations.PointDist(P, (arr[i][0][0][0], arr[i][0][0][1]))
                    if mintmp > minDist :
                        min = mintmp
                        minPind = i
                        lineSt = 0
            if not (arr[i][0][0][0], arr[i][0][0][1]) in visited :
                tmpD = Calculations.PointDist(P, (arr[i][0][0][0], arr[i][0][0][1]))
                if tmpD < min and tmpD > minDist:
                    min = tmpD
                    minPind = i
                    lineSt = 0
            elif visited[(arr[i][0][0][0], arr[i][0][0][1])] == 1 :
                tmpD = Calculations.PointDist(P, (arr[i][0][0][0], arr[i][0][0][1]))
                if tmpD < min and tmpD > minDist:
                    min = tmpD
                    minPind = i
                    lineSt = 0
            if not (arr[i][-1][3][0], arr[i][-1][3][1]) in visited :
                tmpD = Calculations.PointDist(P, (arr[i][-1][3][0], arr[i][-1][3][1]))
                if tmpD < min and tmpD > minDist:
                    min = tmpD
                    minPind = i
                    lineSt = 1
            elif visited[(arr[i][-1][3][0], arr[i][-1][3][1])] == 1 :
                tmpD = Calculations.PointDist(P, (arr[i][-1][3][0], arr[i][-1][3][1]))
                if tmpD < min and tmpD > minDist:
                    min = tmpD
                    minPind = i
                    lineSt = 1
        return (arr[minPind][-lineSt][lineSt * 3][0], arr[minPind][-lineSt][lineSt * 3][1]), min

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

    def DistantPointOnBezier(dist, nr, B, lineSt, C, disttmp, param) :
        t = Symbol('t')
        xtmp = C[-lineSt][lineSt * 3][0]
        ytmp = C[-lineSt][lineSt * 3][1]
        x2tmp = xtmp
        y2tmp = ytmp
        cnt = 1
        m = 0.5
        param = lineSt
        dist1 = disttmp
        lenB = len(B)
        while (nr >= 0) and (nr <= lenB - 1) :
            while (param >= 0 and param <= 1) :
                xtmp = x2tmp
                ytmp = y2tmp
                param = abs(round((lineSt - cnt), 3))
                if (param < 0 or param > 1) :
                    if param < 0 :
                        param = 0
                    else :
                        param = 1
                    break
                x2tmp = B[nr][0].subs(t, param)
                y2tmp = B[nr][1].subs(t, param)
                Pdist = Calculations.PointDist((xtmp, ytmp), (x2tmp, y2tmp))
                disttmp += Pdist
                if param == abs(round((lineSt - 1), 3)) :
                    x2tmp1 = C[nr][lineSt * 3][0]
                    y2tmp1 = C[nr][lineSt * 3][1]
                    Pdist1 = Calculations.PointDist((xtmp, ytmp), (x2tmp1, y2tmp1))
                    disttmp1 = disttmp + Pdist1 - Pdist
                    if (disttmp1 <= dist + 0.05 and disttmp1 >= dist - 0.05) :
                        return x2tmp1, y2tmp1, nr, disttmp1, lineSt
                    if disttmp1 < dist and disttmp1 > disttmp :
                        break

                if (disttmp < dist and param == 1):
                    break
                if (disttmp <= dist + 0.05 and disttmp >= dist - 0.05) :
                    return x2tmp, y2tmp, nr, disttmp, param
                elif m == 0 :
                    break
                elif disttmp > dist + 0.05 :
                    cnt = round((cnt - m), 3)
                    m = round(m * 0.5, 3)
                    disttmp -= Pdist
                    continue
                if (param == 0 or param == 1) :
                    break
                cnt = round((cnt + m), 3)

            if (nr == lenB - 1 and lineSt == 0) :
                break
            if (nr == 0 and lineSt == 1) :
                break
            if lineSt == 0 :
                dist1 += BezierCalculations.CubicBezierLen(C[nr])
                disttmp = dist1
                nr += 1
                param = lineSt
                m = 0.5
                cnt = 1
            else:
                dist1 += BezierCalculations.CubicBezierLen(C[nr])
                disttmp = dist1
                nr -= 1
                param = lineSt
                m = 0.5
                cnt = 1

        x2tmp = B[nr][0].subs(t, param)
        y2tmp = B[nr][1].subs(t, param)
        return x2tmp, y2tmp, nr, disttmp, param
    #returns max bezier distance based on equaly distributed points along both lines
    def MaxBezierTDist(B1, B2, lineSt1, lineSt2) :
        #as the lines may not consist of equal number of segments we have to set t value step for each line independly
        len1 = len(B1)
        len2 = len(B2)
        cnt1 = round(1 / (20 / len1), 5)
        cnts1 = 1
        cnt2 = round(1 / (20 / len2), 5)
        cnts2 = 1
        t1 = lineSt1
        t2 = lineSt2
        nr1 = -lineSt1
        nr2 = -lineSt2
        maxDist = 0
        t = Symbol('t')
        if lineSt1 == 1 :
            cnt1 = -cnt1
            cnts1 = -1
        if lineSt2 == 1 :
            cnt2 = -cnt2
            cnts2 = -1
        for i in range(20) :
            P1 = (B1[nr1][0].subs(t, t1), B1[nr1][1].subs(t, t1))
            P2 = (B2[nr2][0].subs(t, t2), B2[nr2][1].subs(t, t2))
            dist = Calculations.PointDist(P1, P2)
            if dist > maxDist :
                maxDist = dist
            t1 += cnt1
            if t1 > 1 or t1 < 0 :
                nr1 += cnts1
                t1 = t1 - cnts1
            t2 += cnt2
            if t2 > 1 or t2 < 0 :
                nr2 += cnts2
                t2 = t2 - cnts2
        return maxDist

    def MaxBezierDist(B1, B2, CP1, CP2, lineSt1, lineSt2) :
        maxDistParam = BezierCalculations.MaxBezierTDist(B1, B2, lineSt1, lineSt2)
        if maxDistParam <= 1 :
            return maxDistParam, CP1[-lineSt1][lineSt1 * 3], CP2[-lineSt2][lineSt2 * 3]
        t = Symbol('t')
        lenB1 = len(B1)
        lenB2 = len(B2)
        max = -1
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
        nr = 0

        #check distance between endpoints
        x2 = CP2[0][0][0]
        y2 = CP2[0][0][1]
        dist1 = Calculations.PointDist((x2, y2), (CP1[0][0][0], CP1[0][0][1]))
        dist2 = Calculations.PointDist((x2, y2), (CP1[-1][3][0], CP1[-1][3][1]))
        if dist1 < dist2 :
            max = dist1
            P11[0] = CP1[0][0][0]
            P11[1] = CP1[0][0][1]
            P22[0] = x2
            P22[1] = y2
        else :
            max = dist2
            P11[0] = CP1[-1][3][0]
            P11[1] = CP1[-1][3][1]
            P22[0] = x2
            P22[1] = y2

        x2 = CP2[-1][3][0]
        y2 = CP2[-1][3][1]
        dist1 = Calculations.PointDist((x2, y2), (CP1[0][0][0], CP1[0][0][1]))
        dist2 = Calculations.PointDist((x2, y2), (CP1[-1][3][0], CP1[-1][3][1]))
        if dist1 < dist2 and dist1 > max :
            max = dist1
            P11[0] = CP1[0][0][0]
            P11[1] = CP1[0][0][1]
            P22[0] = x2
            P22[1] = y2
        elif dist2 <= dist1 and dist2 > max :
            max = dist2
            P11[0] = CP1[-1][3][0]
            P11[1] = CP1[-1][3][1]
            P22[0] = x2
            P22[1] = y2


        for i in range(lenB2) :
            print(i)
            cnt = 0
            param = 0
            while param <= 1 :
                x2 = B2[i][0].subs(t, param)
                y2 = B2[i][1].subs(t, param)
                min = 10000
                #as the line may be composite, we have to find min distance for every segment form line, and choose the smallest one
                for nr in range(lenB1) :
                    dist, P = BezierCalculations.BezierMinDist(B1[nr][0], B1[nr][1], (x2, y2))
                    if dist < min :
                        P1[0] = P[0]
                        P1[1] = P[1]
                        P2[0] = x2
                        P2[1] = y2
                        min = dist
                    if dist > min:
                        nr -= 1
                        break
                if min > max and min != 10000 :
                    max = min
                    P11[0] = P1[0]
                    P11[1] = P1[1]
                    P22[0] = P2[0]
                    P22[1] = P2[1]
                cnt = round(cnt + 0.1, 1)
                param = cnt
        return max, P11, P22

    #Search for distance between point and Bezier curve
    def BezierMinDist(Bx, By, P) :
        t = Symbol('t')
        param = 0

        x1 = Bx.subs(t, 0)
        y1 = By.subs(t, 0)
        x2 = Bx.subs(t, 1)
        y2 = By.subs(t, 1)
        dist1 = Calculations.PointDist((x1, y1), P)
        dist2 = Calculations.PointDist((x2, y2), P)

        if dist1 <= dist2 :
            min = dist1
            minParam = 0
        elif dist1 > dist2 :
             min = dist2
             minParam = 1
        #check if in segments between previously determined points there is not point with smaller distance to the given point
        param = minParam
        dist1 = 0
        dist2 = 0
        x = Bx.subs(t, param)
        y = By.subs(t, param)
        minP = []
        minP.append(x)
        minP.append(y)
        cnt = 0.1
        while (param >= 0 and param <= 1) :
            if param <= 1 - cnt:
                x1 = Bx.subs(t, param + cnt)
                y1 = By.subs(t, param + cnt)
                dist1 = Calculations.PointDist((x1, y1), P)
            else :
                dist1 = 10000
            if param >= 0 + cnt :
                x2 = Bx.subs(t, param - cnt)
                y2 = By.subs(t, param - cnt)
                dist2 = Calculations.PointDist((x2, y2), P)
            else :
                dist2 = 10000
            if (dist1 < min and dist1 < dist2) :
                param += cnt
                min = dist1
                minP[0] = x1
                minP[1] = y1
                continue
            elif (dist2 < min and dist2 < dist1) :
                param -= cnt
                min = dist2
                minP[0] = x2
                minP[1] = y2
                continue
            if (param <= 0 and param >= 1) :
                return min, minP
            elif cnt > 0.0001:
                cnt /= 2
                continue
            else :
                return min, minP

        #This point technically should never been met
        return min, minP

    def CompositeBezierLen(CP) :
        length = 0
        lenCP = len(CP)
        for i in range(lenCP) :
            length += BezierCalculations.CubicBezierLen(CP[i])
        return length

class ObjectComparison(BezierCalculations, Calculations) :
    def __init__(comp) :
        pass

    #Adds objects to unique object type set
    def OneFileObject(P, line, Bezier, CP, visited, fnr, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2) :
        lenCP = len(CP)
        for i in range(lenCP) :
            if (P == (CP[i][0][0][0], CP[i][0][0][1]) or (P == (CP[i][-1][3][0], CP[i][-1][3][1]))) :
                if fnr == 1 :
    #                PlotBezier(CP[i], '#055583', 1, (0, (3, 5, 1, 5)))
                    if uniqueObjectsBezier1.get(i) == None and ObjectsBezier1.get(i) == None and equalObjectsBezier.get(i) == None :
                        length = BezierCalculations.CompositeBezierLen(CP[i])
                        uniqueObjectsBezier1[i] = length
                if fnr == 2 :
    #                PlotBezier(CP[i], '#bc0e13', 1, (0, (3, 5, 1, 5)))
                    if uniqueObjectsBezier2.get(i) == None and ObjectsBezier2.get(i) == None and equalObjectsBezier2.get(i) == None :
                        length = BezierCalculations.CompositeBezierLen(CP[i])
                        uniqueObjectsBezier2[i] = length
                if not (CP[i][0][0][0], CP[i][0][0][1]) in visited :
                    visited[(CP[i][0][0][0], CP[i][0][0][1])] = 1
                if not (CP[i][-1][3][0], CP[i][-1][3][1]) in visited :
                    visited[(CP[i][-1][3][0], CP[i][-1][3][1])] = 1

        visited[P] = 2
        return visited
    #Calculate LeastSquare value between two given lines
    def LeastSquare(B1, int1, lineSt1, B2, int2, lineSt2, CP1, CP2) :
        dist1 = 0
        dist2 = 0
        value = 0
        lenB1 = len(B1)
        lenB2 = len(B2)
        disttmp1 = 0
        disttmp2 = 0
        if lineSt1 == 1 :
            nr1 = lenB1 - 1
        else :
            nr1 = 0
        if lineSt2 == 1 :
            nr2 = lenB2 - 1
        else :
            nr2 = 0
        param1 = lineSt1
        param2 = lineSt2
        for i in range(11) :
            P1x, P1y, nr1, disttmp1, param1 = BezierCalculations.DistantPointOnBezier(dist1, nr1, B1, lineSt1, CP1, disttmp1, param1)
            P2x, P2y, nr2, disttmp2, param2 = BezierCalculations.DistantPointOnBezier(dist2, nr2, B2, lineSt2, CP2, disttmp2, param2)
            value += (Calculations.PointDist((P1x, P1y), (P2x, P2y)))**2
            dist1 += int1
            dist2 += int2
        t = Symbol('t')
        pos = int(round(lenB1 / 2))
        a = CP1[pos][0][0]
        b = CP1[pos][0][1]
        ax.annotate(value, xy = (a, b), xytext = (a, b), size = 'small', alpha = 0.8,
            color = '#562048', va = 'center', ha = 'center', family = 'monospace')
        return value


    def DiffAllBezier(P1, visited, line1, line2, Bezier1, Bezier2, CP1, CP2, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2, minDist, absMaxDist) :
        #Check if given point isn`t already visited
        if P1 in visited :
            if visited[P1] == 2 :
                return visited, absMaxDist

        #Get coordinates of the closest point in second file
        if len(Bezier2) > 0 :
            print("Closest Point")
            P2, min = BezierCalculations.ClosestPntBezier(P1, CP2, visited, minDist)
            print("Closest Point")
        else :
            min = -1
            visited = ObjectComparison.OneFileObject(P1, line1, Bezier1, CP1, visited, 1, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2)
            return visited, absMaxDist, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2
        #if the point is too far, then probably that is not the respective point
        if min >= 30 :
            print("min")
            visited = ObjectComparison.OneFileObject(P1, line1, Bezier1, CP1, visited, 1, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2)
            print("min")
            return visited, absMaxDist, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2
        #Find all lines with one of the enpoints matching the given point
        lenLine1 = len(line1)
        lenLine2 = len(line2)
        P1rndx = round(P1[0], 3)
        P1rndy = round(P1[1], 3)
        P2rndx = round(P2[0], 3)
        P2rndy = round(P2[1], 3)
        lenCP1 = len(CP1)
        lenCP2 = len(CP2)
        for i in range(lenCP1) :
            P1Ind = -1
            lineSt1 = -1
            t = Symbol('t')
            P1tmprndx0 = round(CP1[i][0][0][0], 3)
            P1tmprndy0 = round(CP1[i][0][0][1], 3)
            P1tmprndx1 = round(CP1[i][-1][3][0], 3)
            P1tmprndy1 = round(CP1[i][-1][3][1], 3)
            if P1rndx == P1tmprndx0 and P1rndy == P1tmprndy0 :
                if not (CP1[i][0][0][0], CP1[i][0][0][1]) in visited :
                    if not (CP1[i][-1][3][0], CP1[i][-1][3][1]) in visited :
                        P1Ind = i
                        lineSt1 = 0
                        visited[(CP1[i][-1][3][0], CP1[i][-1][3][1])] = 1
                    elif visited[(CP1[i][-1][3][0], CP1[i][-1][3][1])] == 1 :
                        P1Ind = i
                        lineSt1 = 0
                elif visited[(CP1[i][0][0][0], CP1[i][0][0][1])] == 1 :
                    if not (CP1[i][-1][3][0], CP1[i][-1][3][1]) in visited :
                        P1Ind = i
                        lineSt1 = 0
                        visited[(CP1[i][-1][3][0], CP1[i][-1][3][1])] = 1
                    elif visited[(CP1[i][-1][3][0], CP1[i][-1][3][1])] == 1 :
                        P1Ind = i
                        lineSt1 = 0
            elif P1rndx == P1tmprndx1 and P1rndy == P1tmprndy1 :
                if not (CP1[i][-1][3][0], CP1[i][-1][3][1]) in visited :
                    if not (CP1[i][0][0][0], CP1[i][0][0][1]) in visited :
                        P1Ind = i
                        lineSt1 = 1
                        visited[(CP1[i][0][0][0], CP1[i][0][0][1])] = 1
                    elif visited[(CP1[i][0][0][0], CP1[i][0][0][1])] == 1 :
                        P1Ind = i
                        lineSt1 = 1
                elif visited[(CP1[i][-1][3][0], CP1[i][-1][3][1])] == 1 :
                    if not (CP1[i][0][0][0], CP1[i][0][0][1]) in visited :
                        P1Ind = i
                        lineSt1 = 1
                        visited[(CP1[i][0][0][0], CP1[i][0][0][1])] = 1
                    elif visited[(CP1[i][0][0][0], CP1[i][0][0][1])] == 1 :
                        P1Ind = i
                        lineSt1 = 1
            else:
                continue
            #This is the case wher the other line end has been "fully visited" or visited[] = 2
            if P1Ind == -1 :
                continue
            a1 = Calculations.LineSlope((CP1[i][0][0][0], CP1[i][0][0][1]), (CP1[i][-1][3][0], CP1[i][-1][3][1]))
            lineSt2 = -1
            P2Ind = -1
            minEndDist = -1

            for j in range(lenCP2):
                lineSt2tmp = -1
                P2tmprndx0 = round(CP2[j][0][0][0], 3)
                P2tmprndy0 = round(CP2[j][0][0][1], 3)
                P2tmprndx1 = round(CP2[j][-1][3][0], 3)
                P2tmprndy1 = round(CP2[j][-1][3][1], 3)
                if P2rndx == P2tmprndx0 and P2rndy == P2tmprndy0 :
                    if not (CP2[j][0][0][0], CP2[j][0][0][1]) in visited :
                        if not (CP2[j][-1][3][0], CP2[j][-1][3][1]) in visited :
                            lineSt2tmp = 0
                        elif visited[(CP2[j][-1][3][0], CP2[j][-1][3][1])] == 1 :
                            lineSt2tmp = 0
                    elif visited[(CP2[j][0][0][0], CP2[j][0][0][1])] == 1 :
                        if not (CP2[j][-1][3][0], CP2[j][-1][3][1]) in visited :
                            lineSt2tmp = 0
                        elif visited[(CP2[j][-1][3][0], CP2[j][-1][3][1])] == 1 :
                            lineSt2tmp = 0
                elif P2rndx == P2tmprndx1 and P2rndy == P2tmprndy1 :
                    if not (CP2[j][-1][3][0], CP2[j][-1][3][1]) in visited :
                        if not (CP2[j][0][0][0], CP2[j][0][0][1]) in visited :
                            lineSt2tmp = 1
                        elif visited[(CP2[j][0][0][0], CP2[j][0][0][1])] == 1 :
                            lineSt2tmp = 1
                    elif visited[(CP2[j][-1][3][0], CP2[j][-1][3][1])] == 1 :
                        if not (CP2[j][0][0][0], CP2[j][0][0][1]) in visited :
                            lineSt2tmp = 1
                        elif visited[(CP2[j][0][0][0], CP2[j][0][0][1])] == 1 :
                            lineSt2tmp = 1
                else:
                    continue
                #This is the case wher the other line end has been "fulli visited" or visited[] = 2
                if lineSt2tmp == -1 :
                    continue
                a2 = Calculations.LineSlope((CP2[j][0][0][0], CP2[j][0][0][1]), (CP2[j][-1][3][0], CP2[j][-1][3][1]))
                angtmp = round((a2 - a1) / (1 + abs(a1 * a2)), 5)
                angle = np.degrees(np.arctan(angtmp))
                if abs(angle) > 5 :
                    #if (i == lenLine1 - 1 and not j in ObjectsBezier2 and not j in equalObjectsBezier2) :

    #                    PlotBezier(CP2[j], '#6c9f92', 1, (0, (3, 5, 1, 5)))
                    #    visited[(Bezier2[j][-lineSt2tmp][0].subs(t, lineSt2tmp), Bezier2[j][-lineSt2tmp][1].subs(t, lineSt2tmp))] = 1
                    #    if not (Bezier2[j][- abs(-lineSt2tmp + 1)][0].subs(t, abs(lineSt2tmp - 1)), Bezier2[j][- abs(-lineSt2tmp + 1)][1].subs(t, abs(lineSt2tmp - 1))) in visited:
                    #        visited[Bezier2[j][- abs(-lineSt2tmp + 1)][0].subs(t, abs(lineSt2tmp - 1)), Bezier2[j][- abs(-lineSt2tmp + 1)][1].subs(t, abs(lineSt2tmp))] = 1
                    #    if uniqueObjectsBezier2.get(j) == None and ObjectsBezier2.get(j) == None:
                    #3        length = CompositeBezierLen(CP2[j])
                    #        uniqueObjectsBezier2[j] = length
                    continue
                if (lineSt1 != -1 and lineSt2tmp != -1) :
                    endDist = Calculations.PointDist((CP1[i][-((lineSt1 + 1) % 2)][((lineSt1 + 1) % 2) * 3][0], CP1[i][-((lineSt1 + 1) % 2)][((lineSt1 + 1) % 2) * 3][1]), (CP2[j][-((lineSt2tmp + 1) % 2)][((lineSt2tmp + 1) % 2) * 3][0], CP2[j][-((lineSt2tmp + 1) % 2)][((lineSt2tmp + 1) % 2) * 3][1]))
                    #->if (i == lenLine1 - 1 and not j in ObjectsBezier2 and not j in equalObjectsBezier2 and endDist >= 30) :
                    #->    visited, absMaxDist = ObjectComparison.DiffAllBezier(P1, visited, line1, line2, Bezier1, Bezier2, CP1, CP2, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, min, absMaxDist)
    #               #->     PlotBezier(CP2[j], '#6c9f92', 1, (0, (3, 5, 1, 5)))
                    #->    visited[(CP2[j][-lineSt2tmp][lineSt2tmp * 3][0], CP2[j][-lineSt2tmp][lineSt2tmp * 3][1])] = 1

                    #->    if not (CP1[i][-((lineSt1 + 1) % 2)][((lineSt1 + 1) % 2) * 3][0], CP1[i][-((lineSt1 + 1) % 2)][((lineSt1 + 1) % 2) * 3][1]) in visited :
                    #->        visited[(CP1[i][-((lineSt1 + 1) % 2)][((lineSt1 + 1) % 2) * 3][0], CP1[i][-((lineSt1 + 1) % 2)][((lineSt1 + 1) % 2) * 3][1])] = 1
                    #->    if uniqueObjectsBezier2.get(i) == None and ObjectsBezier2.get(i) == None:
                    #->        print("comp")
                    #->        length = BezierCalculations.CompositeBezierLen(CP2[i])
                    #->        print("comp")
                    #->        uniqueObjectsBezier2[i] = length
                    if (minEndDist == -1 and endDist < 30):
                        minEndDist = endDist
                        P2Ind = j
                        lineSt2 = lineSt2tmp
                        print("maxDist")
                        absDist, P1max, P2max = BezierCalculations.MaxBezierDist(Bezier1[i], Bezier2[j], CP1[i], CP2[j], lineSt1, lineSt2)
                        print("maxDist")
                    elif endDist <= minEndDist :
                    #    if (i == lenLine1 - 1 and not P2Ind in ObjectsBezier2 and not P2Ind in equalObjectsBezier2) :
    #                #        PlotBezier(CP2[P2Ind], '#6c9f92', 1, (0, (3, 5, 1, 5)), linewidth = 4)
                    #        visited[(Bezier2[P2Ind][-lineSt2][0].subs(t, lineSt2), Bezier2[P2Ind][-lineSt2][1].subs(t, lineSt2))] = 1
                    #        if not (Bezier2[P2Ind][-abs(-lineSt2 + 1)][0].subs(t, abs(lineSt2 - 1)), Bezier2[P2Ind][-abs(-lineSt2 + 1)][1].subs(t, abs(lineSt2 - 1))) :
                    #            visited[(Bezier2[P2Ind][-abs(-lineSt2 + 1)][0].subs(t, abs(lineSt2 - 1)), Bezier2[P2Ind][-abs(-lineSt2 + 1)][1].subs(t, abs(lineSt2 - 1)))] = 1
                    #        if uniqueObjectsBezier2.get(P2Ind) == None and ObjectsBezier2.get(P2Ind) == None :
                    #            length = CompositeBezierLen(CP2[P2Ind])
                    #            uniqueObjectsBezier2[P2Ind] = length
                        if endDist == minEndDist :
                            print("maxDist1")
                            absDisttmp, P1maxtmp, P2maxtmp = BezierCalculations.MaxBezierDist(Bezier1[i], Bezier2[j], CP1[i], CP2[j],  lineSt1, lineSt2)
                            print("maxDist1")
                            if absDisttmp < absDist :
                                minEndDist = endDist
                                P2Ind = j
                                lineSt2 = lineSt2tmp
                                absDist = absDisttmp
                                P1max = P1maxtmp
                                P2max = P2maxtmp
                        else :
                            minEndDist = endDist
                            P2Ind = j
                            lineSt2 = lineSt2tmp
                            print("maxDist3")
                            absDist, P1max, P2max = BezierCalculations.MaxBezierDist(Bezier1[i], Bezier2[j], CP1[i], CP2[j],  lineSt1, lineSt2)
                            print("maxDist3")


            if (P2Ind == -1) :
                visited, absMaxDist, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2 = ObjectComparison.DiffAllBezier(P1, visited, line1, line2, Bezier1, Bezier2, CP1, CP2, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier, ObjectsBezier1, ObjectsBezier2, min, absMaxDist)
    #            PlotBezier(CP1[i], '#055583', 1, (0, (3, 5, 1, 5)))
                t = Symbol('t')
    #            visited[(CP1[i][-lineSt1][lineSt1 * 3][0], CP1[i][-lineSt1][lineSt1 * 3][1])] = 2
                if not (CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][0], CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][1]) :
                    visited[(CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][0], CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][1])] = 1
                if uniqueObjectsBezier1.get(i) == None and ObjectsBezier1.get(i) == None :
                    print("Length")
                    length = BezierCalculations.CompositeBezierLen(CP1[i])
                    print("Length")
                    uniqueObjectsLine1[i] = length

            if (P2Ind != -1 and P1Ind != -1):
                if not (CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][0], CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][1]) in visited :
                    visited[CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][0], CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][1]] = 1
                if not (CP2[P2Ind][-abs(-lineSt2 + 1)][abs(-lineSt2 + 1) * 3][0], CP2[P2Ind][-abs(-lineSt2 + 1)][abs(-lineSt2 + 1) * 3][1]) :
                    visited[CP2[P2Ind][-abs(-lineSt2 + 1)][abs(-lineSt2 + 1) * 3][0], CP2[P2Ind][-abs(-lineSt2 + 1)][abs(-lineSt2 + 1) * 3][1]] = 1

                if not (CP1[i][-lineSt1][lineSt1 * 3][0], CP1[i][-lineSt1][lineSt1 * 3][1]) in visited :
                    visited[CP1[i][-lineSt1][lineSt1 * 3][0], CP1[i][-lineSt1][lineSt1 * 3][1]] = 1
                if not (CP2[P2Ind][-lineSt2][lineSt2 * 3][0], CP2[P2Ind][-lineSt2][lineSt2 * 3][1]) :
                    visited[CP2[P2Ind][-lineSt2][lineSt2 * 3][0], CP2[P2Ind][-lineSt2][lineSt2 * 3][1]] = 1
                print("Length2")
                len1 = BezierCalculations.CompositeBezierLen(CP1[i])
                len2 = BezierCalculations.CompositeBezierLen(CP2[P2Ind])
                print("Length2")
                a = ((P1max[0] + P2max[0]) / 2)
                b = ((P1max[1] + P2max[1]) / 2)
                ax.annotate(round(absDist / 10, 2), xy = (a, b), xytext = (a, b), alpha = 0.5, size = 'small')
                plt.plot([P1max[0], P2max[0]], [P1max[1], P2max[1]], color = 'black', alpha = 0.4)

                #equal objects will be drawn from file 1
                if absDist <= 0.1 :
                    if equalObjectsBezier.get(P1Ind) == None :
                        equalObjectsBezier[P1Ind] = len1
                    if equalObjectsBezier2.get(P2Ind) == None :
                        equalObjectsBezier2[P2Ind] = len2
                    continue

                if absDist > absMaxDist :
                    absMaxDist = absDist

                if ObjectsBezier1.get(i) == None :
                    ObjectsBezier1[i] = len1
                if ObjectsBezier2.get(P2Ind) == None :
                    ObjectsBezier2[P2Ind] = len2



                int1 = len1 / 10
                int2 = len2 / 10
                #print("leastSquare")
                #ObjectComparison.LeastSquare(Bezier1[i], int1, lineSt1, Bezier2[P2Ind], int2, lineSt2, CP1[i], CP2[P2Ind])
                #print("leastSquare")
    #            BezierDiff(Bezier1[i], int1, lineSt1, Bezier2[P2Ind], int2, lineSt2, CP1, CP2)
        visited[P2] = 1
        visited[P1] = 2
        return visited, absMaxDist, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2



class Files(ObjectComparison, Plot):
    def __init(files) :
        pass


    def Compare(points1, line1, Bezier1, CP1, points2, line2, Bezier2, CP2) :


        #Min distance from given point to Bezier curve
        #given: cnt - nr of Bezier curve segment
        #       B - Bezier curve on which we need to find point on
        #       P - given point from which the closest line will be calculated
        #       minDist - that has been found already
        #       len - the length of B array
        #       minPtmp - just some point, that we dont need to create every time the function is called
        #       minP - coordinates of the closest point on Bezier curve to P
        #def MinDistTanBezier(cnt, B, P, minDist, len, minPtmp, minP) :
        #    if cnt > len - 1 :
        #        return len - 1, minP, minDist
        #    dist, minPtmp = BezierMinDist(B[cnt][0], B[cnt][1], P)
        #    if minDist >= dist :
        #        minDist = dist
        #        minP = minPtmp
        #        cnt += 1
        #        return MinDistTanBezier(cnt, B, P, minDist, len)
        #    elif dist > minDist :
        #        cnt -= 1
        #        return cnt, minP, minDist
        #    return cnt, minP, minDist

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

        equalObjectsLine2 = {}
        equalObjectsBezier2 = {}

        lenLine1 = len(line1)
        lenLine2 = len(line2)
        lenBezier1 = len(Bezier1)
        lenBezier2 = len(Bezier2)

        visited = {}
        #    elif visited[i[1]] == 1:
        #        visited = DiffAll(i[1], visited, line1, line2, Bezier1, Bezier2, CP1, CP2, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2)

        absMaxDist = 0
        t = Symbol('t')
        for i in CP1 :
            if not (i[0][0][0], i[0][0][1]) in visited :
                visited, absMaxDist, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2 = ObjectComparison.DiffAllBezier((i[0][0][0], i[0][0][1]), visited, line1, line2, Bezier1, Bezier2, CP1, CP2,uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2, -1, absMaxDist)
            elif visited[(i[0][0][0], i[0][0][1])] == 1 :
                visited, absMaxDist, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2 = ObjectComparison.DiffAllBezier((i[0][0][0], i[0][0][1]), visited, line1, line2, Bezier1, Bezier2, CP1, CP2,uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2, -1, absMaxDist)

        #for i in Bezier1:
        #    if (i[-1][0].subs(t, 1), i[-1][1].subs(t, 1)) in visited :
        #        if visited[(i[-1][0].subs(t, 1), i[-1][1].subs(t, 1))] == 1 :
        #            visited = DiffAll((i[-1][0].subs(t, 1), i[-1][1].subs(t, 1)), visited, line1, line2, Bezier1, Bezier2, CP1, CP2,uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2)

        for i in range(lenLine2) :
            if not i in ObjectsLine2 and not i in equalObjectsLine2 :
                visited = ObjectComparison.OneFileObject(line2[i][0], line2, Bezier2, CP2, visited, 2, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2)

        for i in range(lenBezier2) :
            if not i in ObjectsBezier2 and not i in equalObjectsBezier2:
                visited = ObjectComparison.OneFileObject((CP2[i][0][0][0], CP2[i][0][0][1]), line2, Bezier2, CP2, visited, 2, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2)

        return visited, absMaxDist, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2, uniqueObjectsBezier1, uniqueObjectsBezier2

        #new color of the first file : #20456d
    def PlotFiles(absMaxDist, Bezier1, CP1, Bezier2, CP2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2, uniqueObjectsBezier1, uniqueObjectsBezier2) :
        nreqf = 0.0
        nr1f = 1.0
        nr2f = 2.0
        nrf = 3.0
        t = Symbol('t')
        print()
        lenBezier1 = len(Bezier1)
        for i in range(lenBezier1) :
            if i in equalObjectsBezier :
                Plot.PlotBezier(CP1[i], '#60635d', 0.3, 'solid')
                nreqf += 0.01
                pos = int(round(len(Bezier1[i]) / 5))
                a = Bezier1[i][pos][0].subs(t, 0.5)
                b = Bezier1[i][pos][1].subs(t, 0.5)
                ax.annotate(round(nreqf, 2), xy = (a, b), xytext = (a, b), family = 'monospace',
                    color = '#60635d', ha = 'center', va = 'center', weight = 'bold')
                print(round(nreqf, 2), ' - ', equalObjectsBezier[i] / 10, 'cm')
            if i in ObjectsBezier1 :
                Plot.PlotBezier(CP1[i], firstFileCol, 0.5, 'solid')
                nr1f += 0.01
                pos = int(round(len(Bezier1[i]) / 5))
                a = Bezier1[i][pos][0].subs(t, 0.5)
                b = Bezier1[i][pos][1].subs(t, 0.5)
                ax.annotate(round(nr1f, 2), xy = (a, b), xytext = (a, b), family = 'monospace',
                    color = firstFileCol, ha = 'center', va = 'center', weight = 'bold')
                print(round(nr1f, 2), ' - ', ObjectsBezier1[i] / 10, 'cm')
            if i in uniqueObjectsBezier1 :
                Plot.PlotBezier(CP1[i], firstFileCol, 0.5, (0, (3, 5, 1, 5)))
                nr1f += 0.01
                pos = int(round(len(Bezier1[i]) / 2))
                a = Bezier1[i][pos][0].subs(t, 0.5)
                b = Bezier1[i][pos][1].subs(t, 0.5)
                ax.annotate(round(nr1f, 2), xy = (a, b), xytext = (a, b), family = 'monospace',
                    color = firstFileCol, ha = 'center', va = 'center', weight = 'bold')
                print(round(nr1f, 2), ' - ', uniqueObjectsBezier1[i] / 10, 'cm')

        print()
        lenBezier2 = len(Bezier2)
        for i in range(lenBezier2) :
            if i in ObjectsBezier2 :
                Plot.PlotBezier(CP2[i], secondFileCol, 0.5, 'solid')
                nr2f += 0.01
                pos = int(round(len(Bezier2[i]) / 4))
                a = Bezier2[i][pos][0].subs(t, 0.5)
                b = Bezier2[i][pos][1].subs(t, 0.5)
                ax.annotate(round(nr2f, 2), xy = (a, b), xytext = (a, b), family = 'monospace',
                    color = '#59080a', ha = 'center', va = 'center', weight = 'bold')
                print(round(nr2f, 2), ' - ', ObjectsBezier2[i] / 10, 'cm')
            if i in uniqueObjectsBezier2 :
                Plot.PlotBezier(CP2[i], secondFileCol, 0.5, (0, (3, 5, 1, 5)))
                nr2f += 0.01
                pos = int(round(len(Bezier2[i]) / 2))
                a = Bezier2[i][pos][0].subs(t, 0.5)
                b = Bezier2[i][pos][1].subs(t, 0.5)
                ax.annotate(round(nr2f, 2), xy = (a, b), xytext = (a, b), family = 'monospace',
                    color = '#59080a', ha = 'center', va = 'center', weight = 'bold')
                print(round(nr2f, 2), ' - ', uniqueObjectsBezier2[i] / 10, 'cm')

        t = Symbol('t')
        for i in CP1 :
            Plot.PlotBezier(i, firstFileCol, 0.2, 'solid')
            plt.plot(i[0][0][0], i[0][0][1], 'o', ms = 1, color = firstFileCol, alpha = 0.7)
            plt.plot(i[-1][3][0], i[-1][3][1], 'o', ms = 1, color = firstFileCol, alpha = 0.7)
        for i in CP2 :
            Plot.PlotBezier(i, secondFileCol, 0.2, 'solid')
            plt.plot(i[0][0][0], i[0][0][1], 'o', ms = 1, color = secondFileCol, alpha = 0.7)
            plt.plot(i[-1][3][0], i[-1][3][1], 'o', ms = 1, color = secondFileCol, alpha = 0.7)

        print()
        print("MAX Distance")
        print(absMaxDist / 10)

        return

    def WriteLog(nr, fileName1, fileName2, absMaxDist, maxOffset, uniqueObjectsBezier1, uniqueObjectsBezier2, ObjectsBezier1) :
        resultlog.info('{0}--- {1}   {2}'.format(nr, fileName1, fileName2))
        if absMaxDist / 10 > maxOffset :
            resultlog.info('  Max Distance: {0} cm'.format(absMaxDist / 10))
            resultlog.info('  Unique Object Sum: {0} + {1}'.format(len(uniqueObjectsBezier1), len(uniqueObjectsBezier2)))
            resultlog.info('  Mismatching Objects Sum: {0}'.format(len(ObjectsBezier1)))
            return 0
        else:
            resultlog.info('OK')
            return 1


def GetObjects(fileName1, fileName2) :
    print("FILE1")
    if fileName1[-3] == 'd' and fileName1[-2] == 'x' and fileName1[-1] == 'f' :
        file1 = File1(fileName1)
        points1, line1, Bezier1, CP1 = file1.CalculateObjects()
        fileName1 = fileName1.replace('.dxf', '')
    else :
        file1 = File2(fileName1)
        Files()
        points1, line1, Bezier1, CP1 = file1.CalculateObjects()
        fileName1 = fileName1.replace('.svg', '')
    print("FILE1")
    print("FILE2")
    if fileName2[-3] == 'd' and fileName2[-2] == 'x' and fileName2[-1] == 'f' :
        file2 = File1(fileName2)
        points2, line2, Bezier2, CP2 = file2.CalculateObjects()
        fileName2 = fileName1.replace('.dxf', '')
    else :
        file2 = File2(fileName2)
        points2, line2, Bezier2, CP2 = file2.CalculateObjects()
        fileName2 = fileName1.replace('.svg', '')
    print("FILE2")

    visited, absMaxDist, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2, uniqueObjectsBezier1, uniqueObjectsBezier2 = Files.Compare(points1, line1, Bezier1, CP1, points2, line2, Bezier2, CP2)
    Files.PlotFiles(absMaxDist, Bezier1, CP1, Bezier2, CP2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2, uniqueObjectsBezier1, uniqueObjectsBezier2)
    return absMaxDist, uniqueObjectsBezier1, uniqueObjectsBezier2, ObjectsBezier1, fileName1, fileName2

def MyForm(x, p):
    return x/10


if not os.path.exists(filename1):
    infolog.error('ERROR:     {0} not found'.format(filename1))
if not os.path.exists(filename2):
    infolog.error('ERROR:     {0} not found'.format(filename2))


nr = 1
isOkCnt = 0
mismatchFiles = {}

if (filename1.endswith(".dxf") or filename1.endswith(".svg")) and (filename2.endswith(".dxf") or filename2.endswith(".svg")) :
    fileName1 = filename1
    fileName2 = filename2
    #Graph formatting
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

    absMaxDist, uniqueObjectsBezier1, uniqueObjectsBezier2, ObjectsBezier1, fileName1, fileName2 = GetObjects(fileName1, fileName2)
    isOk = Files.WriteLog(nr, fileName1, fileName2, absMaxDist, maxOffset, uniqueObjectsBezier1, uniqueObjectsBezier2, ObjectsBezier1)
    if showPlot == True :
        plt.show()
elif not (filename1.endswith(".dxf") or filename1.endswith(".svg")) or (filename2.endswith(".dxf") or filename2.endswith(".svg")) :
    unfoundcnt = 0
    unfoundPair = []
    for filename in os.listdir(filename1) :
        if filename.endswith(".dxf") or filename.endswith(".svg") :
            fileName1 = filename1 + '/' + filename
            existsFile2 = False
            for filenametmp in os.listdir(filename2) :
                if (filenametmp.endswith(".dxf") or filenametmp.endswith(".svg")) and (filenametmp.startswith(filename)):

                    #Graph formatting
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

                    fileName2 = filename2 + '/' + filenametmp
                    absMaxDist, uniqueObjectsBezier1, uniqueObjectsBezier2, ObjectsBezier1, fileName1, fileName2 = GetObjects(fileName1, fileName2)
                    isOk = Files.WriteLog(nr, fileName1, fileName2, absMaxDist, maxOffset, uniqueObjectsBezier1, uniqueObjectsBezier2, ObjectsBezier1)
                    isOkCnt += isOk
                    if isOk == 0:
                        mismatchFiles[(filename, filenametmp)] = absMaxDist / 10
                    nr += 1
                    existsFile2 = True
                    if showPlot == True :
                        plt.show()
            if existsFile2 == False :
                unfoundcnt += 1
                unfoundPair.append(filename)
                resultlog.error('ERROR:     matching file not found for {0}'.format(filename))
        else :
            continue
else :
    resultlog.error('ERROR:     given file types don´t match')


if os.path.exists(fileName1 + "_" + fileName2 + ".log") :
    os.remove(fileName1 + "_" + fileName2 + ".log")

handler = logging.FileHandler(fileName1 + "_" + fileName2 + ".log")
handler.setFormatter(logging.Formatter('%(message)s'))
resultlog = logging.getLogger("resultlog")
resultlog.setLevel(logging.INFO)
resultlog.addHandler(handler)


resultlog.info('Number of Matching files (max offset: {2} cm):      {0} / {1}'.format(isOkCnt, nr - 1, maxOffset))
if isOkCnt != nr - 1 :
    resultlog.info('')
    resultlog.info('files that don`t match:'.format(fileName2))
    cnt = 0
    for filenames in mismatchFiles :
        cnt += 1
        resultlog.info('{0}. max offset: {1} cm   filename: {2}'.format(cnt, mismatchFiles[filenames], filenames))
if unfoundcnt > 0:
    resultlog.info('')
    resultlog.info("{0} files (in {1}) with unfound comparable file: ".format(unfoundcnt, filename1))
    cnt = 0
    for filename in unfoundPair :
        cnt += 1
        resultlog.info('    {0}    {1}'.format(cnt, filename))

stat = open(statFile, "r")
textLines = stat.readlines()
linecnt = len(textLines)

if len(textLines) == 0:
    stat.close()
    stat = open(statFile, "w")
    stat.write("{0} / {1} compared files match (max offset {2} cm)".format(isOkCnt, nr - 1, maxOffset))
else :
    stat.close()
    difOffset = True
    print(len(textLines))
    stat = fileinput.input(statFile, inplace = True)
    for n, textline in enumerate(stat, start = 0):
        if len(textline) < 5 :
            continue
        statOk = ''
        i = 0
        while textline[i] != ' ':
            statOk += textline[i]
            i += 1
        statOk = int(statOk)

        statcnt = ''
        i += 3
        while textline[i] != ' ':
            statcnt += textline[i]
            i += 1
        statcnt = int(statcnt)

        pastOffset = ''
        i = len(textline) - 2
        tmp = "aa"
        while tmp != "cm" :
            tmp = textline[i] + tmp[0]
            i -= 1
        i -= 1
        while textline[i] != ' ':
            pastOffset = textline[i] + pastOffset
            i -= 1
        pastOffset = round(float(pastOffset), 2)
        if pastOffset == maxOffset :
            difOffset = False
        #stat.close()
        #stat = open(statFile, "w")
            if n == 0 :
                print("{0} / {1} compared files match (max offset {2} cm)   {3}% OK".format(statOk + isOkCnt, statcnt + nr - 1, maxOffset, (100 / (statcnt + nr - 1) * (statOk + isOkCnt))), end = '')
            else :
                print("\n{0} / {1} compared files match (max offset {2} cm) {3}% OK".format(statOk + isOkCnt, statcnt + nr - 1, maxOffset, (100 / (statcnt + nr - 1) * (statOk + isOkCnt))), end = '')
        #textline in fileinput.input(statFile, inplace = True)
            #textline = ("{0} / {1} compared files match (max offset {2} cm)\n".format(statOk + isOkCnt, statcnt + nr - 1, maxOffset)) + textline.rstrip('\n')
            continue

        if pastOffset > maxOffset and difOffset == True:
            difOffset = False
            if n == 0 :
                print("{0} / {1} compared files match (max offset {2} cm)   {3}% OK".format(isOkCnt, nr - 1, maxOffset, 100 / (nr - 1) * isOkCnt))
                print(textline.replace('\n', ''), end = '')
            else :
                print("\n{0} / {1} compared files match (max offset {2} cm) {3}% OK".format(isOkCnt, nr - 1, maxOffset, 100 / (nr - 1) * isOkCnt))
                print(textline.replace('\n', ''), end = '')
            #print(textline, end = '')
            continue
        if difOffset == True and n == linecnt - 1:
            #print(textline, end = '')
            if n == 0 :
                print(textline.replace('\n', ''))
                print("{0} / {1} compared files match (max offset {2} cm)   {3}% OK".format(isOkCnt, nr - 1, maxOffset, 100 / (nr - 1) * isOkCnt), end = '')
            else :
                print("\n{0}".format(textline.replace('\n', '')))
                print("{0} / {1} compared files match (max offset {2} cm)   {3}% OK".format(isOkCnt, nr - 1, maxOffset, 100 / (nr - 1) * isOkCnt), end = '')
            continue
        if n == 0 :
            print(textline.replace('\n', ''), end = '')
        else :
            print("\n{0}".format(textline.replace('\n', '')), end = '')
stat.close()
