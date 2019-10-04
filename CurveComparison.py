import dxfgrabber
import numpy
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
from pathlib import Path

#toCython
import pyximport; pyximport.install(language_level=3, pyimport=True)

#modules
import calculations
import calculateBezier
import transformation
import bezierCalculations
import DXF
import SVG

#Read parameters from cmd file
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

#remove extention for further
if filename1.endswith(".dxf") :
    fileName1 = filename1.replace(".dxf", '')
elif filename1.endswith(".svg") :
    fileName1 = filename1.replace(".svg", '')
else :
    fileName1 = filename1
if filename2.endswith(".dxf") :
    fileName2 = filename2.replace(".dxf", '')
elif filename2.endswith(".svg") :
    fileName2 = filename2.replace(".svg", '')
else :
    fileName2 = filename2
#write new log file for given files
if os.path.exists(fileName1 + "_" + fileName2 + ".log") :
    os.remove(fileName1 + "_" + fileName2 + ".log")

handler = logging.FileHandler(fileName1 + "_" + fileName2 + ".log")
handler.setFormatter(logging.Formatter('%(message)s'))
resultlog = logging.getLogger("resultlog")
resultlog.setLevel(logging.INFO)
resultlog.addHandler(handler)

#set color for file plotting
firstFileCol = '#000ec7'
secondFileCol = '#bc0e13'

class Plot() :
    def __init__(plot) :
        pass

    #plots Bezier lines
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

#Comparison of objects
#we will go through all objects in both files by picking the first object in first files
#then in second file we will find the most similar one (by type and coordinates)
#we also know that some of the objects are connected which we are going to use to specify two 'equal' objects in both files
#to determine how similar are both files we will calculate max distance between respective objects where:
#   *points - distance between points (but this will not be calculated because for our target it`s not important)
#   lines - distance between end-points
#   Bezier curves - distance for set of points on Bezier (with t parameter change of 0.1)
#       for Bezier we need to iterate both lines because one of them can end before the other one and actual max length can be larger

#C - C[composite line nr][controlpoint nr][x, y], similar for Bezier curve`s list

#for second type of comparison we will compare respective Bezier curves by getting distance between points with equal parameter value
class ObjectComparison() :
    def __init__(comp) :
        pass

    #Adds objects to unique object type set (objects that are found in only one file)
    def OneFileObject(P, line, Bezier, CP, visited, fnr, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2) :
        lenCP = len(CP)
        for i in range(lenCP) :
            if (P == (CP[i][0][0][0], CP[i][0][0][1]) or (P == (CP[i][-1][3][0], CP[i][-1][3][1]))) :
                if fnr == 1 :
    #                PlotBezier(CP[i], '#055583', 1, (0, (3, 5, 1, 5)))
                    if uniqueObjectsBezier1.get(i) == None and ObjectsBezier1.get(i) == None and equalObjectsBezier.get(i) == None :
                        length = bezierCalculations.CompositeBezierLen(CP[i])
                        uniqueObjectsBezier1[i] = length
                if fnr == 2 :
    #                PlotBezier(CP[i], '#bc0e13', 1, (0, (3, 5, 1, 5)))
                    if uniqueObjectsBezier2.get(i) == None and ObjectsBezier2.get(i) == None and equalObjectsBezier2.get(i) == None :
                        length = bezierCalculations.CompositeBezierLen(CP[i])
                        uniqueObjectsBezier2[i] = length
                if not (CP[i][0][0][0], CP[i][0][0][1]) in visited :
                    visited[(CP[i][0][0][0], CP[i][0][0][1])] = 1
                if not (CP[i][-1][3][0], CP[i][-1][3][1]) in visited :
                    visited[(CP[i][-1][3][0], CP[i][-1][3][1])] = 1

        visited[P] = 2
        return visited
    #Calculate squared distance sum value between two given lines (there is a line conecting both Bezier curves in each parameter growth by set value (int1(2)))
    def SquaredDist(B1, int1, lineSt1, B2, int2, lineSt2, CP1, CP2) :
        dist1 = 0
        dist2 = 0
        value = 0
        lenB1 = len(B1)
        lenB2 = len(B2)
        disttmp1 = 0
        disttmp2 = 0
        #take the respective endpoints of Bezier curves to start growing parameter value from the same endpoints
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
        #add squared distance to squared distance sum(value) in each int1(2) step
        for i in range(11) :
            P1x, P1y, nr1, disttmp1, param1 = bezierCalculations.DistantPointOnBezier(dist1, nr1, B1, lineSt1, CP1, disttmp1, param1)
            P2x, P2y, nr2, disttmp2, param2 = bezierCalculations.DistantPointOnBezier(dist2, nr2, B2, lineSt2, CP2, disttmp2, param2)
            value += (calculations.PointDist((P1x, P1y), (P2x, P2y)))**2
            dist1 += int1
            dist2 += int2
        t = Symbol('t')
        pos = int(round(lenB1 / 2))
        a = CP1[pos][0][0]
        b = CP1[pos][0][1]
        ax.annotate(value, xy = (a, b), xytext = (a, b), size = 'small', alpha = 0.8,
            color = '#562048', va = 'center', ha = 'center', family = 'monospace')
        return value

    #Find diferences between given files
    def DiffAllBezier(P1, visited, line1, line2, Bezier1, Bezier2, CP1, CP2, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2, minDist, absMaxDist, showPlot) :
        #Check if given point isn`t already visited
        if P1 in visited :
            if visited[P1] == 2 :
                return visited, absMaxDist

        #Get coordinates of the closest point in second file
        if len(Bezier2) > 0 :
            print("Closest Point")
            P2, min = bezierCalculations.ClosestPntBezier(P1, CP2, visited, minDist)
            print("Closest Point")
        else :
            min = -1
            visited = ObjectComparison.OneFileObject(P1, line1, Bezier1, CP1, visited, 1, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2)
            return visited, absMaxDist, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2
        #if the point is too far, then probably that is not the respective point (now distance = 30 mm)
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
            #we round point coordinates because there may be minimal(unimportant) offset
            P1tmprndx0 = round(CP1[i][0][0][0], 3)
            P1tmprndy0 = round(CP1[i][0][0][1], 3)
            P1tmprndx1 = round(CP1[i][-1][3][0], 3)
            P1tmprndy1 = round(CP1[i][-1][3][1], 3)
            if P1rndx == P1tmprndx0 and P1rndy == P1tmprndy0 :
                #we have to check if the point is already visited (visited[] = 1) (if so then skip) and if it matches chosen point
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
            #This is the case where the line end has been visited[] = 2
            if P1Ind == -1 :
                continue

            a1 = calculations.LineSlope((CP1[i][0][0][0], CP1[i][0][0][1]), (CP1[i][-1][3][0], CP1[i][-1][3][1]))
            lineSt2 = -1
            P2Ind = -1
            minEndDist = -1
            for j in range(lenCP2):
                lineSt2tmp = -1
                #we round point coordinates because there may be minimal(unimportant) offset
                P2tmprndx0 = round(CP2[j][0][0][0], 3)
                P2tmprndy0 = round(CP2[j][0][0][1], 3)
                P2tmprndx1 = round(CP2[j][-1][3][0], 3)
                P2tmprndy1 = round(CP2[j][-1][3][1], 3)
                if P2rndx == P2tmprndx0 and P2rndy == P2tmprndy0 :
                    #we have to check if the point is already visited (visited[] = 1) (if so then skip) and if it matches chosen point
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
                #This is the case where the line end has been visited[] = 2
                if lineSt2tmp == -1 :
                    continue
                #Compare line`s that connects curves endpoints slopes, and check if it doesn`t exceed max offset angle and may still be seen as "equal"
                a2 = calculations.LineSlope((CP2[j][0][0][0], CP2[j][0][0][1]), (CP2[j][-1][3][0], CP2[j][-1][3][1]))
                angtmp = round((a2 - a1) / (1 + abs(a1 * a2)), 5)
                angle = numpy.degrees(numpy.arctan(angtmp))
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
                #if both lines, which are considered equal are found
                if (lineSt1 != -1 and lineSt2tmp != -1) :
                    #calculate distance between endpoints
                    endDist = calculations.PointDist((CP1[i][-((lineSt1 + 1) % 2)][((lineSt1 + 1) % 2) * 3][0], CP1[i][-((lineSt1 + 1) % 2)][((lineSt1 + 1) % 2) * 3][1]), (CP2[j][-((lineSt2tmp + 1) % 2)][((lineSt2tmp + 1) % 2) * 3][0], CP2[j][-((lineSt2tmp + 1) % 2)][((lineSt2tmp + 1) % 2) * 3][1]))
                    #->if (i == lenLine1 - 1 and not j in ObjectsBezier2 and not j in equalObjectsBezier2 and endDist >= 30) :
                    #->    visited, absMaxDist = ObjectComparison.DiffAllBezier(P1, visited, line1, line2, Bezier1, Bezier2, CP1, CP2, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, min, absMaxDist)
    #               #->     PlotBezier(CP2[j], '#6c9f92', 1, (0, (3, 5, 1, 5)))
                    #->    visited[(CP2[j][-lineSt2tmp][lineSt2tmp * 3][0], CP2[j][-lineSt2tmp][lineSt2tmp * 3][1])] = 1

                    #->    if not (CP1[i][-((lineSt1 + 1) % 2)][((lineSt1 + 1) % 2) * 3][0], CP1[i][-((lineSt1 + 1) % 2)][((lineSt1 + 1) % 2) * 3][1]) in visited :
                    #->        visited[(CP1[i][-((lineSt1 + 1) % 2)][((lineSt1 + 1) % 2) * 3][0], CP1[i][-((lineSt1 + 1) % 2)][((lineSt1 + 1) % 2) * 3][1])] = 1
                    #->    if uniqueObjectsBezier2.get(i) == None and ObjectsBezier2.get(i) == None:
                    #->        print("comp")
                    #->        length = bezierCalculations.CompositeBezierLen(CP2[i])
                    #->        print("comp")
                    #->        uniqueObjectsBezier2[i] = length
                    #check if distance is considered possible for equal objects
                    if (minEndDist == -1 and endDist < 30):
                        minEndDist = endDist
                        P2Ind = j
                        lineSt2 = lineSt2tmp
                        #calculate max distance between found objects
                        print("maxDist")
                        absDist, P1max, P2max = bezierCalculations.MaxBezierDist(Bezier1[i], Bezier2[j], CP1[i], CP2[j], lineSt1, lineSt2)
                        print("maxDist")
                    #check if endpoint distance is smaller than endpoint distance between previous respective objects to grow possibility that we have found the most "equal" objects between given files
                    elif endDist <= minEndDist :
                    #    if (i == lenLine1 - 1 and not P2Ind in ObjectsBezier2 and not P2Ind in equalObjectsBezier2) :
    #                #        PlotBezier(CP2[P2Ind], '#6c9f92', 1, (0, (3, 5, 1, 5)), linewidth = 4)
                    #        visited[(Bezier2[P2Ind][-lineSt2][0].subs(t, lineSt2), Bezier2[P2Ind][-lineSt2][1].subs(t, lineSt2))] = 1
                    #        if not (Bezier2[P2Ind][-abs(-lineSt2 + 1)][0].subs(t, abs(lineSt2 - 1)), Bezier2[P2Ind][-abs(-lineSt2 + 1)][1].subs(t, abs(lineSt2 - 1))) :
                    #            visited[(Bezier2[P2Ind][-abs(-lineSt2 + 1)][0].subs(t, abs(lineSt2 - 1)), Bezier2[P2Ind][-abs(-lineSt2 + 1)][1].subs(t, abs(lineSt2 - 1)))] = 1
                    #        if uniqueObjectsBezier2.get(P2Ind) == None and ObjectsBezier2.get(P2Ind) == None :
                    #            length = CompositeBezierLen(CP2[P2Ind])
                    #            uniqueObjectsBezier2[P2Ind] = length
                        #if endpoint distances are equal, then lets check between whitch files absolute distance between objects is smaller, therefore we choose the "slosest" objects
                        if endDist == minEndDist :
                            print("maxDist1")
                            absDisttmp, P1maxtmp, P2maxtmp = bezierCalculations.MaxBezierDist(Bezier1[i], Bezier2[j], CP1[i], CP2[j],  lineSt1, lineSt2)
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
                            absDist, P1max, P2max = bezierCalculations.MaxBezierDist(Bezier1[i], Bezier2[j], CP1[i], CP2[j],  lineSt1, lineSt2)
                            print("maxDist3")
            #condition where respective equal object is not found in second file
            if (P2Ind == -1) :

                visited, absMaxDist, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2 = ObjectComparison.DiffAllBezier(P1, visited, line1, line2, Bezier1, Bezier2, CP1, CP2, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier, ObjectsBezier1, ObjectsBezier2, min, absMaxDist, showPlot)
                t = Symbol('t')
    #            visited[(CP1[i][-lineSt1][lineSt1 * 3][0], CP1[i][-lineSt1][lineSt1 * 3][1])] = 2
                if not (CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][0], CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][1]) in visited :
                    visited[(CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][0], CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][1])] = 1
                #add length to curve in library
                if uniqueObjectsBezier1.get(i) == None and ObjectsBezier1.get(i) == None :
                    print("Length")
                    length = bezierCalculations.CompositeBezierLen(CP1[i])
                    print("Length")
                    uniqueObjectsLine1[i] = length
            #condition where respective equal objects are found
            if (P2Ind != -1 and P1Ind != -1):
                #add points to visited library
                if not (CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][0], CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][1]) in visited :
                    visited[CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][0], CP1[i][-abs(-lineSt1 + 1)][abs(-lineSt1 + 1) * 3][1]] = 1
                if not (CP2[P2Ind][-abs(-lineSt2 + 1)][abs(-lineSt2 + 1) * 3][0], CP2[P2Ind][-abs(-lineSt2 + 1)][abs(-lineSt2 + 1) * 3][1]) in visited :
                    visited[CP2[P2Ind][-abs(-lineSt2 + 1)][abs(-lineSt2 + 1) * 3][0], CP2[P2Ind][-abs(-lineSt2 + 1)][abs(-lineSt2 + 1) * 3][1]] = 1
                #add the other endpoints to visited library
                if not (CP1[i][-lineSt1][lineSt1 * 3][0], CP1[i][-lineSt1][lineSt1 * 3][1]) in visited :
                    visited[CP1[i][-lineSt1][lineSt1 * 3][0], CP1[i][-lineSt1][lineSt1 * 3][1]] = 1
                if not (CP2[P2Ind][-lineSt2][lineSt2 * 3][0], CP2[P2Ind][-lineSt2][lineSt2 * 3][1]) in visited:
                    visited[CP2[P2Ind][-lineSt2][lineSt2 * 3][0], CP2[P2Ind][-lineSt2][lineSt2 * 3][1]] = 1
                print("Length2")
                len1 = bezierCalculations.CompositeBezierLen(CP1[i])
                len2 = bezierCalculations.CompositeBezierLen(CP2[P2Ind])
                print("Length2")
                if showPlot == True :
                    a = ((P1max[0] + P2max[0]) / 2)
                    b = ((P1max[1] + P2max[1]) / 2)
                    ax.annotate(round(absDist / 10, 2), xy = (a, b), xytext = (a, b), alpha = 0.5, size = 'small')
                    plt.plot([P1max[0], P2max[0]], [P1max[1], P2max[1]], color = 'black', alpha = 0.4)

                #if max offset is less or equal to 0.1, then we consider objects as one (equal)
                #equal objects will be drawn from file 1
                if absDist <= 0.1 :
                    if equalObjectsBezier.get(P1Ind) == None :
                        equalObjectsBezier[P1Ind] = len1
                    if equalObjectsBezier2.get(P2Ind) == None :
                        equalObjectsBezier2[P2Ind] = len2
                    continue
                #if max offset is larger than 0.1, then we consider respective objects different
                if absDist > absMaxDist :
                    absMaxDist = absDist
                if ObjectsBezier1.get(i) == None :
                    ObjectsBezier1[i] = len1
                if ObjectsBezier2.get(P2Ind) == None :
                    ObjectsBezier2[P2Ind] = len2


                #int1 = len1 / 10
                #int2 = len2 / 10
                #print("SquaredDist")
                #ObjectComparison.SquaredDist(Bezier1[i], int1, lineSt1, Bezier2[P2Ind], int2, lineSt2, CP1[i], CP2[P2Ind])
                #print("SquaredDist")
    #            BezierDiff(Bezier1[i], int1, lineSt1, Bezier2[P2Ind], int2, lineSt2, CP1, CP2)
        #set found points as visited
        visited[P2] = 1
        visited[P1] = 2
        return visited, absMaxDist, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2



class Files(ObjectComparison, Plot):
    def __init(files) :
        pass

    #function that calls which objects to compare
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
        absMaxDist = 0
        t = Symbol('t')
        #check elements from 1. file that are not visited[] = 2
        for i in CP1 :
            if not (i[0][0][0], i[0][0][1]) in visited :
                visited, absMaxDist, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2 = ObjectComparison.DiffAllBezier((i[0][0][0], i[0][0][1]), visited, line1, line2, Bezier1, Bezier2, CP1, CP2,uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2, -1, absMaxDist, showPlot)
            elif visited[(i[0][0][0], i[0][0][1])] == 1 :
                visited, absMaxDist, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2 = ObjectComparison.DiffAllBezier((i[0][0][0], i[0][0][1]), visited, line1, line2, Bezier1, Bezier2, CP1, CP2,uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2, -1, absMaxDist, showPlot)

        #for i in Bezier1:
        #    if (i[-1][0].subs(t, 1), i[-1][1].subs(t, 1)) in visited :
        #        if visited[(i[-1][0].subs(t, 1), i[-1][1].subs(t, 1))] == 1 :
        #            visited = DiffAll((i[-1][0].subs(t, 1), i[-1][1].subs(t, 1)), visited, line1, line2, Bezier1, Bezier2, CP1, CP2,uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2)

        #for i in range(lenLine2) :
        #    if not i in ObjectsLine2 and not i in equalObjectsLine2 :
        #        visited = ObjectComparison.OneFileObject(line2[i][0], line2, Bezier2, CP2, visited, 2, uniqueObjectsLine1, uniqueObjectsLine2, uniqueObjectsBezier1, uniqueObjectsBezier2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2)

        #check the rest objects - the objects from 2. file
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
        if showPlot == True :
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

    def WriteLog(nr, f1, f2, absMaxDist, maxOffset, uniqueObjectsBezier1, uniqueObjectsBezier2, ObjectsBezier1) :
        resultlog.info('\nMax acceptable offset:   {0}'.format(maxOffset))
        resultlog.info('{0}--- {1}   {2}'.format(nr, f1, f2))
        if absMaxDist / 10 > maxOffset :
            resultlog.info('  Max Distance: {0} cm'.format(absMaxDist / 10))
            resultlog.info('  Unique Object Sum: {0} + {1}'.format(len(uniqueObjectsBezier1), len(uniqueObjectsBezier2)))
            resultlog.info('  Mismatching Objects Sum: {0}'.format(len(ObjectsBezier1)))
            return 0
        else:
            resultlog.info('    OK')
            return 1

#check which file format is given for further calculations
def GetObjects(fileName1, fileName2) :
    print("FILE1")
    if fileName1[-3] == 'd' and fileName1[-2] == 'x' and fileName1[-1] == 'f' :
        #file1 = DXF(fileName1)
        points1, line1, Bezier1, CP1 = DXF.CalculateObjects(fileName1)
        fileName1 = fileName1.replace('.dxf', '')
    else :
        #file1 = SVG(fileName1)
        points1, line1, Bezier1, CP1 = SVG.CalculateObjects(fileName1)
        fileName1 = fileName1.replace('.svg', '')
    print("FILE1")
    print("FILE2")
    if fileName2[-3] == 'd' and fileName2[-2] == 'x' and fileName2[-1] == 'f' :
        #file2 = DXF.__init__(fileName2)
        points2, line2, Bezier2, CP2 = DXF.CalculateObjects(fileName2)
        fileName2 = fileName1.replace('.dxf', '')
    else :
        #file2 = SVG(fileName2)
        points2, line2, Bezier2, CP2 = SVG.CalculateObjects(fileName2)
        fileName2 = fileName1.replace('.svg', '')
    print("FILE2")

    visited, absMaxDist, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2, uniqueObjectsBezier1, uniqueObjectsBezier2 = Files.Compare(points1, line1, Bezier1, CP1, points2, line2, Bezier2, CP2)
    if showPlot == True :
        Files.PlotFiles(absMaxDist, Bezier1, CP1, Bezier2, CP2, equalObjectsBezier, equalObjectsBezier2, ObjectsBezier1, ObjectsBezier2, uniqueObjectsBezier1, uniqueObjectsBezier2)
    return absMaxDist, uniqueObjectsBezier1, uniqueObjectsBezier2, ObjectsBezier1, fileName1, fileName2

def MyForm(x, p):
    return x/10


f1 = str(Path(filename1).resolve())
f2 = str(Path(filename2).resolve())
resultlog.info('{0}       {1}'.format(f1, f2))

if not os.path.exists(filename1):
    resultlog.error('ERROR:     {0} not found'.format(str(Path('CurveComparison').resolve()) + '\\' + filename1))
if not os.path.exists(filename2):
    resultlog.error('ERROR:     {0} not found'.format(str(Path('CurveComparison').resolve()) + '\\' + filename2))


nr = 1
isOkCnt = 0
mismatchFiles = {}
unfoundcnt = 0

if (filename1.endswith(".dxf") or filename1.endswith(".svg")) and (filename2.endswith(".dxf") or filename2.endswith(".svg")) :
    fileName1 = filename1
    fileName2 = filename2
    if showPlot == True :
        #Graph formatting
        plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
        plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        major_ticks = numpy.arange(-1000, 1000, 100)
        minor_ticks = numpy.arange(-1000, 1000, 10)
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
    isOk = Files.WriteLog(nr, f1, f2, absMaxDist, maxOffset, uniqueObjectsBezier1, uniqueObjectsBezier2, ObjectsBezier1)
    isOkCnt += isOk
    if isOk == 0:
        mismatchFiles[(f1, f2)] = absMaxDist / 10
    nr += 1
    if showPlot == True :
        plt.show()
elif not (filename1.endswith(".dxf") or filename1.endswith(".svg")) or (filename2.endswith(".dxf") or filename2.endswith(".svg")) :
    unfoundPair = []
    for filename in os.listdir(filename1) :
        if filename.endswith(".dxf") or filename.endswith(".svg") :
            fileName1tmp = str(Path(filename1).resolve()) + '\\' + filename
            if filename.endswith('.dxf') :
                filename = filename.replace('.dxf', '')
            if filename.endswith('.svg') :
                filename = filename.replace('.svg', '')
            existsFile2 = False
            for filenametmp in os.listdir(filename2) :
                if (filenametmp.endswith(".dxf") or filenametmp.endswith(".svg")) and (filenametmp.startswith(filename + '_')):
                    if showPlot == True :
                        #Graph formatting
                        plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
                        plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
                        fig = plt.figure()
                        ax = fig.add_subplot(1, 1, 1)
                        major_ticks = numpy.arange(-1000, 1000, 100)
                        minor_ticks = numpy.arange(-1000, 1000, 10)
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

                    fileName2tmp = filename2 + '\\' + filenametmp
                    absMaxDist, uniqueObjectsBezier1, uniqueObjectsBezier2, ObjectsBezier1, fileName1, fileName2 = GetObjects(fileName1tmp, fileName2tmp)
                    isOk = Files.WriteLog(nr, fileName1tmp, fileName2tmp, absMaxDist, maxOffset, uniqueObjectsBezier1, uniqueObjectsBezier2, ObjectsBezier1)
                    isOkCnt += isOk
                    if isOk == 0:
                        mismatchFiles[(fileName1tmp, fileName2tmp)] = absMaxDist / 10
                    nr += 1
                    existsFile2 = True
                    if showPlot == True :
                        plt.show()
            if existsFile2 == False :
                unfoundcnt += 1
                unfoundPair.append(filename)
                resultlog.error('ERROR:     matching file not found for {0}'.format(fileName1tmp))
        else :
            continue
else :
    resultlog.error('ERROR:     given file types don´t match')


if isOkCnt != nr - 1 :
    resultlog.info('')
    resultlog.info('Number of Matching files (max offset: {2} cm):      {0} / {1}'.format(isOkCnt, nr - 1, maxOffset))
    resultlog.info('')
    resultlog.info('files that don`t match:')
    cnt = 0
    for filenames in mismatchFiles :
        cnt += 1
        resultlog.info('{0}. max offset: {1} cm   files: {2}'.format(cnt, mismatchFiles[filenames], filenames))
if unfoundcnt > 0:
    resultlog.info('')
    resultlog.info("{0} files (in {1}) with unfound comparable file: ".format(unfoundcnt, Path(filename1).resolve()))
    cnt = 0
    for filename in unfoundPair :
        cnt += 1
        resultlog.info('    {0}    {1}'.format(cnt, filename))

stat = open(statFile, "r")
textLines = stat.readlines()
linecnt = len(textLines)

if len(textLines) == 0 and nr - 1 != 0:
    stat.close()
    stat = open(statFile, "w")
    print(nr-1)
    stat.write("{0} / {1} compared files match (max offset {2} cm)  {3}% OK".format(isOkCnt, nr - 1, maxOffset, 100 / (nr - 1) * isOkCnt))
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

statInfo = open(statFile.replace('.log', '_INFO.log'), "r")
textLines = statInfo.readlines()
linecnt = len(textLines)
nr = 0
if linecnt > 0 :
    i = 0
    nr = ''
    while textLines[-1][i] != '.' :
        nr = nr + textLines[-1][i]
        i += 1
    nr = int(nr)
statInfo.close()

statInfo = open(statFile.replace('.log', '_INFO.log'), "a")
if len(mismatchFiles) > 0 :
    statInfo.write("-------------------------------------------------\n")
    statInfo.write("{0}. {1}_{2}.log   {3} file pairs offsets exceeds maximum acceptable offset {4}cm\n".format(nr + 1, fileName1, fileName2, len(mismatchFiles), maxOffset))
