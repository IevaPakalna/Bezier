import calculations
import numpy



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
            dist1 = calculations.PointDist(CP[i][0][0], CP[j][0][0])
            dist2 = calculations.PointDist(CP[i][0][0], CP[j][-1][3])
            dist3 = calculations.PointDist(CP[i][-1][3], CP[j][-1][3])

            a1 = calculations.LineSlope(P11, P12)
            a2 = calculations.LineSlope(P21, P22)

            #find the largest distance
            if dist1 >= dist2 and dist1 >= dist3 and dist1 > maxDist1 :
                maxDist1 = dist1
                slope = calculations.LineSlope(CP[i][0][0], CP[j][0][0])
                alpha1 = numpy.degrees(numpy.arctan((slope - a1) / (1 + a1 * slope)))
                alpha2 = numpy.degrees(numpy.arctan((slope - a2) / (1 + a2 * slope)))

                angle = numpy.degrees(numpy.arctan((slope) / (1)))
            if dist2 >= dist1 and dist2 >= dist3 and dist2 > maxDist2:
                maxDist2 = dist2
                slope = calculations.LineSlope(CP[i][0][0], CP[j][-1][3])
                alpha1 = numpy.degrees(numpy.arctan((slope - a1) / (1 + a1 * slope)))
                alpha2 = numpy.degrees(numpy.arctan((slope - a2) / (1 + a2 * slope)))

                angle = numpy.degrees(numpy.arctan((slope) / (1)))
            if dist3 >= dist1 and dist3 >= dist2 and dist3 > maxDist3:
                maxDist = dist3
                slope = calculations.LineSlope(CP[i][-1][3], CP[j][-1][3])
                alpha1 = numpy.degrees(numpy.arctan((slope - a1) / (1 + a1 * slope)))
                alpha2 = numpy.degrees(numpy.arctan((slope - a2) / (1 + a2 * slope)))

                angle = numpy.degrees(numpy.arctan((slope) / (1)))

            angle1 = numpy.degrees(numpy.arctan((a1) / (1)))
            angle2 = numpy.degrees(numpy.arctan((a2) / (1)))

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
#        #a1 = calculations.LineSlope(P11, P12)
#        #a2 = calculations.LineSlope(P21, P22)
#        #alpha = numpy.degrees(numpy.arctan((a2 - a1) / (1 + a1 * a2)))
#        #dist1tmp = calculations.PointDist(P11, P12)
#        #dist2tmp = calculations.PointDist(P21, P22)
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
#            points[i][0] = numpy.cos(alpha) * (points[i][0] - P11[0]) - numpy.sin(alpha) * (points[i][1] - P11[1]) + P11[0]
#            points[i][1] = numpy.sin(alpha) * (points[i][0] - P11[0]) + numpy.cos(alpha) * (points[i][1] - P11[1]) + P11[1]

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
#                    CP[i][j][k][0] = numpy.cos(alpha) * (CP[i][j][k][0] - P11[0]) - numpy.sin(alpha) * (CP[i][j][k][1] - P11[1]) + P11[0]
#                    CP[i][j][k][1] = numpy.sin(alpha) * (CP[i][j][k][0] - P11[0]) + numpy.cos(alpha) * (CP[i][j][k][1] - P11[1]) + P11[1]

                CP[i][j][k] = tuple(CP[i][j][k])
    CP2.clear()
    return points, line2, CP
