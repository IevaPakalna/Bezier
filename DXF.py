import dxfgrabber
from sympy import init_printing, Symbol, UnevaluatedExpr, expand, pretty_print as pprint, latex

import calculations, calculateBezier

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
        return SortInsertPosLines(P, CP, l, med)
    if CP[med][0][0][0] < P[0] :
        return SortInsertPosLines(P, CP, med + 1, r)
    if CP[med][0][0][0] == P[0] :
        if CP[med][0][0][1] == P[1] :
            return -1 #there already exists identical point, therefore we will not save it
        if CP[med][0][0][1] > P[1] :
            return SortInsertPosLines(P, CP, l, med)
        if CP[med][0][0][1] < P[1] :
            return SortInsertPosLines(P, CP, med + 1, r)
#arranges lines in way that lines start is always closer to (0,0) point than end
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
    pos = SortInsertPosLines(lineStart1, CP1, 0, linesSum)
    if pos == -1 :
        pos = linesSum

    Bezier1.insert(pos, [calculateBezier.BezierFormulaComp([lineStart1, lineStart1, lineEnd1, lineEnd1])])
    CP1.insert(pos, [[lineStart1, lineStart1, lineEnd1, lineEnd1]])

    linesSum += 1

    x = [lineStart1[0], lineEnd1[0]]
    y = [lineStart1[1], lineEnd1[1]]
    return Bezier1, CP1, linesSum

#read dxf file
def CalculateObjects(file1) :
    file = dxfgrabber.readfile(file1)
    type1 = [entity.dxftype for entity in file.entities]
    output1 = [entity for entity in file.entities]
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
            pos = SortInsertPos(point1, points1, 0, r)
            if pos != -1 :
                points1.insert(pos, [point1[0], point1[1]])
                x = [point1[0]]
                y = [point1[1]]

        #Line
        if entity.dxftype == 'LINE':
            P1 = entity.start
            P2 = entity.end
            Bezier1, CP1, linesSum = Lines(P1, P2, Bezier1, CP1, linesSum)
#        #Circle
#        if entity.dxftype == 'CIRCLE':
#            centerPoints1 = entity.center
#            radius1 = entity.radius
#            lenCirclePnts = len(centerPoints1)
#            for i in lenCirclePnts :
#                print("==========================++++==++++=+++++==========")
#                print(centerPoints1[i], radius1[i])
#                circle1.append(centerPoints1[i], radius1[i])
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
                Bezier1, CP1, linesSum = Lines(PolylinePoints1[0], PolylinePoints1[1], Bezier1, CP1, linesSum)
                continue
            if len(PolylinePoints1) == 3 :
                PolylinePoints1.append(PolylinePoints1[-1])
            Beziertmp, CP = calculateBezier.CompositeBezier(PolylinePoints1, 1)
            Bezier1.append(Beziertmp)
            CP1.append(CP)
        #LWPolyline
        if entity.dxftype == 'LWPOLYLINE':
            LWPolylinePoints1 = entity.points
            Beziertmp, CP = calculateBezier.CompositeBezier(LWPolylinePoints1, 1)
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
                    a1 = calculations.LineSlope(CP1[j][-1][0], CP1[j][-1][3])
                    a3 = calculations.LineSlope(CP1[i][0][0], CP1[i][0][3])
                    if a1 == a3 :
                        t = Symbol('t')
                        a1Fx, a1Fy = calculateBezier.ParamLineFormula(CP1[j][-1][0], CP1[j][-1][3])
                        a3Fx, a3Fy = calculateBezier.ParamLineFormula(CP1[i][0][0], CP1[i][0][3])
                        y1 = a1Fx.subs(t, 0)
                        y2 = a1Fx.subs(t, 0)
                        if y1 == y2:
                            if CP1[i][0][0][0] >= CP1[j][-1][0][0] and CP1[i][0][3][0] >= CP1[j][-1][3][0] and CP1[i][0][0][0] <= CP1[j][-1][3][0] and CP1[j][0][0][1] == CP1[i][0][3][1]:
                                if CP1[i][0][0][0] == CP1[j][-1][0][0] and CP1[i][0][3][0] == CP1[j][-1][3][0] and CP1[i][0][0][1] == CP1[j][-1][0][1] and CP1[i][0][3][1] == CP1[j][-1][3][1]:
                                    i += 1
                                    continue
                                CP1[j][-1][3] = CP1[i][0][3]
                                CP1[j][-1][2] = CP1[i][0][3]
                                Bezier1[j][-1] = calculateBezier.BezierFormulaComp(CP1[j][-1])
                                Bezier1.pop(i)
                                CP1.pop(i)
                                i -= 1
                            elif CP1[i][0][0][1] >= CP1[j][-1][0][1] and CP1[i][0][3][1] >= CP1[j][-1][3][1] and CP1[i][0][0][1] <= CP1[j][-1][3][1] and CP1[j][0][0][0] == CP1[i][0][3][0]:
                                if CP1[i][0][0][0] == CP1[j][-1][0][0] and CP1[i][0][3][0] == CP1[j][-1][3][0] and CP1[i][0][0][1] == CP1[j][-1][0][1] and CP1[i][0][3][1] == CP1[j][-1][3][1]:
                                    i += 1
                                    continue
                                CP1[j][-1][3] = CP1[i][0][3]
                                CP1[j][-1][2] = CP1[i][0][3]
                                Bezier1[j][-1] = calculateBezier.BezierFormulaComp(CP1[j][-1])
                                Bezier1.pop(i)
                                CP1.pop(i)
                                i -= 1
                i += 1

        j += 1
    return points1, line1, Bezier1, CP1
