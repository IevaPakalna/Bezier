import calculations
import calculateBezier
import transformation
from sympy import init_printing, Symbol, UnevaluatedExpr, expand, pretty_print as pprint, latex
import numpy

#arranges lines in a way that line start is closer to (0,0) point than end
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

#reads value of parameter
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
#reads parameter (name)
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
            #based on radius, and half of radical axis, will use Pithagorean theorem
            midx = (P1x + P2x) / 2
            midy = (P1y + P2y) / 2
            dist = calculations.PointDist((P1x, P1y), (midx, midy))

            edge3 = numpy.sqrt(round((rx)**2 - (dist)**2, 5))
            #get perpendicular formula which also will coincide with 3rd edge
            PF = calculations.PerpFormula((P1x, P1y), (midx, midy), (midx, midy))

            x = Symbol('x')
            tmpx = -1
            tmpy = PF.subs(x, tmpx)

            #get two possible circle centers
            c1 = []
            c2 = []
            c1 = calculations.DistantPoint((midx, midy), (midx, midy), (tmpx, tmpy), edge3)
            c2 = calculations.DistantPoint((midx, midy), (midx, midy), (tmpx, tmpy), - edge3)
            c = []
            #get the correct arc based on largeArcFlag and sweepFlag from 4 possible options
            radToDeg = 57.2957795

            a1 = calculations.LineSlope(c1, (P1x, P1y))
            a2 = calculations.LineSlope(c1, (P2x, P2y))
            t1 = numpy.arctan(a1 / 1)
            if P1x < c1[0] :
                t1 = ((t1 + 180 / radToDeg) % (2 * numpy.pi))
            a2 = calculations.LineSlope(c1, (P2x, P2y))
            t2 = numpy.arctan(a2 / 1)
            if P2x < c1[0] :
                t2 = ((t2 + 180 / radToDeg) % (2 * numpy.pi))
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

            a1 = calculations.LineSlope(c, (P1x, P1y))
            a2 = calculations.LineSlope(c, (P2x, P2y))
            t1 = numpy.arctan(a1 / 1)
            if P1x < c[0] :
                t1 = ((t1 + 180 / radToDeg) % (2 * numpy.pi))
            a2 = calculations.LineSlope(c, (P2x, P2y))
            t2 = numpy.arctan(a2 / 1)
            if P2x < c[0] :
                t2 = ((t2 + 180 / radToDeg) % (2 * numpy.pi))
            t = t1

            if sweepFlag == 0 :
                cnt = -2 * numpy.pi / 18
            else :
                cnt = 2 * numpy.pi / 18

            if (t1 < t2 and sweepFlag == 1) or (t1 > t2 and sweepFlag == 0):
                while t < t2 :
                    xprim = rx * numpy.cos(t) + c[0]
                    yprim = rx * numpy.sin(t) + c[1]
                    polyPoints.append((xprim, yprim))
                    t += cnt
                if t - cnt == t1 :
                    t += cnt / 2
                    xprim = rx * numpy.cos(t) + c[0]
                    yprim = rx * numpy.sin(t) + c[1]
                    polyPoints.append((xprim, yprim))
                xprim = rx * numpy.cos(t2) + c[0]
                yprim = rx * numpy.sin(t2) + c[1]
                polyPoints.append((xprim, yprim))
            else :
                t = t2
                while t < 360 / radToDeg :
                    xprim = rx * numpy.cos(t) + c[0]
                    yprim = rx * numpy.sin(t) + c[1]
                    polyPoints.append((xprim, yprim))
                    t += cnt
                xprim = rx * numpy.cos(0) + c[0]
                yprim = rx * numpy.sin(0) + c[1]
                polyPoints.append((xprim, yprim))
                t = 0
                while t < t1 :
                    xprim = rx * numpy.cos(t) + c[0]
                    yprim = rx * numpy.sin(t) + c[1]
                    polyPoints.append((xprim, yprim))
                    t += cnt
                xprim = rx * numpy.cos(t1) + c[0]
                yprim = rx * numpy.sin(t1) + c[1]
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
                cx, pos = GetValue(text, pos, length)
                continue
            if attrib == 'cy' :
                cy, pos = GetValue(text, pos, length)
                continue
            if attrib == 'r' :
                r, pos = GetValue(text, pos, length)
                continue
            if attrib == 'pathLength' :

                continue
            if attrib == 'stroke' :
                color, pos = GetStrValue(text, pos, length)
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
            x = r * numpy.cos(t)
            y = r * numpy.sin(t)
            polylinePoints.append((x, y))
            t += 36

    return points, polylinePoints, isShape

            #these we will add later...or not (seems that we dont need them)
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
                x, pos = GetValue(text, pos, length)
            if attrib == 'y' :
                y, pos = GetValue(text, pos, length)
            if attrib == 'dx' :
                dx, pos = GetValue(text, pos, length)
            if attrib == 'dy' :
                dy, pos = GetValue(text, pos, length)
            if attrib == 'rotate' :
                alpha, pos = GetValue(text, pos, length)
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
                x1, pos = GetValue(text, pos + 1, length)
            if attrib == 'y1' :
                y1, pos = GetValue(text, pos + 1, length)
            if attrib == 'x2' :
                x2, pos = GetValue(text, pos + 1, length)
            if attrib == 'y2' :
                y2, pos = GetValue(text, pos + 1, length)
#            if attrib == 'pathLength' :
            if attrib == 'stroke' :
                color, pos = GetStrValue(text, pos + 1, length)
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
                CPoints, polyPoints, pos = GetControlPoints(text, pos + 2, length)
                continue
            if attrib == 'pathLength' :

                continue
            if attrib == 'stroke' :
                color, pos = GetStrValue(text, pos + 1, length)
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
    file2 = open(file2)
    textLines = file2.readlines()

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
                    points, polylinePoints, isShape = WriteCircle(textLine, i + 1, lenText, points, polylinePoints)
                if element == 'text' :
                    pointTags = WriteText(textLine, i + 1, lenText, pointTags, group)
                if element == 'ellipse' :
                    WriteEllipse(textLine, i + 1, lenText)
                if element == 'line' :
                    CPtmp, isShape = WriteLine(textLine, i + 1, lenText)
                    if isShape == False :
                        continue
#                        try :
#                            CP2.index([CPtmp])
#                        except :
                    pos = SortInsertPosLines(CPtmp[0], CP2, 0, linesSum)
                    CP2.insert(pos, [CPtmp])
                    linesSum += 1
                if element == 'path' :
                    CPtmp, polylinePoints, isShape = WritePath(textLine, i + 1, lenText, polylinePoints)
                    if isShape == False :
                        continue
                    if len(CPtmp) != 0 :
                        #for sorting we need to make sure Bezier is formed in the right way, meaning that the segments are corresponding
                        if CPtmp[0][0] > CPtmp[3][0] :
                            CPtmp.reverse()
                        elif CPtmp[0][1] > CPtmp[3][1] :
                            CPtmp.reverse()
                        pos = SortInsertPosLines(CPtmp[0], CP2, 0, linesSum)
                        CP2.append([CPtmp])
                        linesSum += 1
                if element == 'polygon' :
                    WritePolygon(textLine, i + 1, lenText)
                if element == 'polyline' :
                    WritePolyline(textLine, i + 1, lenText)
                if element == 'rect' :
                    WriteRect(textLine, i + 1, lenText)
                break
    #compress lines, so there are not seperate line segments
    lenCP2 = len(CP2)
    j = 0
    while j < len(CP2) :
        #length of CP1 is changing as lines are being compressed
        i = j + 1
        while i < len(CP2) :
            a1 = calculations.PointDist(CP2[i][0][0], CP2[i][0][1])
            a2 = calculations.PointDist(CP2[i][0][1], CP2[i][0][2])
            a3 = calculations.PointDist(CP2[i][0][2], CP2[i][0][3])

            a4 = calculations.PointDist(CP2[j][-1][0], CP2[j][-1][1])
            a5 = calculations.PointDist(CP2[j][-1][1], CP2[j][-1][2])
            a6 = calculations.PointDist(CP2[j][-1][2], CP2[j][-1][3])
            if (((a1 == 0 and a2 == 0) or (a2 == 0 and a3 == 0) or (a1 == 0 and a3 == 0)) and ((a4 == 0 and a5 == 0) or (a5 == 0 and a6 == 0) or (a4 == 0 and a6 == 0))) :
                a1 = calculations.LineSlope(CP2[j][-1][0], CP2[j][-1][3])
                a3 = calculations.LineSlope(CP2[i][0][0], CP2[i][0][3])
                if a1 == a3 :
                    t = Symbol('t')
                    a1Fx, a1Fy = calculateBezier.ParamLineFormula(CP2[j][-1][0], CP2[j][-1][3])
                    a3Fx, a3Fy = calculateBezier.ParamLineFormula(CP2[i][0][0], CP2[i][0][3])
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
                dist1 = calculations.PointDist(CP2[j][-1][0], CP2[j][-1][1])
                dist2 = calculations.PointDist(CP2[j][-1][1], CP2[j][-1][2])
                dist3 = calculations.PointDist(CP2[j][-1][2], CP2[j][-1][3])
                if dist3 != 0 :
                    P11 = CP2[j][-1][2]
                    P12 = CP2[j][-1][3]
                    a1 = round(calculations.LineSlope(CP2[j][-1][2], CP2[j][-1][3]), 3)
                    a1Fx, a1Fy = calculateBezier.ParamLineFormula(CP2[j][-1][2], CP2[j][-1][3])
                elif dist2 != 0 :
                    P11 = CP2[j][-1][1]
                    P12 = CP2[j][-1][3]
                    a1 = round(calculations.LineSlope(CP2[j][-1][1], CP2[j][-1][3]), 3)
                    a1Fx, a1Fy = calculateBezier.ParamLineFormula(CP2[j][-1][1], CP2[j][-1][3])
                else :
                    P11 = CP2[j][-1][0]
                    P12 = CP2[j][-1][3]
                    a1 = round(calculations.LineSlope(CP2[j][-1][0], CP2[j][-1][3]), 3)
                    a1Fx, a1Fy = calculateBezier.ParamLineFormula(CP2[j][-1][0], CP2[j][-1][3])

                dist1 = calculations.PointDist(CP2[j][0][0], CP2[j][0][1])
                dist2 = calculations.PointDist(CP2[j][0][1], CP2[j][0][2])
                dist3 = calculations.PointDist(CP2[j][0][2], CP2[j][0][3])

                if dist1 != 0 :
                    P21 = CP2[j][0][0]
                    P22 = CP2[j][0][1]
                    a2 = round(calculations.LineSlope(CP2[j][0][0], CP2[j][0][1]), 3)
                    a2Fx, a2Fy = calculateBezier.ParamLineFormula(CP2[j][0][0], CP2[j][0][1])
                elif dist2 != 0 :
                    P21 = CP2[j][0][0]
                    P22 = CP2[j][0][2]
                    a2 = round(calculations.LineSlope(CP2[j][0][0], CP2[j][0][2]), 3)
                    a2Fx, a2Fy = calculateBezier.ParamLineFormula(CP2[j][0][0], CP2[j][0][2])
                else :
                    P21 = CP2[j][0][0]
                    P22 = CP2[j][0][3]
                    a2 = round(calculations.LineSlope(CP2[j][0][0], CP2[j][0][3]), 3)
                    a2Fx, a2Fy = calculateBezier.ParamLineFormula(CP2[j][0][0], CP2[j][0][3])

                dist1 = calculations.PointDist(CP2[i][-1][0], CP2[i][-1][1])
                dist2 = calculations.PointDist(CP2[i][-1][1], CP2[i][-1][2])
                dist3 = calculations.PointDist(CP2[i][-1][2], CP2[i][-1][3])
                if dist3 != 0 :
                    P31 = CP2[i][-1][2]
                    P32 = CP2[i][-1][3]
                    a3 = round(calculations.LineSlope(CP2[i][-1][2], CP2[i][-1][3]), 3)
                    a3Fx, a3Fy = calculateBezier.ParamLineFormula(CP2[i][-1][3], CP2[i][-1][2])
                elif dist2 != 0 :
                    P31 = CP2[i][-1][1]
                    P32 = CP2[i][-1][3]
                    a3 = round(calculations.LineSlope(CP2[i][-1][1], CP2[i][-1][3]), 3)
                    a3Fx, a3Fy = calculateBezier.ParamLineFormula(CP2[i][-1][3], CP2[i][-1][1])
                else :
                    P31 = CP2[i][-1][0]
                    P32 = CP2[i][-1][3]
                    a3 = round(calculations.LineSlope(CP2[i][-1][0], CP2[i][-1][3]), 3)
                    a3Fx, a3Fy = calculateBezier.ParamLineFormula(CP2[i][-1][3], CP2[i][-1][0])

                dist1 = calculations.PointDist(CP2[i][0][0], CP2[i][0][1])
                dist2 = calculations.PointDist(CP2[i][0][1], CP2[i][0][2])
                dist3 = calculations.PointDist(CP2[i][0][2], CP2[i][0][3])
                if dist1 != 0 :
                    P41 = CP2[i][0][0]
                    P42 = CP2[i][0][1]
                    a4 = round(calculations.LineSlope(CP2[i][0][0], CP2[i][0][1]), 3)
                    a4Fx, a4Fy = calculateBezier.ParamLineFormula(CP2[i][0][0], CP2[i][0][1])
                elif dist2 != 0 :
                    P41 = CP2[i][0][0]
                    P42 = CP2[i][0][2]
                    a4 = round(calculations.LineSlope(CP2[i][0][0], CP2[i][0][2]), 3)
                    a4Fx, a4Fy = calculateBezier.ParamLineFormula(CP2[i][0][0], CP2[i][0][2])
                else :
                    P41 = CP2[i][0][0]
                    P42 = CP2[i][0][3]
                    a4 = round(calculations.LineSlope(CP2[i][0][0], CP2[i][0][3]), 3)
                    a4Fx, a4Fy = calculateBezier.ParamLineFormula(CP2[i][0][0], CP2[i][0][3])

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
            Beziertmp, CPtmp = Lines(polylinePoints[0], polylinePoints[1], Bezier, CP2, linesSum)
            continue
        Beziertmp, CPtmp = calculateBezier.CompositeBezier(i, 2)
        CP2.append(CPtmp)
    points, line, CP2 = transformation.Transformation(points, line, CP2)

    Beziertmp = [[]]
    for i in CP2 :
        for j in i:
            Bx, By = calculateBezier.BezierFormulaComp(j)
            Beziertmp[0].append([Bx, By])
        Bezier.extend(Beziertmp)

        Beziertmp.clear()
        Beziertmp.append([])

#        for i in polylinePoints :
#            Beziertmp, CPtmp = calculateBezier.CompositeBezier(i, 2)
#            CP2.append(CPtmp)
#            Bezier.append(Beziertmp)
    return points, line, Bezier, CP2
