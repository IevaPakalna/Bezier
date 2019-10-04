from sympy import init_printing, Symbol, UnevaluatedExpr, expand, pretty_print as pprint, latex
import numpy

import calculations

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
                mintmp = calculations.PointDist(P, (arr[i][0][0][0], arr[i][0][0][1]))
                if mintmp > minDist :
                    min = mintmp
                    minPind = i
                    lineSt = 0
            elif visited[(arr[i][0][0][0], arr[i][0][0][1])] == 1 :
                mintmp = calculations.PointDist(P, (arr[i][0][0][0], arr[i][0][0][1]))
                if mintmp > minDist :
                    min = mintmp
                    minPind = i
                    lineSt = 0
        if not (arr[i][0][0][0], arr[i][0][0][1]) in visited :
            tmpD = calculations.PointDist(P, (arr[i][0][0][0], arr[i][0][0][1]))
            if tmpD < min and tmpD > minDist:
                min = tmpD
                minPind = i
                lineSt = 0
        elif visited[(arr[i][0][0][0], arr[i][0][0][1])] == 1 :
            tmpD = calculations.PointDist(P, (arr[i][0][0][0], arr[i][0][0][1]))
            if tmpD < min and tmpD > minDist:
                min = tmpD
                minPind = i
                lineSt = 0
        if not (arr[i][-1][3][0], arr[i][-1][3][1]) in visited :
            tmpD = calculations.PointDist(P, (arr[i][-1][3][0], arr[i][-1][3][1]))
            if tmpD < min and tmpD > minDist:
                min = tmpD
                minPind = i
                lineSt = 1
        elif visited[(arr[i][-1][3][0], arr[i][-1][3][1])] == 1 :
            tmpD = calculations.PointDist(P, (arr[i][-1][3][0], arr[i][-1][3][1]))
            if tmpD < min and tmpD > minDist:
                min = tmpD
                minPind = i
                lineSt = 1
    return (arr[minPind][-lineSt][lineSt * 3][0], arr[minPind][-lineSt][lineSt * 3][1]), min

#returns length of cubic Bezier curve
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
        ft = numpy.sqrt(fttmp)
        L += 0.5 * (w[i] * ft)
    return L

#returns point coordinates, given composite Bezier curve, Bezier curves startpoint and distance
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
    #iterate trough every segment till searched point is found
    while (nr >= 0) and (nr <= lenB - 1) :
        #check that parameter is in its interval
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
            #get new point coordinates
            x2tmp = B[nr][0].subs(t, param)
            y2tmp = B[nr][1].subs(t, param)
            Pdist = calculations.PointDist((xtmp, ytmp), (x2tmp, y2tmp))
            disttmp += Pdist
            #if parameter value is equal to segments second endpoint
            if param == abs(round((lineSt - 1), 3)) :
                x2tmp1 = C[nr][lineSt * 3][0]
                y2tmp1 = C[nr][lineSt * 3][1]
                Pdist1 = calculations.PointDist((xtmp, ytmp), (x2tmp1, y2tmp1))
                disttmp1 = disttmp + Pdist1 - Pdist
                #if distance is close enough, then return point coordinates
                if (disttmp1 <= dist + 0.05 and disttmp1 >= dist - 0.05) :
                    return x2tmp1, y2tmp1, nr, disttmp1, lineSt
                if disttmp1 < dist and disttmp1 > disttmp :
                    break

            if (disttmp < dist and param == abs(round((lineSt - 1), 3))):
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
            dist1 += CubicBezierLen(C[nr])
            disttmp = dist1
            nr += 1
            param = lineSt
            m = 0.5
            cnt = 1
        else:
            dist1 += CubicBezierLen(C[nr])
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
        dist = calculations.PointDist(P1, P2)
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
    maxDistParam = MaxBezierTDist(B1, B2, lineSt1, lineSt2)
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
    dist1 = calculations.PointDist((x2, y2), (CP1[0][0][0], CP1[0][0][1]))
    dist2 = calculations.PointDist((x2, y2), (CP1[-1][3][0], CP1[-1][3][1]))
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
    dist1 = calculations.PointDist((x2, y2), (CP1[0][0][0], CP1[0][0][1]))
    dist2 = calculations.PointDist((x2, y2), (CP1[-1][3][0], CP1[-1][3][1]))
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
                dist, P = BezierMinDist(B1[nr][0], B1[nr][1], (x2, y2))
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
    dist1 = calculations.PointDist((x1, y1), P)
    dist2 = calculations.PointDist((x2, y2), P)

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
            dist1 = calculations.PointDist((x1, y1), P)
        else :
            dist1 = 10000
        if param >= 0 + cnt :
            x2 = Bx.subs(t, param - cnt)
            y2 = By.subs(t, param - cnt)
            dist2 = calculations.PointDist((x2, y2), P)
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
        length += CubicBezierLen(CP[i])
    return length
