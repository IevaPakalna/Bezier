import dxfgrabber

dxf = dxfgrabber.readfile("svarki_002.dxf")
print("DXF version: {}".format(dxf.dxfversion))
header_var_count = len(dxf.header) # dict of dxf header vars
layer_count = len(dxf.layers) # collection of layer definitions
block_definition_count = len(dxf.blocks) #  dict like collection of block definitions
entity_count = len(dxf.entities) # list like collection of entities


print("header")
print(header_var_count)
print("layer cnt")
print(layer_count)
print("block")
print(block_definition_count)
print("entity")
print(entity_count)

EntitySection = dxf.entities
print("-------")
print("len")
print(EntitySection.__len__())

print("-------")
print("iter")

print(EntitySection.__iter__())

print("-------")
print("get item")
print(EntitySection.__getitem__(0))

print("==================")
#type of objects in file
print("type of objects in file")
type=[entity.dxftype for entity in dxf.entities]
print(type)


output = [entity for entity in dxf.entities]

Color = [entity.color for entity in output]
print(Color)
Linetype = [entity.linetype for entity in output if entity.dxftype != 'POINT']
print(Linetype)

#get parameters of objects
PointCoord = [entity.point for entity in output if entity.dxftype == 'POINT']

LineStart = [entity.start for entity in output if entity.dxftype == 'LINE']
LineEnd = [entity.end for entity in output if entity.dxftype == 'LINE']

CenterPoints = [entity.center for entity in output if entity.dxftype == 'CIRCLE']
Radius = [entity.radius for entity in output if entity.dxftype == 'CIRCLE']

PolylineIsClosed = [entity.is_closed for entity in output if entity.dxftype == 'POLYLINE']
PolylineSplineType = [entity.spline_type for entity in output if entity.dxftype == 'POLYLINE']
PolylinePoints = [entity.points for entity in output if entity.dxftype == 'POLYLINE']
PolylineControlPoints = [entity.control_points for entity in output if entity.dxftype == 'POLYLINE']
PolylineBulge = [entity.bulge for entity in output if entity.dxftype == 'POLYLINE']
PolylineVertexCount = [entity.__len__() for entity in output if entity.dxftype == 'POLYLINE']

LWPolylinePoints = [entity.points for entity in output if entity.dxftype == 'LWPOLYLINE']
LWPolylineIsClosed = [entity.is_closed for entity in output if entity.dxftype == 'LWPOLYLINE']
LWPolylinePointsCount = [entity.__len__() for entity in output if entity.dxftype == 'LWPOLYLINE']

EllipseCenter = [entity.center for entity in output if entity.dxftype == 'ELLIPSE']
EllipseMajorAxis = [entity.major_axis for entity in output if entity.dxftype == 'ELLIPSE']

ArcCenter = [entity.center for entity in output if entity.dxftype == 'ARC']
ArcRadius = [entity.radius for entity in output if entity.dxftype == 'ARC']
ArcStartAngle = [entity.start_angle for entity in output if entity.dxftype == 'ARC']
ArcEndAngle = [entity.end_angle for entity in output if entity.dxftype == 'ARC']

SplineDegree = [entity.degree for entity in output if entity.dxftype == 'SPLINE']
SplineStartTangent = [entity.start_tangent for entity in output if entity.dxftype == 'SPLINE']
SplineEndTangent = [entity.end_tangent for entity in output if entity.dxftype == 'SPLINE']
SplineControlPoints = [entity.control_points for entity in output if entity.dxftype == 'SPLINE']
SplineFitPoints = [entity.fit_points for entity in output if entity.dxftype == 'SPLINE']
SplineKnots = [entity.knots for entity in output if entity.dxftype == 'SPLINE']
SplineIsClosed = [entity.is_closed for entity in output if entity.dxftype == 'SPLINE']


#Point
for x in range(len(PointCoord)):
    print(" POINT\n")

    print(len(PointCoord))
    print(x)
    print("^")
    print(PointCoord[x])
#Line
for x in range(len(LineStart)):

    print(" LINE\n")

    print("Line start")
    print(LineStart[x])

    print("Line end")
    print(LineEnd[x])
#Circle
for x in range(len(CenterPoints)):

    print(" CIRCLE \n")

    print("Center Point")
    print(CenterPoints[x])

    print("Radius")
    print(Radius[x])
#Polyline
for x in range(len(PolylineIsClosed)):

    print(" POLYLINE\n")

    print("is closed")
    print(PolylineIsClosed[x])


    print("Spline type")
    print(PolylineSplineType[x])


    print("line points")
    print(PolylinePoints[x])


    print("Control points")
    print(PolylineControlPoints[x])


    print("bulge")
    print(PolylineBulge[x])


    print("vertex count:")
    print(PolylineVertexCount[x])
#LWPolyline
for x in range(len(LWPolylineIsClosed)):

    print(" LWPOLYLINE\n")

    print(LWPolylinePoints)

    print(LWPolylineIsClosed)

    print(LWPolylinePointsCount)
#Ellipse
for x in range(len(EllipseCenter)):

    print(" ELLIPSE\n")

    print(EllipseCenter)

    print(EllipseMajorAxis)
#Arc
for x in range(len(ArcCenter)):

    print(" ARC\n")

    print(ArcCenter)

    print(ArcRadius)

    print(ArcStartAngle)

    print(ArcEndAngle)
#Spline
for x in range(len(SplineIsClosed)):

    print(" SPLINE\n")
    SplineDegree = [entity.degree for entity in output if entity.dxftype == 'SPLINE']
    print(SplineDegree)
    SplineStartTangent = [entity.start_tangent for entity in output if entity.dxftype == 'SPLINE']
    print(SplineStartTangent)
    SplineEndTangent = [entity.end_tangent for entity in output if entity.dxftype == 'SPLINE']
    print(SplineEndTangent)
    SplineControlPoints = [entity.control_points for entity in output if entity.dxftype == 'SPLINE']
    print(SplineControlPoints)
    SplineFitPoints = [entity.fit_points for entity in output if entity.dxftype == 'SPLINE']
    print(SplineFitPoints)
    SplineKnots = [entity.knots for entity in output if entity.dxftype == 'SPLINE']
    print(SplineKnots)
    SplineIsClosed = [entity.is_closed for entity in output if entity.dxftype == 'SPLINE']
    print(SplineIsClosed)




print("=============")
#Block
#print("BLOCKS")
#BlockBasepoint= [Block.BlockBasepoint for Block in dxf.entities]
#print(BlockBasepoint)
