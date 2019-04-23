import dxfgrabber

dxf = dxfgrabber.readfile("kvadrats_002.dxf")
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

for x in type:
#Circle
    if x == 'CIRCLE':
        print("CIRCLE")
        CenterPoints = [entity.center for entity in dxf.entities]
        print("Center Point")
        print(CenterPoints)
        Radius = [entity.radius for entity in dxf.entities]
        print("Radius")
        print(Radius)
#Polyline
    if x == 'POLYLINE':
        print("POLYLINE")
        PolylineIsClosed = [Polyline.is_closed for Polyline in dxf.entities]
        print("is closed")
        print(PolylineIsClosed)

        PolylineSplineType = [Polyline.spline_type for Polyline in dxf.entities]
        print("Spline type")
        print(PolylineSplineType)

        PolylinePoints = [Polyline.points for Polyline in dxf.entities]
        print("line points")
        print(PolylinePoints)

        PolylineControlPoints = [Polyline.control_points for Polyline in dxf.entities]
        print("Control points")
        print(PolylineControlPoints)

        PolylineBulge = [Polyline.bulge for Polyline in dxf.entities]
        print("bulge")
        print(PolylineBulge)

        PolylineVertexCount = [Polyline.__len__() for Polyline in dxf.entities]
        print("vertex count:")
        print(PolylineVertexCount)
#Point
    if x == 'POINT':
        print("POINT")
        PointCoord = [entity.point for entity in dxf.entities]
        print(PointCoord)
#Line
    if x == 'LINE':
        print("LINE")
        LineStart = [Line.start for Line in dxf.blocks]
        print("Line start")
        print(LineStart)
        LineEnd = [Line.end for Line in dxf.blocks]
        print("Line end")
        print(LineEnd)


CenterPoints = [entity.center for entity in dxf.entities]




print("=============")
#Block
print("BLOCKS")
BlockBasepoint= [block.BlockBasepoint for block in dxf.blocks]
print(BlockBasepoint)
