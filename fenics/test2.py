from fenics import *
from CSVtoNodalValues import *
import math
xmin = -0.00175
ymin = -0.00175
zmin = -0.0002
stencil   = 0.00005
xmax = 0.00175
ymax = 0.00175
zmax = 0.0

mesh = BoxMesh(Point(xmin,ymin,zmin),Point(xmax,ymax,zmax),1,1,1)
#
# # print "dof to vertex map - ME"
V = FunctionSpace(mesh,'CG',1)
#
for v in vertices(mesh):
    print(v.point().x(),v.point().y(),v.point().z())
#
print(vertex_to_dof_map(V))

# f = Function(V)
# f.vector()[1] =1
