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
mesh = BoxMesh(Point(xmin,ymin,zmin),Point(xmax,ymax,zmax),70,70,4)


E0 = 10e5
nu = 0.28
alpha = 5e-6
Tref = 20
E  = E0
mu = E/(2*(1+nu))
lamb = E*nu/((1+nu)*(1-2*nu))

V = VectorFunctionSpace(mesh,'CG',1)
Vs = FunctionSpace(mesh,'CG',1)
u = TrialFunction(V)
d = u.geometric_dimension()
v = TestFunction(V)
data = getTimestep(0,99)


# Assign the Temperature term manually
delta_T = Function(Vs)
vtodmap = vertex_to_dof_map(Vs)
for vert in vertices(mesh):
    idx = (vert.index())
    hashTup = (vert.point().x(),vert.point().y(),vert.point().z())
    delta_T.vector()[vtodmap[idx]]   = data[hashTup].T - Tref


# Boundary Conditions
def centerPinned(x,on_boundary):
    isOnBound = ( math.sqrt(x[0] ** 2 + x[1] ** 2) < 10.0 * stencil)
    if isOnBound:
        if abs(x[2] - zmin) > stencil * 0.1:
            isOnBound = False
    return isOnBound

def clamped_boundary(x,on_boundary):
    foo = (x[0] > xmax-stencil*0.5) or (x[0] < xmin+stencil*0.5)
    bar = (x[1] > ymax-stencil*0.5) or (x[1] < ymin+stencil*0.5)
    return on_boundary and (foo or bar)

bc = DirichletBC(V,Constant((0,0,0)),clamped_boundary)
# define strain-displacement
def epsilon(u):
    return 0.5*(nabla_grad(u)+nabla_grad(u).T)
# define stifness tensor
# def sigma(u):
    # return lamb*((nabla_div(u)-alpha*delta_T*(3+2*mu/lamb))*Identity(d))+2*mu*epsilon(u)
def sigma(u):
    return lamb*nabla_div(u)*Identity(d)+2*mu*epsilon(u)

# Assign the source term manually
f = Function(V)
vtodmap = vertex_to_dof_map(V)
for vert in vertices(mesh):
    idx = 3 * (vert.index())
    hashTup = (vert.point().x(),vert.point().y(),vert.point().z())
    f.vector()[vtodmap[idx]]   = data[hashTup].f[0]
    f.vector()[vtodmap[idx+1]] = data[hashTup].f[1]
    f.vector()[vtodmap[idx+2]] = data[hashTup].f[2]

# fff = interpolate(f, V)


form = inner(sigma(u),epsilon(v))*dx - dot(f, v)*dx
a = lhs(form)
L = rhs(form)
# a = inner(sigma(u),epsilon(v))*dx
# L = dot(f, v)*dx

print("-------------SOLVING-------------")
u = Function(V)
solve(a == L,u,bc)

vtkfile = File('./thermoelasticity.pvd')
vtkfile << u

# vvv = VectorFunctionSpace(mesh, "CG", 1)
#
vtkfile = File('./thermoelasticity_f.pvd')
vtkfile << f

vtkfile = File('./thermoelasticity_T.pvd')
vtkfile << delta_T
