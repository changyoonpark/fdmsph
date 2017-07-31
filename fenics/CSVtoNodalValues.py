import math

class Mesh:

    def hashFunc(self,a):
        return (format(a[0],'.10f'),format(a[1],'.10f'),format(a[2],'.10f'))

    def norm2(self,pos1,pos2):
        foo = (pos1[0]-pos2[0])**2 + (pos1[1]-pos2[1])**2 + (pos1[2]-pos2[2])**2
        return math.sqrt(foo)

    def __init__(self):
        self.nodes = dict()

    def __setitem__(self,idx,item):
        self.nodes[self.hashFunc(idx)] = item

    def __getitem__(self,idx):
        return self.nodes[self.hashFunc(idx)]

    # def findNode(self,pos,tol):
    #     for node in self.nodes:
    #         if self.norm2(node.pos,pos) < tol:
    #             return node
    #     return None


class Node:
    def __init__(self,pos,T,f):
        self.pos = pos
        self.T = T
        self.f = f
    def __str__(self):
        string =  "pos : {}, {}, {}\n".format(self.pos[0],self.pos[1],self.pos[2])
        string += "T : {}\n".format(self.T)
        string += "f : {}, {}, {}".format(self.f[0],self.f[1],self.f[2])
        return string

def getTimestep(i,j):

    f = open('../outputs/out{}/out.{}.csv'.format(i,j), 'r')

    dl = dict()
    n = 0
    mesh = Mesh()

    for line in f:
        lineSplit = line.split(",")
        if n == 0:
            l = 0
            for entry in lineSplit:
                dl[entry.replace('"','')] = l
                l += 1
        else :
            a = [float(lineSplit[dl["x"]]),float(lineSplit[dl["y"]]),float(lineSplit[dl["z"]])]
            b = float(lineSplit[dl["temp"]])
            c = [float(lineSplit[dl["fxSensed"]]),float(lineSplit[dl["fySensed"]]),float(lineSplit[dl["fzSensed"]])]
            node = Node(a,b,c)
            mesh[a] = node

        n += 1

    return mesh
