import msprime
import numpy as np
import sys
import timeit

nodes = msprime.NodeTable()
edges = msprime.EdgeTable()

g = []
with open(sys.argv[1],"r") as f:
    for line in f:
        l = line.rstrip().split(" ")
        g.append(float(l[1]))

time = np.array(g)
time -= time.max()
time *= -1.0

nodes.append_columns(time=time,flags=[-1]*len(time))

p = []
c = []
l = []
r = []

with open(sys.argv[2],"r") as f:
    for line in f:
        li = line.rstrip().split(" ")
        p.append(int(li[0]))
        c.append(int(li[1]))
        l.append(float(li[2]))
        r.append(float(li[3]))

edges.set_columns(parent=p,child=c,left=l,right=r)

msprime.sort_tables(nodes=nodes,edges=edges)

N=int(sys.argv[3])
samples=[i for i in range(len(time)-2*N,len(time))] 
n=nodes
e=edges
def doit():
    msprime.simplify_tables(nodes=n,edges=e,samples=samples)

time = timeit.timeit(doit,number=1)
#ts = msprime.simplify_tables(nodes=nodes,edges=edges,samples=samples)
print(time)
#for i in edges:
#    print(i.parent,i.child,i.left,i.right,nodes[i.parent].time)


