import msprime
import numpy as np
import sys
import timeit
import struct

nodes = msprime.NodeTable()
edges = msprime.EdgeTable()

g = []
with open(sys.argv[1],"rb") as f:
    while True:
        a=struct.unpack('i',f.read(4))
        if a[0]==-1:
            break
        t=struct.unpack('d',f.read(8))
        g.append(t[0])
    # for line in f:
    #     l = line.rstrip().split(" ")
    #     g.append(float(l[1]))
time = np.array(g)
time -= time.max()
time *= -1.0

nodes.append_columns(time=time,flags=[-1]*len(time))

p = []
c = []
l = []
r = []

with open(sys.argv[2],"rb") as f:
    while True:
        pi=struct.unpack('i',f.read(4))
        if pi[0] == -1:
            break
        ci=struct.unpack('i',f.read(4))
        li=struct.unpack('d',f.read(8))
        ri=struct.unpack('d',f.read(8))
        p.append(pi[0])
        c.append(ci[0])
        l.append(li[0])
        r.append(ri[0])

edges.set_columns(parent=p,child=c,left=l,right=r)

msprime.sort_tables(nodes=nodes,edges=edges)

N=int(sys.argv[3])
samples=[i for i in range(len(time)-2*N,len(time))] 
n=nodes
e=edges
ts=None
def doit():
    ts=msprime.simplify_tables(nodes=n,edges=e,samples=samples)

time = timeit.timeit(doit,number=1)
#ts = msprime.simplify_tables(nodes=nodes,edges=edges,samples=samples)
print(time)
with open(sys.argv[4],'w') as f:
    for i in edges:
        f.write("{} {} {} {}\n".format(i.parent,i.child,i.left,i.right,nodes[i.parent].time))


