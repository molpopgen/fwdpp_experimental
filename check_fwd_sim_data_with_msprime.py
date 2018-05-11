import msprime
import numpy as np
import sys
import timeit
import time
import struct

print(msprime.__file__,msprime.__version__)

nodes = msprime.NodeTable()
edges = msprime.EdgeTable()

g = []
with open(sys.argv[1],"rb") as f:
    print(sys.argv[1])
    while True:
        a=struct.unpack('i',f.read(4))
        if a[0]==-1:
            break
        t=struct.unpack('d',f.read(8))
        g.append(t[0])
    # for line in f:
    #     l = line.rstrip().split(" ")
    #     g.append(float(l[1]))
times = np.array(g)
times -= times.max()
times *= -1.0

nodes.append_columns(time=times,flags=[-1]*len(times))

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


N=int(sys.argv[3])
samples=[i for i in range(len(times)-2*N,len(times))] 
ts=None

A=time.time()
msprime.sort_tables(nodes=nodes,edges=edges)
B=time.time()
ts=msprime.simplify_tables(nodes=nodes,edges=edges,samples=samples)
C=time.time()

print("Sorting: ",B-A,"seconds")
print("Simplifying: ",C-B,"seconds")


with open(sys.argv[4],'w') as f:
    for i in edges:
        f.write("{} {} {} {}\n".format(i.parent,i.child,i.left,i.right,nodes[i.parent].time))
with open(sys.argv[5],'w') as f:
    for i in nodes:
        f.write("{}\n".format(i.time))


