import msprime
import numpy as np

nodes = msprime.NodeTable()
edges = msprime.EdgeTable()

g = []
with open("test_nodes.txt","r") as f:
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

with open("test_edges.txt","r") as f:
    for line in f:
        li = line.rstrip().split(" ")
        p.append(int(li[0]))
        c.append(int(li[1]))
        l.append(float(li[2]))
        r.append(float(li[3]))

edges.set_columns(parent=p,child=c,left=l,right=r)

msprime.sort_tables(nodes=nodes,edges=edges)

for i in edges:
    print(i.parent,i.child,i.left,i.right,nodes[i.parent].time)
