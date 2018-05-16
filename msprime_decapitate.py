import msprime
import numpy as np
import struct
import sys

n = int(sys.argv[1])
nedges = int(sys.argv[2])

ts = msprime.simulate(n,mutation_rate=2500,recombination_rate=2500)
t = ts.dump_tables()

print(len(t.nodes),len(t.edges),len(t.mutations))
t.edges.set_columns(parent=t.edges.parent[:nedges],
    child=t.edges.child[:nedges],
    left=t.edges.left[:nedges],
    right=t.edges.right[:nedges])

print(len(t.nodes),len(t.edges),len(t.mutations))
with open("decap_nodes.bin","wb") as f:
    f.write(struct.pack('i',len(t.nodes)))
    for no in t.nodes:
        f.write(struct.pack('d',no.time))
        
with open("decap_edges.bin","wb") as f:
    f.write(struct.pack('i',len(t.edges)))
    for e in t.edges:
        f.write(struct.pack('i',e.parent))
        f.write(struct.pack('i',e.child))
        f.write(struct.pack('d',e.left))
        f.write(struct.pack('d',e.right))

with open("decap_mutations.bin","wb") as f:
    f.write(struct.pack('i',len(t.mutations)))
    for i,j in zip(t.mutations, t.sites):
        f.write(struct.pack('i',i.node))
        f.write(struct.pack('d',j.position))
        
# finally, let's simplify things down so that we can look at
# what msprime thinks

msprime.simplify_tables(nodes=t.nodes,edges=t.edges,sites=t.sites,mutations=t.mutations, samples=np.where(t.nodes.time == 0)[0].astype(np.int32))
print(len(t.nodes),len(t.edges),len(t.mutations))
