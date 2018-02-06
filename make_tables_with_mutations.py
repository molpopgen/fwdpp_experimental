import msprime
import struct
import numpy as np

ts = msprime.simulate(20, mutation_rate=10., recombination_rate=20.,random_seed=42)

print(len(ts.tables.nodes),len(ts.tables.edges),len(ts.tables.sites),len(ts.tables.mutations))
print(ts.tables.edges)
print(ts.tables.sites)
print(ts.tables.mutations)

with open('test_nodes.bin', 'wb') as f:
    s = struct.pack('i', len(ts.tables.nodes))
    f.write(s)
    for n in ts.tables.nodes:
        s = struct.pack('=d', n.time)
        f.write(s)

with open('test_edges.bin', 'wb') as f:
    s = struct.pack('i', len(ts.tables.edges))
    f.write(s)
    for e in ts.tables.edges:
        s = struct.pack('=d', e.left)
        f.write(s)
        s = struct.pack('=d', e.right)
        f.write(s)
        s = struct.pack('i', e.parent)
        f.write(s)
        s = struct.pack('i', e.child)
        f.write(s)

with open('test_sites.bin', 'wb') as f:
    s = struct.pack('i', len(ts.tables.sites))
    f.write(s)
    for site, mut in zip(ts.tables.sites,ts.tables.mutations):
        s = struct.pack('i', mut.node)
        f.write(s)
        s = struct.pack('=d', site.position)
        f.write(s)

with open('test_variant_counts.txt','w') as f:
    for i in ts.variants():
        pos=i.position
        n=len(np.where(i.genotypes==1)[0])
        print(pos,n,i.genotypes)
        f.write("{} {}\n".format(pos,len))
