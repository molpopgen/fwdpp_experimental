import msprime

n=msprime.NodeTable()
e=msprime.EdgeTable()
s=msprime.SiteTable()
m=msprime.MutationTable()

e.add_row(parent=0,child=1,left=0,right=0.4)
e.add_row(parent=0,child=1,left=0.6,right=0.8)
e.add_row(parent=0,child=2,left=0,right=0.8)
e.add_row(parent=3,child=0,left=0.,right=0.4)
e.add_row(parent=3,child=0,left=0.6,right=0.8)

n.add_row(time=1,flags=msprime.NODE_IS_SAMPLE)
n.add_row(time=0,flags=msprime.NODE_IS_SAMPLE)
n.add_row(time=0,flags=msprime.NODE_IS_SAMPLE)
n.add_row(time=2,flags=msprime.NODE_IS_SAMPLE)

s.add_row(position=0.4,ancestral_state='0')
s.add_row(position=0.8,ancestral_state='0')
m.add_row(site=0,node=0,derived_state='1')
m.add_row(site=1,node=3,derived_state='1')

msprime.sort_tables(nodes=n,edges=e,sites=s,mutations=m)

print(e)

idmap = msprime.simplify_tables(nodes=n,edges=e,sites=s,mutations=m,samples=[1,2], sequence_length = 1.0)

print(idmap)
print(e)
print(s)
print(m)

ts=msprime.load_tables(nodes=n,edges=e,sites=s,mutations=m)
print(ts.genotype_matrix())
