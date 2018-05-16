import unittest
import msprime

class Tests(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.n=msprime.NodeTable()
        self.e=msprime.EdgeTable()
        self.s=msprime.SiteTable()
        self.m=msprime.MutationTable()
        self.L=1.0  
    def test1(self):

        self.n.set_columns(time=[1,0,0],flags=[msprime.NODE_IS_SAMPLE]*3)

        self.e.add_row(parent=0,child=1,left=0,right=0.4)
        self.e.add_row(parent=0,child=1,left=0.6,right=1.0)
        self.e.add_row(parent=0,child=2,left=0,right=1)

        self.s.add_row(position=0.5,ancestral_state='0')
        self.m.add_row(site=0,node=1,derived_state='1')

        msprime.sort_tables(nodes=self.n,edges=self.e,
                sites=self.s,mutations=self.m)
        idmap = msprime.simplify_tables(nodes=self.n,edges=self.e,
                sites=self.s,mutations=self.m,samples=[1,2])
        ts = msprime.load_tables(nodes=self.n,edges=self.e,sites=self.s,
                mutations=self.m)
        m = ts.genotype_matrix()
        self.assertEqual(m[0:].sum(),1)

    def test2(self):
        self.n.set_columns(time=[1,0,0,2],flags=[msprime.NODE_IS_SAMPLE]*4)
        
        self.e.add_row(parent=0,child=1,left=0,right=0.4)
        self.e.add_row(parent=0,child=1,left=0.6,right=1.0)
        self.e.add_row(parent=0,child=2,left=0,right=1)
        self.e.add_row(parent=3,child=0,left=0,right=0.4    )

        self.s.add_row(position=0.5,ancestral_state='0')
        self.m.add_row(site=0,node=3,derived_state='1')

        msprime.sort_tables(nodes=self.n,edges=self.e,
                sites=self.s,mutations=self.m)
        idmap = msprime.simplify_tables(nodes=self.n,edges=self.e,
                sites=self.s,mutations=self.m,samples=[1,2])
        ts = msprime.load_tables(nodes=self.n,edges=self.e,sites=self.s,
                mutations=self.m)
        m = ts.genotype_matrix()
        self.assertEqual(m[0:].sum(),0)

    def test3(self):

        self.n.set_columns(time=[1,0,0],flags=[msprime.NODE_IS_SAMPLE]*3)

        self.e.add_row(parent=0,child=1,left=0,right=0.4)
        self.e.add_row(parent=0,child=1,left=0.6,right=1.0)
        self.e.add_row(parent=0,child=2,left=0,right=1)

        self.s.add_row(position=0.4,ancestral_state='0')
        self.m.add_row(site=0,node=1,derived_state='1')

        msprime.sort_tables(nodes=self.n,edges=self.e,
                sites=self.s,mutations=self.m)
        idmap = msprime.simplify_tables(nodes=self.n,edges=self.e,
                sites=self.s,mutations=self.m,samples=[1,2])
        ts = msprime.load_tables(nodes=self.n,edges=self.e,sites=self.s,
                mutations=self.m)
        m = ts.genotype_matrix()
        self.assertEqual(m[0:].sum(),1)

    def test4(self):
        self.n.set_columns(time=[1,0,0,2],flags=[msprime.NODE_IS_SAMPLE]*4)
        
        self.e.add_row(parent=0,child=1,left=0,right=0.4)
        self.e.add_row(parent=0,child=1,left=0.6,right=1.0)
        self.e.add_row(parent=0,child=2,left=0,right=1)
        self.e.add_row(parent=3,child=0,left=0,right=0.4)

        self.s.add_row(position=0.4,ancestral_state='0')
        self.m.add_row(site=0,node=3,derived_state='1')

        msprime.sort_tables(nodes=self.n,edges=self.e,
                sites=self.s,mutations=self.m)
        idmap = msprime.simplify_tables(nodes=self.n,edges=self.e,
                sites=self.s,mutations=self.m,samples=[1,2])
        ts = msprime.load_tables(nodes=self.n,edges=self.e,sites=self.s,
                mutations=self.m)
        m = ts.genotype_matrix()
        self.assertEqual(m[0:].sum(),0)

    def test5(self):
        self.n.set_columns(time=[1,0,0,2],flags=[msprime.NODE_IS_SAMPLE]*4)
        
        self.e.add_row(parent=0,child=1,left=0,right=0.4)
        self.e.add_row(parent=0,child=1,left=0.6,right=1.0)
        self.e.add_row(parent=0,child=2,left=0,right=1)
        self.e.add_row(parent=3,child=0,left=0,right=0.4)
        self.e.add_row(parent=3,child=0,left=0.7,right=0.8)

        self.s.add_row(position=0.7,ancestral_state='0')
        self.m.add_row(site=0,node=3,derived_state='1')

        msprime.sort_tables(nodes=self.n,edges=self.e,
                sites=self.s,mutations=self.m)
        idmap = msprime.simplify_tables(nodes=self.n,edges=self.e,
                sites=self.s,mutations=self.m,samples=[1,2])
        ts = msprime.load_tables(nodes=self.n,edges=self.e,sites=self.s,
                mutations=self.m)
        m = ts.genotype_matrix()
        self.assertEqual(m[0:].sum(),2)

    def test6(self):
        self.e.add_row(left=0,right=0.4,parent=0,child=1)
        self.e.add_row(left=0.6,right=0.8,parent=0,child=1)
        self.e.add_row(left=0.,right=1.,parent=0,child=2)
        self.e.add_row(left=0.,right=0.4,parent=3,child=0)
        self.e.add_row(left=0.4,right=0.8,parent=3,child=0)
        self.n.set_columns(time=[1,0,0,2],flags=[msprime.NODE_IS_SAMPLE]*4)

        self.s.add_row(position=0.4,ancestral_state='0')
        self.m.add_row(site=0,node=3,derived_state='1')

        msprime.sort_tables(nodes=self.n,edges=self.e,
                sites=self.s,mutations=self.m)
        idmap = msprime.simplify_tables(nodes=self.n,edges=self.e,
                sites=self.s,mutations=self.m,samples=[1,2])
        ts = msprime.load_tables(nodes=self.n,edges=self.e,sites=self.s,
                mutations=self.m)
        m = ts.genotype_matrix()
        self.assertEqual(m[0:].sum(),1)

if __name__=="__main__":
    unittest.main();
