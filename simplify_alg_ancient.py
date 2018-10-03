"""
Updated Python implementation of what is 
presented in the Kelleher et al paper.
Logic updated to msprime 0.6.x and 
ancient samples are handled.
"""
import heapq
import numpy as np
import sys
import msprime
import simplify as mspsimplify
from collections import namedtuple

#Edge = namedtuple('Edge', ['left', 'right', 'parent', 'child'])


class Segment(object):
    """
    An ancestral segment mapping a given node ID to a half-open genomic
    interval [left, right).
    """

    def __init__(self, left, right, node):
        assert left < right
        self.left = left
        self.right = right
        self.node = node

    def __lt__(self, other):
        # We implement the < operator so that we can use the heapq directly.
        # We are only interested in the left coordinate.
        return self.left < other.left

    def __repr__(self):
        return repr((self.left, self.right, self.node))


def overlapping_segments(segments):
    """
    Returns an iterator over the (left, right, X) tuples describing the
    distinct overlapping segments in the specified set.
    """
    S = sorted(segments, key=lambda x: x.left)
    print(S)
    n = len(S)
    # Insert a sentinel at the end for convenience.
    S.append(Segment(sys.float_info.max/2, sys.float_info.max, 0))
    right = S[0].left
    X = []
    j = 0
    print('right = ',right)
    while j < n:
        # Remove any elements of X with right <= left
        left = right
        X = [x for x in X if x.right > left]
        print(left,X)
        if len(X) == 0:
            print('yes')
            left = S[j].left
        while j < n and S[j].left == left:
            X.append(S[j])
            j += 1
        j -= 1
        right = min(x.right for x in X)
        right = min(right, S[j + 1].left)
        print("yielding ",left,right,X)
        yield left, right, X
        j += 1

    while len(X) > 0:
        left = right
        X = [x for x in X if x.right > left]
        if len(X) > 0:
            right = min(x.right for x in X)
            print("yielding ",left,right,X)
            yield left, right, X


def add_ancestry(A, input_id, left, right, node):
    if len(A[input_id]) == 0:
        assert left<right,"{} {}".format(left,right)
        x = Segment(left, right, node)
        A[input_id].append(x)
    else:
        tail = A[input_id][-1]
        if tail.right == left and tail.node == node:
            tail.right = right
        else:
            assert left<right,"{} {}".format(left,right)
            x = Segment(left, right, node)
            A[input_id].append(x)


def flush_edges(temp_edges, Eo):
    n = 0
    for child in sorted(temp_edges.keys()):
        for edge in temp_edges[child]:
            Eo.add_row(edge.left,edge.right,edge.parent,edge.child)
            n+=1

    # ute = []
    # for e in temp_edges:
    #     if e not in ute:
    #         ute.append(e)
    # for e in ute:
    #     #print(e)
    #     Eo.add_row(left=e.left,right=e.right,parent=e.parent,child=e.child)
    #     n+=1
    return n

def buffer_edge(temp_edges,left,right,parent,child):
    if child not in temp_edges:
        temp_edges[child]=[msprime.Edge(left,right,parent,child)]
    else:
        l=temp_edges[child][-1]
        if l.right == left:
            l.right=right
        else:
            temp_edges[child].append(msprime.Edge(left,right,parent,child))


def simplify(S, Ni, Ei, L):
    """
    This is an implementation of the simplify algorithm described in Appendix A
    of the paper.
    """
    No = msprime.NodeTable()
    Eo = msprime.EdgeTable()
    A = [[] for _ in range(len(Ni))]
    Q = []

    idmap = [-1 for i in range(len(Ni))]
    ancient_nodes = []
    for u in S:
        v = No.add_row(time=Ni.time[u], flags=1)
        idmap[u] = v
        if Ni.time[u] != 0.0:
            ancient_nodes.append(u)
        assert(v == len(No)-1)
        A[u] = [Segment(0, L, v)]

    # for u in S:
    #     print(u, A[u])
    # print("ancient nodes = ", ancient_nodes)

    # These changes make sure that
    # we collect edges for merging
    # in proper time order.
    # inodes = [i for i in range(len(Ni))]
    # inodes = sorted(inodes,key=lambda x:Ni.time[x])
    # for u in range(len(Ni)):
    # for u in inodes:
    #     for e in [e for e in Ei if e.parent == u]:
    #         for x in A[e.child]:
    #             if x.right > e.left and e.right > x.left:
    #                 y = Segment(max(x.left, e.left), min(
    #                     x.right, e.right), x.node)
    #                 heapq.heappush(Q, y)
    ei = 0
    while ei < len(Ei):
        u = Ei.parent[ei]
        edges = []
        while ei < len(Ei) and Ei.parent[ei] == u:
            edges.append(Ei[ei])
            # for x in A[e.child]:
            #     if x.right > e.left and e.right > x.left:
            #         y = Segment(max(x.left, e.left), min(
            #             x.right, e.right), x.node)
            #         # print(y)
            #         segs.append(y)
            #         heapq.heappush(Q, y)
            ei += 1
        segs=[]
        for edge in edges:
            for x in A[edge.child]:     
                if x.right > edge.left and edge.right > x.left:
                    y = Segment(max(x.left, edge.left), min(x.right, edge.right), x.node)
                    segs.append(y)
        
        with open('me.txt','a') as f:
            f.write(str(segs))
            f.write('\n')
        prev_right = 0

        temp_edges = {}
        v = idmap[u]
        is_sample = v != -1
        if is_sample:
            x=A[u][-1]
            assert x.left==0 and x.right==L
            A[u]=[]
        for left, right, X in overlapping_segments(segs):
            if len(X) == 1:
                ancestry_node = X[0].node
                if is_sample:
                    assert left < right
                    buffer_edge(temp_edges,left,right,v,ancestry_node)
                    #temp_edges.append(msprime.Edge(left, right, v, ancestry_node))
                    ancestry_node = v
            else:
                if v == -1:
                    v = No.add_row(time=Ni.time[u])
                    idmap[u] = v
                ancestry_node = v
                for x in X:
                    assert left < right
                    buffer_edge(temp_edges,left,right,v,x.node)
                    #temp_edges.append(msprime.Edge(left, right, v, x.node))
            if is_sample and left != prev_right: 
                add_ancestry(A, u, prev_right, left, v)
            assert left<right
            assert ancestry_node!=-1
            add_ancestry(A, u, left, right, ancestry_node)
            prev_right = right
        if is_sample and prev_right != L:
            add_ancestry(A,u,prev_right,L,v)
        if v != -1:
            n = flush_edges(temp_edges, Eo)
            if n == 0 and not is_sample:
                No.truncate(v)
                idmap[u] = -1
    for i in Eo:
        assert i.parent < len(No), "bad parent index"
        assert i.child < len(No), "bad child index"
    return msprime.load_tables(nodes=No, edges=Eo)


def verify():
    """
    Checks that simplify() does the right thing, by comparing to the implementation
    in msprime.
    """
    for n in [10, 100, 1000]:
        ts = msprime.simulate(n, recombination_rate=1, random_seed=1)
        nodes = ts.tables.nodes
        edges = ts.tables.edges
        print("simulated for ", n)

        for N in range(2, 10):
            sample = list(range(N))
            ts1 = simplify(sample, nodes, edges, ts.sequence_length)
            ts2 = ts.simplify(sample)

            n1 = ts1.tables.nodes
            n2 = ts2.tables.nodes
            assert np.array_equal(n1.time, n2.time)
            assert np.array_equal(n1.flags, n2.flags)
            e1 = ts1.tables.edges
            e2 = ts2.tables.edges
            assert np.array_equal(e1.left, e2.left)
            assert np.array_equal(e1.right, e2.right)
            assert np.array_equal(e1.parent, e1.parent)
            assert np.array_equal(e1.child, e1.child)


if __name__ == "__main__":
    # verify()

    # np.random.seed(666)
    # Generate initial TreeSequence
    ts = msprime.simulate(5000, recombination_rate=10) #, random_seed=1)
    nodes = ts.tables.nodes
    edges = ts.tables.edges

    with open("test_nodes.txt", "w") as f:
        for i in range(len(nodes)):
            f.write("{} {}\n".format(i, nodes[i].time))

    with open("test_edges.txt", "w") as f:
        for i in range(len(edges)):
            f.write("{} {} {} {}\n".format(
                edges[i].parent, edges[i].child, edges[i].left, edges[i].right))

    # Simplify nodes and edges with respect to the following samples:
    sample = np.array(np.random.choice(len(ts.tables.nodes), 500,
                                     replace=False),dtype=np.int32)
    fc = np.zeros(len(ts.tables.nodes), dtype=np.uint32)
    fc[sample] = 1
    tt = ts.dump_tables()
    tt.nodes.set_columns(flags=fc, time=ts.tables.nodes.time)
    ts = ts.load_tables(tt)
    msts = ts.simplify(sample)
    S=mspsimplify.Simplifier(ts,sample)
    Sts, sidmap = S.simplify()

    ts1 = simplify(sample, nodes, edges, ts.sequence_length)

    for i, j, rec in zip(msts.tables.edges, ts1.tables.edges, range(len(ts1.tables.edges))):
        # print(i, j)
        assert i.parent == j.parent, "parent error {} {} {}".format(
            i.parent, j.parent, rec)
        assert i.child == j.child, "child error"
        assert i.left == j.left, "left error"
        assert i.right == j.right, "right error"

