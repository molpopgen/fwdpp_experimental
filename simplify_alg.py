"""
Python implementation of Algorithm S.

The example works by generating an initial TreeSequence
for a sample of 10 haplotypes using msprime.  We then
simplify the node/edge table in that TreeSequence with
respect to the first three samples.

This is Jerome Kelleher's implementation.  File is copy/pasted
from petrelharp/ftprime_ms
"""
import heapq
import numpy as np

import msprime


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


def simplify(S, Ni, Ei, L):
    """
    This is an implementation of the simplify algorithm described in Appendix A
    of the paper.
    """
    No = msprime.NodeTable()
    Eo = msprime.EdgeTable()
    A = [[] for _ in range(len(Ni))]
    Q = []

    ancient_nodes = []
    for u in S:
        v = No.add_row(time=Ni.time[u], flags=1)
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
        while ei < len(Ei) and Ei.parent[ei] == u:
            e = Ei[ei]
            for x in A[e.child]:
                if x.right > e.left and e.right > x.left:
                    y = Segment(max(x.left, e.left), min(
                        x.right, e.right), x.node)
                    heapq.heappush(Q, y)
            ei += 1

        v = -1
        while len(Q) > 0:
            l = Q[0].left
            r = L
            X = []
            while len(Q) > 0 and Q[0].left == l:
                x = heapq.heappop(Q)
                X.append(x)
                r = min(r, x.right)
            if len(Q) > 0:
                r = min(r, Q[0].left)

            if len(X) == 1:
                x = X[0]
                alpha = x
                if len(Q) > 0 and Q[0].left < x.right:
                    alpha = Segment(x.left, Q[0].left, x.node)
                    x.left = Q[0].left
                    heapq.heappush(Q, x)
            else:
                if v == -1:
                    v = No.add_row(time=Ni.time[u])
                alpha = Segment(l, r, v)
                for x in X:
                    Eo.add_row(l, r, v, x.node)
                    if x.right > r:
                        x.left = r
                        heapq.heappush(Q, x)

            A[u].append(alpha)

    # Sort the output edges and compact them as much as possible into
    # the output table. We skip this for the algorithm listing as it's pretty mundane.
    # TODO replace this with a calls to squash_edges() and sort_tables()
    E = list(Eo)
    Eo.clear()
    E.sort(key=lambda e: (e.parent, e.child, e.right, e.left))
    start = 0
    for j in range(1, len(E)):
        condition = (
            E[j - 1].right != E[j].left or
            E[j - 1].parent != E[j].parent or
            E[j - 1].child != E[j].child)
        if condition:
            Eo.add_row(E[start].left, E[j - 1].right,
                       E[j - 1].parent, E[j - 1].child)
            start = j
    j = len(E)
    Eo.add_row(E[start].left, E[j - 1].right, E[j - 1].parent, E[j - 1].child)

    # for i in Eo:
    #     print(i.left, i.right, i.parent, i.child,
    #           No.time[i.parent], No.time[i.child])
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

    # Generate initial TreeSequence
    ts = msprime.simulate(5000, recombination_rate=10, random_seed=1)
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
    sample = [0, 1, 2, 19, 33, 11, 12]
    ts1 = simplify(sample, nodes, edges, ts.sequence_length)

    msts = ts.simplify(sample)


    for i,j, in zip(msts.tables.edges, ts1.tables.edges):
        assert i.parent == j.parent, "parent error {} {}".format(i.parent,j.parent)
        assert i.child == j.child, "child error"
        assert i.left == j.left, "left error"
        assert i.right == j.right, "right error"

    # print(ts1.tables.nodes)
    # print(ts1.tables.edges)
