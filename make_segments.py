import sys
import numpy as np
from collections import namedtuple

Segment = namedtuple('Segment', ['left', 'right', 'node'])


def overlapping_segments(segments):
    """
    Returns an iterator over the (left, right, X) tuples describing the
    distinct overlapping segments in the specified set.

    From msprime/tests/simplify.py by Jerome Kelleher
    """
    S = sorted(segments, key=lambda x: x.left)
    n = len(S)
    # Insert a sentinel at the end for convenience.
    S.append(Segment(sys.float_info.max, sys.float_info.max, 0))
    right = S[0].left
    X = []
    j = 0
    while j < n:
        # Remove any elements of X with right <= left
        left = right
        X = [x for x in X if x.right > left]
        if len(X) == 0:
            left = S[j].left
        while j < n and S[j].left == left:
            X.append(S[j])
            j += 1
        j -= 1
        right = min(x.right for x in X)
        right = min(right, S[j + 1].left)
        yield left, right, X
        j += 1

    while len(X) > 0:
        left = right
        X = [x for x in X if x.right > left]
        if len(X) > 0:
            right = min(x.right for x in X)
            yield left, right, X


nsegs = int(sys.argv[1])
ofname = sys.argv[2]
ofname2 = sys.argv[3]
segs = []
with open(ofname, 'w') as f:
    for i in range(nsegs):
        l = np.random.uniform(0, 1)
        r = np.random.uniform(l, 1)
        n = np.random.randint(0, 1000)
        f.write("{} {} {}\n".format(l, r, n))
        segs.append(Segment(l, r, n))

with open(ofname2, 'w') as f:
    for left, right, X in overlapping_segments(segs):
        #f.write("{} {} {}\n".format(left,right,X))
        f.write("{} {}-> ".format(left, right))
        for x in X:
            f.write("{},{},{} |".format(x.left, x.right, int(x.node)))
        f.write("\n")
