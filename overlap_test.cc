#include <vector>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <cassert>

using namespace std;

struct segment
{
    double left, right;
    int32_t node;
};

void
calculate_overlaps(vector<segment>& segs, char* outfile)
{
    sort(segs.begin(), segs.end(),
         [](segment& a, segment& b) { return a.left < b.left; });
    // "cap" segs for convenienct
    size_t nsegs = segs.size();
    size_t noverlapping = 0, index = 0, j;
    vector<segment> overlapping;
    vector<segment>::iterator overlapping_end = overlapping.end();
    segs.emplace_back(segment{ numeric_limits<double>::max(),
                               numeric_limits<double>::max(), 0 });
    double right = segs[0].left;
    double left = right;
    while (index < nsegs)
        {
            left = right;
            overlapping_end = stable_partition(
                overlapping.begin(), overlapping_end,
                [left](const segment& seg) { return seg.right > left; });
            noverlapping = distance(overlapping.begin(), overlapping_end);
            if (noverlapping == 0)
                {
                    left = segs[index].left;
                }
            while (index < nsegs && segs[index].left == left)
                {
                    overlapping_end
                        = overlapping.insert(
                              overlapping.begin() + noverlapping, segs[index])
                          + 1;
                    ++index;
                }
            noverlapping = distance(overlapping.begin(), overlapping_end);
            right = min(segs[index].left,
                        min_element(overlapping.begin(), overlapping_end,
                                    [](const segment& a, const segment& b) {
                                        return a.right < b.right;
                                    })
                            ->right);
            assert(left < right);
            cout << left << ' ' << right << "-> ";
            for (j = 0; j < noverlapping; ++j)
                {
                    cout << overlapping[j].left << ',' << overlapping[j].right
                         << ',' << overlapping[j].node << " | ";
                }
            cout << '\n';
        }
    while (overlapping_end > overlapping.begin())
        {
            left = right;
            overlapping_end = stable_partition(
                overlapping.begin(), overlapping_end,
                [left](const segment& seg) { return seg.right > left; });
            if (overlapping_end > overlapping.begin())
                {
                    right
                        = min_element(overlapping.begin(), overlapping_end,
                                      [](const segment& a, const segment& b) {
                                          return a.right < b.right;
                                      })
                              ->right;

                    cout << left << ' ' << right << "-> ";
                    for (auto i = overlapping.begin(); i < overlapping_end;
                         ++i)
                        {
                            cout << i->left << ',' << i->right << ','
                                 << i->node << " | ";
                        }
                    cout << '\n';
                }
        }
}

int
main(int argc, char** argv)
{
    int argn = 1;
    char* infile = argv[argn++];
    char* outfile = argv[argn++];

    ifstream in(infile);

    vector<segment> segs;
    double left, right;
    int32_t node;
    while (!in.eof())
        {
            in >> left >> right >> node >> ws;
            segs.emplace_back(segment{ left, right, node });
        }
    calculate_overlaps(segs, outfile);
}
