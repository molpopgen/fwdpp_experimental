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
    size_t noverlapping = 0, index = 0, j, k;
    vector<segment*> overlapping;
    segs.emplace_back(segment{ numeric_limits<double>::max(),
                               numeric_limits<double>::max(), 0 });
    double right = segs[0].left;
    double left = right;
    while (index < nsegs)
        {
            left = right;
            k = 0;
            for (j = 0; j < noverlapping; ++j)
                {
                    if (overlapping[j]->right > left)
                        {
                            overlapping[k] = overlapping[j];
                            ++k;
                        }
                }
            noverlapping = k;
            if (k == 0)
                {
                    left = segs[index].left;
                }
            while (index < nsegs && segs[index].left == left)
                {
                    overlapping.insert(overlapping.begin() + noverlapping,
                                       &segs[index]);
                    ++noverlapping;
                    ++index;
                }
            right = segs[index].left;
            for (j = 0; j < noverlapping; ++j)
                {
                    right = min(right, overlapping[j]->right);
                }
            assert(left < right);
            cout << left << ' ' << right << "-> ";
            for (j = 0; j < noverlapping; ++j)
                {
                    cout << overlapping[j]->left << ','
                         << overlapping[j]->right << ','
                         << overlapping[j]->node << " | ";
                }
            cout << '\n';
        }
    while (noverlapping > 0)
        {
            left = right;
            k = 0;
            for (j = 0; j < noverlapping; ++j)
                {
                    if (overlapping[j]->right > left)
                        {
                            overlapping[k] = overlapping[j];
                            ++k;
                        }
                }
            noverlapping = k;
            if (noverlapping > 0)
                {
                    right = numeric_limits<double>::max();
                    for (j = 0; j < noverlapping; ++j)
                        {
                            right = min(right, overlapping[j]->right);
                        }
                    cout << left << ' ' << right << "-> ";
                    for (j = 0; j < noverlapping; ++j)
                        {
                            cout << overlapping[j]->left << ','
                                 << overlapping[j]->right << ','
                                 << overlapping[j]->node << " | ";
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
