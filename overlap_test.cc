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

class segment_overlapper
{
  private:
    vector<segment>::const_iterator sbeg, send;

    inline double
    set_partition()
    {
        double tright = numeric_limits<double>::max();
        auto b = overlapping.begin();
        for (auto i = overlapping.begin(); i < overlapping_end; ++i)
            {
                if (i->right > left)
                    {
                        *b = *i;
                        tright = min(tright, b->right);
                        ++b;
                    }
            }
        overlapping_end = b;
        return tright;
    }

    //inline double
    //min_right_overlap()
    //{
    //    return min_element(overlapping.begin(), overlapping_end,
    //                       [](const segment& a, const segment& b) {
    //                           return a.right < b.right;
    //                       })
    //        ->right;
    //}

  public:
    vector<segment> overlapping;
    vector<segment>::iterator overlapping_end;
    double left, right;
    segment_overlapper(const vector<segment>& segs)
        // The - 1 for send assumes a "cap"/sentinel value.
        : sbeg(segs.begin()), send(segs.end() - 1), overlapping{},
          overlapping_end(overlapping.end()), left(0),
          right(numeric_limits<double>::max())
    {
    }
    bool
    operator()()
    {
        bool rv = 0;
        if (sbeg < send)
            {
                left = right;
                auto tright = set_partition();
                if (num_overlaps() == 0)
                    {
                        left = sbeg->left;
                    }
                while (sbeg < send && sbeg->left == left)
                    {
                        tright = std::min(tright, sbeg->right);
                        overlapping_end
                            = overlapping.insert(overlapping_end, *sbeg) + 1;
                        ++sbeg;
                    }
                right = min(sbeg->left, tright);
                rv = true;
            }
        else
            {
                left = right;
                right = numeric_limits<double>::max();
                auto tright = set_partition();
                if (num_overlaps() > 0)
                    {
                        right = tright;
                        rv = true;
                    }
            }
        return rv;
    }
    std::int64_t
    num_overlaps()
    {
        return std::distance(overlapping.begin(), overlapping_end);
    }
};

void
calculate_overlaps(vector<segment>& segs, char* outfile)
{
    sort(segs.begin(), segs.end(),
         [](segment& a, segment& b) { return a.left < b.left; });
    size_t nsegs = segs.size();
    size_t index = 0;
    vector<segment> overlapping;
    vector<segment>::iterator overlapping_end = overlapping.end();
    // "cap" segs for convenienct
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
            if (distance(overlapping.begin(), overlapping_end) == 0)
                {
                    left = segs[index].left;
                }
            while (index < nsegs && segs[index].left == left)
                {
                    overlapping_end
                        = overlapping.insert(overlapping_end, segs[index]) + 1;
                    ++index;
                }
            right = min(segs[index].left,
                        min_element(overlapping.begin(), overlapping_end,
                                    [](const segment& a, const segment& b) {
                                        return a.right < b.right;
                                    })
                            ->right);
            assert(left < right);
            cout << left << ' ' << right << "-> ";
            for (auto j = overlapping.begin(); j < overlapping_end; ++j)
                {
                    cout << j->left << ',' << j->right << ',' << j->node
                         << " | ";
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
    sort(segs.begin(), segs.end(),
         [](segment& a, segment& b) { return a.left < b.left; });
    segs.emplace_back(segment{ numeric_limits<double>::max(),
                               numeric_limits<double>::max(), -1 });
    segment_overlapper o(segs);
    ofstream out(outfile);
    out.setf(std::ios_base::fixed, std::ios_base::floatfield);
    out.precision(6);
    while (o())
        {
            out << o.left << ' ' << o.right << "-> ";
            for (auto i = o.overlapping.begin(); i < o.overlapping_end; ++i)
                {
                    out << i->left << ',' << i->right << ',' << i->node
                         << " |";
                }
            out << '\n';
        }
    //calculate_overlaps(segs, outfile);
}
