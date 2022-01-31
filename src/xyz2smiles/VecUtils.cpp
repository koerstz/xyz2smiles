#include <vector>
#include "VecUtils.h"

ListofListInts IterProduct(ListofListInts list) {
    // Not ideal makes all combinations in memory.
    ListofListInts result = {{}};
    for (std::vector<int> pool: list) {
        ListofListInts tmp = {};
        for (int y: pool) {
            for (std::vector<int> x: result) {
                x.push_back(y);
                tmp.push_back(x);
            }
        }
        result = tmp;
    }
    return result;
}


int SumVector (std::vector<int> values) {
    int s = 0;
    for (int i : values) {
        s += i;
    }
    return s;
}