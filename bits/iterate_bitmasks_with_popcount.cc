#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

template<typename F>
void iterate_bitmasks_with_popcount(int n, int k, F &&f) {
    if (k == 0) {
        f(0);
        return;
    }

    int64_t mask = (1LL << k) - 1;

    while (mask < 1LL << n) {
        f(mask);
        int zeros = __builtin_ctzll(mask);
        int ones = __builtin_ctzll(~mask >> zeros);
        mask += (1LL << zeros) + (1LL << (ones - 1)) - 1;
    }
}

int main() {
    for (int n = 0; n <= 20; n++) {
        vector<vector<int>> masks_by_count(n + 1);

        for (int mask = 0; mask < 1 << n; mask++)
            masks_by_count[__builtin_popcount(mask)].push_back(mask);

        for (int k = 0; k <= n; k++) {
            vector<int> k_masks;

            iterate_bitmasks_with_popcount(n, k, [&](int64_t mask) -> void {
                k_masks.push_back(int(mask));
            });

            assert(k_masks == masks_by_count[k]);
        }
    }

    const int maximum = 64;
    vector<vector<uint64_t>> choose(maximum + 1);

    for (int n = 0; n <= maximum; n++) {
        choose[n].resize(n + 1);
        choose[n][0] = choose[n][n] = 1;

        for (int r = 1; r < n; r++)
            choose[n][r] = choose[n - 1][r - 1] + choose[n - 1][r];
    }

    auto test_masks = [&](int n, int k) -> void {
        int64_t previous = -1;
        uint64_t count = 0;

        iterate_bitmasks_with_popcount(n, k, [&](int64_t mask) -> void {
            assert(__builtin_popcountll(mask) == k);
            assert(mask > previous);
            previous = mask;
            count++;
        });

        assert(count == choose[n][k]);
    };

    test_masks(60, 5);
    test_masks(60, 55);
}
