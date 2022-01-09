#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

void output(int mask, int n) {
    for (int i = 0; i < n; i++)
        cout << (mask >> i & 1);

    cout << '\n';
}

int main() {
    int n, mask;
    cin >> n >> mask;

    for (int super = mask; super < 1 << n; super = (super + 1) | mask) {
        printf("%3d: ", super);
        output(super, n);
    }
}
