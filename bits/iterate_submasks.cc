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

    for (int sub = mask; ; sub = (sub - 1) & mask) {
        printf("%3d: ", sub);
        output(sub, n);
        if (sub == 0) break;
    }
}
