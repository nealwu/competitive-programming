#include <iostream>
using namespace std;

int64_t floor_div(int64_t a, int64_t b) {
    return a / b - ((a ^ b) < 0 && a % b != 0);
}

int64_t ceil_div(int64_t a, int64_t b) {
    return a / b + ((a ^ b) > 0 && a % b != 0);
}


int main() {
    int64_t a, b;

    while (cin >> a >> b)
        cout << floor_div(a, b) << ' ' << ceil_div(a, b) << '\n';
}
