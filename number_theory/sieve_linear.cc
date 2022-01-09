#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

vector<int> smallest_factor;
vector<bool> prime;
vector<int> primes;

// Note: this sieve is O(n), but the constant factor is worse than the O(n log log n) sieve due to the multiplication.
void sieve(int maximum) {
    maximum = max(maximum, 1);
    smallest_factor.assign(maximum + 1, 0);
    prime.assign(maximum + 1, true);
    prime[0] = prime[1] = false;
    primes = {};

    for (int i = 2; i <= maximum; i++) {
        if (prime[i]) {
            smallest_factor[i] = i;
            primes.push_back(i);
        }

        for (int p : primes) {
            if (p > smallest_factor[i] || int64_t(i) * p > maximum)
                break;

            prime[i * p] = false;
            smallest_factor[i * p] = p;
        }
    }
}


#include <numeric>

int main() {
    int n;
    scanf("%d", &n);
    sieve(n);
    printf("%lld = sum of primes\n", accumulate(primes.begin(), primes.end(), 0LL));
    printf("%d = number of primes\n", accumulate(prime.begin(), prime.end(), 0));
    printf("%lld = sum of smallest_factor\n", accumulate(smallest_factor.begin(), smallest_factor.end(), 0LL));
}
