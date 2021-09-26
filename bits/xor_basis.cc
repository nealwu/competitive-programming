
const int BITS = 30;

template<typename T>
struct xor_basis {
    // A list of basis values sorted in decreasing order, where each value has a unique highest bit.
    // We use a static array instead of a vector for better performance.
    T basis[BITS];
    int n = 0;

    T min_value(T start) const {
        if (n == BITS)
            return 0;

        for (int i = 0; i < n; i++)
            start = min(start, start ^ basis[i]);

        return start;
    }

    T max_value(T start = 0) const {
        if (n == BITS)
            return (T(1) << BITS) - 1;

        for (int i = 0; i < n; i++)
            start = max(start, start ^ basis[i]);

        return start;
    }

    bool add(T x) {
        x = min_value(x);

        if (x == 0)
            return false;

        basis[n++] = x;
        int k = n - 1;

        // Insertion sort.
        while (k > 0 && basis[k] > basis[k - 1]) {
            swap(basis[k], basis[k - 1]);
            k--;
        }

        // Remove the highest bit of x from other basis elements.
        // TODO: this can be removed for speed if desired.
        for (int i = k - 1; i >= 0; i--)
            basis[i] = min(basis[i], basis[i] ^ x);

        return true;
    }

    void merge(const xor_basis<T> &other) {
        for (int i = 0; i < other.n && n < BITS; i++)
            add(other.basis[i]);
    }

    void merge(const xor_basis<T> &a, const xor_basis<T> &b) {
        if (a.n > b.n) {
            *this = a;
            merge(b);
        } else {
            *this = b;
            merge(a);
        }
    }
};
