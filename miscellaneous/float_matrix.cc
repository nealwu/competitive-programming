#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

// TODO: switch to `double` if `long double` is unnecessary and the time limit is tight.
// Using `long double` is more accurate, but it can be 50-60% slower than `double`.
using matrix_float = long double;

// TODO: if using float_column_vector, we can write the float_matrix in the format matrix[x] = a row of coefficients
// used to build the x-th element of the float_column_vector. So matrix[0][2] is the coefficient that element 2
// contributes to the next element 0.
// The other option is to take a single-row 1 * n float_matrix and multiply it by the n * n float_matrix. Then
// matrix[0][2] is the coefficient that 0 contributes to the next element 2.
struct float_column_vector {
    int rows;
    vector<matrix_float> values;

    float_column_vector(int _rows = 0) {
        init(_rows);
    }

    template<typename T>
    float_column_vector(const vector<T> &v) {
        init(v);
    }

    void init(int _rows) {
        rows = _rows;
        values.assign(rows, 0);
    }

    template<typename T>
    void init(const vector<T> &v) {
        rows = int(v.size());
        values = vector<matrix_float>(v.begin(), v.end());
    }

    matrix_float& operator[](int index) { return values[index]; }
    const matrix_float& operator[](int index) const { return values[index]; }
};

// Warning: very inefficient for many small matrices of fixed size. For that, use float_matrix_fixed_size.cc instead.
struct float_matrix {
    static float_matrix IDENTITY(int n) {
        float_matrix identity(n);

        for (int i = 0; i < n; i++)
            identity[i][i] = 1;

        return identity;
    }

    int rows, cols;
    vector<vector<matrix_float>> values;

    float_matrix(int _rows = 0, int _cols = -1) {
        init(_rows, _cols);
    }

    template<typename T>
    float_matrix(const vector<vector<T>> &v) {
        init(v);
    }

    void init(int _rows, int _cols = -1) {
        rows = _rows;
        cols = _cols < 0 ? rows : _cols;
        values.assign(rows, vector<matrix_float>(cols, 0));
    }

    template<typename T>
    void init(const vector<vector<T>> &v) {
        rows = int(v.size());
        cols = v.empty() ? 0 : int(v[0].size());
        values.assign(rows, vector<matrix_float>(cols, 0));

        for (int i = 0; i < rows; i++) {
            assert(int(v[i].size()) == cols);
            copy(v[i].begin(), v[i].end(), values[i].begin());
        }
    }

    vector<matrix_float>& operator[](int index) { return values[index]; }
    const vector<matrix_float>& operator[](int index) const { return values[index]; }

    bool is_square() const {
        return rows == cols;
    }

    float_matrix operator*(const float_matrix &other) const {
        assert(cols == other.rows);
        float_matrix product(rows, other.cols);

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                if (values[i][j] != 0)
                    for (int k = 0; k < other.cols; k++)
                        product[i][k] += values[i][j] * other[j][k];

        return product;
    }

    float_matrix& operator*=(const float_matrix &other) {
        return *this = *this * other;
    }

    float_column_vector operator*(const float_column_vector &column) const {
        assert(cols == column.rows);
        float_column_vector product(rows);

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                product[i] += values[i][j] * column[j];

        return product;
    }

    float_matrix& operator*=(matrix_float mult) {
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                values[i][j] *= mult;

        return *this;
    }

    float_matrix operator*(matrix_float mult) const {
        return float_matrix(*this) *= mult;
    }

    float_matrix& operator+=(const float_matrix &other) {
        assert(rows == other.rows && cols == other.cols);

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                values[i][j] += other[i][j];

        return *this;
    }

    float_matrix operator+(const float_matrix &other) const {
        return float_matrix(*this) += other;
    }

    float_matrix& operator-=(const float_matrix &other) {
        assert(rows == other.rows && cols == other.cols);

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                values[i][j] -= other[i][j];

        return *this;
    }

    float_matrix operator-(const float_matrix &other) const {
        return float_matrix(*this) -= other;
    }

    float_matrix pow(int64_t p) const {
        assert(p >= 0);
        assert(is_square());
        float_matrix m = *this, result = IDENTITY(rows);

        while (p > 0) {
            if (p & 1)
                result *= m;

            p >>= 1;

            if (p > 0)
                m *= m;
        }

        return result;
    }

    void print(ostream &os) const {
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                os << values[i][j] << (j < cols - 1 ? ' ' : '\n');

        os << '\n';
    }
};


#include <iomanip>

void read_matrix(float_matrix &m) {
    int r, c;
    cin >> r >> c;
    m = float_matrix(r, c);

    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++) {
            double x;
            cin >> x;
            m[i][j] = x;
        }
}

int main() {
    cout << setprecision(16);

    float_matrix m1, m2;
    read_matrix(m1);
    read_matrix(m2);
    (m1 + m1).print(cout);
    (m2 - m2).print(cout);
    (m1 * m2).print(cout);

    read_matrix(m1);
    int64_t p;
    cin >> p;
    (m1 * p).print(cout);
    m1.pow(p).print(cout);
}
