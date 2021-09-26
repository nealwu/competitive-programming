// Solution to https://www.spoj.com/problems/ORDERSET/
#include <iostream>
using namespace std;

#include <ext/pb_ds/assoc_container.hpp>
using namespace __gnu_pbds;

// WARNING: functions as a set (doesn't allow duplicates); insert pairs instead if duplicates are needed.
// Consider using splay_tree instead if constant factor is an issue (e.g., log^2 solutions), especially with duplicates.
template<typename T>
using ordered_set = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;

int main() {
    int Q;
    scanf("%d", &Q);
    ordered_set<int> values;

    for (int q = 0; q < Q; q++) {
        char op;
        int x;
        scanf(" %c %d", &op, &x);

        if (op == 'I') {
            values.insert(x);
        } else if (op == 'D') {
            values.erase(x);
        } else if (op == 'K') {
            x--;

            if (x >= int(values.size()))
                puts("invalid");
            else
                // find_by_order(x) gives the x-th element in sorted order; if x = 2, gives the third smallest value.
                printf("%d\n", *values.find_by_order(x));
        } else if (op == 'C') {
            // order_of_key(x) returns the count of elements < x.
            printf("%d\n", int(values.order_of_key(x)));
        }
    }
}
