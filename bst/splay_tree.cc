#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>
using namespace std;

// TODO: if a map is needed, search for instances of "splay_key" and add "splay_value" as well where appropriate.
using splay_key = int;

struct splay_node {
    splay_node *parent = nullptr, *child[2] = {nullptr, nullptr};
    splay_key key;
    int size = 1;

    // int64_t sum = 0;

    friend int get_size(splay_node *x) {
        return x == nullptr ? 0 : x->size;
    }

    // friend int64_t get_sum(splay_node *x) {
    //     return x == nullptr ? 0 : x->sum;
    // }

    int parent_index() const {
        return parent == nullptr ? -1 : int(this == parent->child[1]);
    }

    void set_child(int index, splay_node *x) {
        child[index] = x;

        if (x != nullptr)
            x->parent = this;
    }

    void join() {
        size = get_size(child[0]) + get_size(child[1]) + 1;
        // TODO: Add other operations here. Remember to store my_value in addition to, e.g., min_value or max_value.
        // sum = key + get_sum(child[0]) + get_sum(child[1]);
    }
};

auto random_address = [] { char *p = new char; delete p; return uint64_t(p); };
mt19937_64 splay_rng(chrono::steady_clock::now().time_since_epoch().count() * (random_address() | 1));
int64_t splay_count = 0;

struct splay_tree {
    static const int POOL_SIZE = 10000;
    static vector<splay_node*> node_pool;
    static vector<splay_node*> pointers_to_delete;

    static splay_node* new_node(const splay_key &key) {
        if (node_pool.empty()) {
            splay_node *ptr = new splay_node[POOL_SIZE];
            pointers_to_delete.push_back(ptr);
            node_pool.reserve(POOL_SIZE);

            for (int i = POOL_SIZE - 1; i >= 0; i--)
                node_pool.push_back(ptr + i);
        }

        splay_node *node = node_pool.back();
        node_pool.pop_back();
        *node = splay_node();
        node->key = key;
        // TODO: add more logic for creating a new node here.
        node->join();
        return node;
    }

    static bool _exit_delete_setup;

    static void _delete_pointers() {
        for (splay_node *node : pointers_to_delete)
            delete[] node;

        pointers_to_delete.clear();
    }

    ~splay_tree() {
        if (!_exit_delete_setup) {
            atexit(_delete_pointers);
            _exit_delete_setup = true;
        }
    }

    splay_node *root = nullptr;

    splay_tree(const vector<splay_key> &keys = {}) {
        init(keys);
    }

    splay_node *construct(const vector<splay_key> &keys, int start, int end) {
        if (start >= end)
            return nullptr;

        if (end - start == 1)
            return new_node(keys[start]);

        int mid = (start + end) / 2;
        splay_node *current = new_node(keys[mid]);
        current->set_child(0, construct(keys, start, mid));
        current->set_child(1, construct(keys, mid + 1, end));
        current->join();
        return current;
    }

    // Constructs the splay tree in linear time if the keys are already sorted.
    void init(vector<splay_key> keys) {
        if (!is_sorted(keys.begin(), keys.end()))
            sort(keys.begin(), keys.end());

        set_root(construct(keys, 0, int(keys.size())));
    }

    bool empty() const {
        return root == nullptr;
    }

    int size() const {
        return get_size(root);
    }

    splay_node *set_root(splay_node *x) {
        if (x != nullptr)
            x->parent = nullptr;

        return root = x;
    }

    void rotate_up(splay_node *x, bool x_join = true) {
        splay_node *p = x->parent, *gp = p->parent;
        int index = x->parent_index();

        if (gp == nullptr)
            set_root(x);
        else
            gp->set_child(p->parent_index(), x);

        p->set_child(index, x->child[!index]);
        x->set_child(!index, p);
        p->join();

        if (x_join)
            x->join();
    }

    void splay(splay_node *x) {
        splay_count++;

        while (x != root) {
            splay_node *p = x->parent;

            if (p != root)
                rotate_up(x->parent_index() == p->parent_index() ? p : x, false);

            rotate_up(x, false);
        }

        x->join();
    }

    static constexpr double SPLAY_RNG_RANGE = double(splay_rng.max()) + 1.0;
    static constexpr double LOG_CONSTANT = 2.0;
    static constexpr double SPLAY_PROB = 0.25;
    static const int SIZE_CUTOFF = 400;

    void check_splay(splay_node *x, int depth) {
        assert(x != nullptr);
        int n = size(), log_n = 32 - __builtin_clz(n);

        // Splay when deep or with a certain random chance when small.
        if (depth > LOG_CONSTANT * log_n || (n <= SIZE_CUTOFF && double(splay_rng()) < SPLAY_PROB * SPLAY_RNG_RANGE))
            splay(x);
    }

    pair<splay_node*, int> insert(const splay_key &key, bool require_unique = false) {
        return insert(new_node(key), require_unique);
    }

    // Returns {new node pointer, index (number of existing elements that are strictly less)}
    pair<splay_node*, int> insert(splay_node *x, bool require_unique = false) {
        if (root == nullptr)
            return {set_root(x), 0};

        splay_node *current = root, *previous = nullptr;
        int below = 0, depth = 0;

        while (current != nullptr) {
            previous = current;
            depth++;

            if (current->key < x->key) {
                below += get_size(current->child[0]) + 1;
                current = current->child[1];
            } else {
                // Check for equal keys.
                if (require_unique && !(x->key < current->key)) {
                    below += get_size(current->child[0]);
                    check_splay(current, depth);
                    return {current, below};
                }

                current = current->child[0];
            }
        }

        previous->set_child(int(previous->key < x->key), x);
        check_splay(x, depth);

        for (splay_node *node = x; node != nullptr; node = node->parent)
            node->join();

        return {x, below};
    }

    splay_node *begin() {
        if (root == nullptr)
            return nullptr;

        splay_node *x = root;
        int depth = 0;

        while (x->child[0] != nullptr) {
            x = x->child[0];
            depth++;
        }

        check_splay(x, depth);
        return x;
    }

    // To iterate through all nodes in order:
    // for (splay_node *node = tree.begin(); node != nullptr; node = tree.successor(node))
    splay_node *successor(splay_node *x) const {
        if (x == nullptr)
            return nullptr;

        if (x->child[1] != nullptr) {
            x = x->child[1];

            while (x->child[0] != nullptr)
                x = x->child[0];

            return x;
        }

        while (x->parent_index() == 1)
            x = x->parent;

        return x->parent;
    }

    splay_node *predecessor(splay_node *x) const {
        if (x == nullptr)
            return nullptr;

        if (x->child[0] != nullptr) {
            x = x->child[0];

            while (x->child[1] != nullptr)
                x = x->child[1];

            return x;
        }

        while (x->parent_index() == 0)
            x = x->parent;

        return x->parent;
    }

    splay_node *last() {
        if (root == nullptr)
            return nullptr;

        splay_node *x = root;
        int depth = 0;

        while (x->child[1] != nullptr) {
            x = x->child[1];
            depth++;
        }

        check_splay(x, depth);
        return x;
    }

    void clear() {
        vector<splay_node*> nodes;
        nodes.reserve(size());

        for (splay_node *node = begin(); node != nullptr; node = successor(node))
            nodes.push_back(node);

        for (splay_node *node : nodes) {
            // Instead of deleting, add `node` back to `node_pool`.
            *node = splay_node();
            node_pool.push_back(node);
        }

        set_root(nullptr);
    }

    void erase(splay_node *x) {
        splay_node *new_x = nullptr, *fix_node = nullptr;

        if (x->child[0] == nullptr || x->child[1] == nullptr) {
            new_x = x->child[int(x->child[0] == nullptr)];
            fix_node = x->parent;
        } else {
            splay_node *next = successor(x);
            assert(next != nullptr && next->child[0] == nullptr);
            new_x = next;
            fix_node = next->parent == x ? next : next->parent;

            next->parent->set_child(next->parent_index(), next->child[1]);
            next->set_child(0, x->child[0]);
            next->set_child(1, x->child[1]);
        }

        if (x == root)
            set_root(new_x);
        else
            x->parent->set_child(x->parent_index(), new_x);

        int depth = 0;

        for (splay_node *node = fix_node; node != nullptr; node = node->parent) {
            node->join();
            depth++;
        }

        if (fix_node != nullptr)
            check_splay(fix_node, depth);

        // Instead of deleting, add `x` back to `node_pool`.
        *x = splay_node();
        node_pool.push_back(x);
    }

    // Returns {node pointer, index (number of existing elements that are strictly less)}
    pair<splay_node*, int> lower_bound(const splay_key &key) {
        splay_node *current = root, *previous = nullptr, *answer = nullptr;
        int below = 0, depth = 0;

        while (current != nullptr) {
            previous = current;
            depth++;

            if (current->key < key) {
                below += get_size(current->child[0]) + 1;
                current = current->child[1];
            } else {
                answer = current;
                current = current->child[0];
            }
        }

        if (previous != nullptr)
            check_splay(previous, depth);

        return make_pair(answer, below);
    }

    bool contains(const splay_key &key) {
        splay_node *node = lower_bound(key).first;
        return node != nullptr && node->key == key;
    }

    bool erase(const splay_key &key) {
        splay_node *x = lower_bound(key).first;

        if (x == nullptr || x->key != key)
            return false;

        erase(x);
        return true;
    }

    splay_node *node_at_index(int index) {
        if (index < 0 || index >= size())
            return nullptr;

        splay_node *current = root;
        int depth = 0;

        while (current != nullptr) {
            int left_size = get_size(current->child[0]);
            depth++;

            if (index == left_size) {
                check_splay(current, depth);
                return current;
            }

            if (index < left_size) {
                current = current->child[0];
            } else {
                current = current->child[1];
                index -= left_size + 1;
            }
        }

        assert(false);
    }

    // Returns a splay_node pointer representing all nodes with key less than `key`. If none, returns `nullptr`.
    splay_node *query_prefix_key(const splay_key &key) {
        splay_node *node = lower_bound(key).first;

        if (node == nullptr)
            return root;

        splay(node);
        return node->child[0];
    }

    // Returns a splay_node pointer representing the first `count` nodes. If none, returns `nullptr`.
    splay_node *query_prefix_count(int count) {
        if (count <= 0)
            return nullptr;
        else if (count >= size())
            return root;

        splay_node *node = node_at_index(count);
        splay(node);
        return node->child[0];
    }

    // Returns a splay_node pointer representing all nodes with key >= `key`. If none, returns `nullptr`.
    splay_node *query_suffix_key(const splay_key &key) {
        splay_node *node = lower_bound(key).first;

        if (node == nullptr)
            return nullptr;

        node = predecessor(node);

        if (node == nullptr)
            return root;

        splay(node);
        return node->child[1];
    }

    // Returns a splay_node pointer representing the last `count` nodes. If none, returns `nullptr`.
    splay_node *query_suffix_count(int count) {
        if (count >= size())
            return root;
        else if (count <= 0)
            return nullptr;

        int index = size() - count;
        splay_node *node = node_at_index(index - 1);
        splay(node);
        return node->child[1];
    }

    // Returns a splay_node pointer representing the range [start, end). If none, returns `nullptr`.
    splay_node *query_range(int start, int end) {
        if (start >= end)
            return nullptr;
        else if (start <= 0)
            return query_prefix_count(end);
        else if (end >= size())
            return query_suffix_count(size() - start);

        splay_node *before = node_at_index(start - 1);
        splay_node *after = node_at_index(end);
        splay(after);
        splay(before);

        if (after->parent != before)
            rotate_up(after);

        assert(before->child[1] == after);
        return after->child[0];
    }

    // Returns a splay_node pointer representing the range whose key is in [lower, upper). If none, returns `nullptr`.
    splay_node *query_range_key(const splay_key &lower, const splay_key &upper) {
        return query_range(lower_bound(lower).second, lower_bound(upper).second);
    }
};

bool splay_tree::_exit_delete_setup = false;
vector<splay_node*> splay_tree::node_pool;
vector<splay_node*> splay_tree::pointers_to_delete;



void print_tree(splay_node *x, int depth = 0) {
    cerr << string(depth, ' ');

    if (x == nullptr) {
        cerr << "null" << endl;
        return;
    }

    cerr << x->key << " (" << x->size << ")\n";

    if (x->child[0] != nullptr || x->child[1] != nullptr) {
        print_tree(x->child[0], depth + 1);
        print_tree(x->child[1], depth + 1);
    }
}

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    splay_tree tree;

    int require_unique;
    cin >> require_unique;

    string task;
    int x;

    while (cin >> task >> x) {
        if (task == "insert") {
            cout << tree.insert(x, require_unique).second << '\n';
        } else if (task == "index") {
            splay_node *node = tree.node_at_index(x);

            if (node == nullptr)
                cout << "none" << '\n';
            else
                cout << node->key << '\n';

            auto &&normalize_size = [&](int size) {
                return min(max(size, 0), tree.size());
            };

            assert(get_size(tree.query_prefix_count(x)) == normalize_size(x));
            assert(get_size(tree.query_suffix_count(x)) == normalize_size(x));
        } else if (task == "less_than") {
            int answer = tree.lower_bound(x).second;
            cout << answer << '\n';

            assert(get_size(tree.query_prefix_key(x)) == answer);
            assert(get_size(tree.query_suffix_key(x)) == tree.size() - answer);
        } else if (task == "erase") {
            cout << tree.erase(x) << '\n';
        } else {
            assert(false);
        }

        // cerr << "tree:" << endl;
        // print_tree(tree.root);
    }

    int start = int(splay_rng() % (tree.size() + 1));
    int end = int(splay_rng() % (tree.size() + 1));

    if (start > end)
        swap(start, end);

    assert(get_size(tree.query_range(start, end)) == end - start);
    splay_node *start_node = tree.node_at_index(start);
    splay_node *end_node = tree.node_at_index(end);
    splay_key lower = start_node == nullptr ? numeric_limits<splay_key>::max() : start_node->key;
    splay_key upper = end_node == nullptr ? numeric_limits<splay_key>::max() : end_node->key;
    int range_size = get_size(tree.query_range_key(lower, upper));
    assert(range_size == tree.lower_bound(upper).second - tree.lower_bound(lower).second);
    cerr << splay_count << " splays" << '\n';
}
