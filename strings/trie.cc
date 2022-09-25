#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>
using namespace std;

template<char MIN_CHAR = 'a', int ALPHABET = 26>
struct array_trie {
    struct trie_node {
        array<int, ALPHABET> child;
        int words_here = 0, starting_with = 0;

        trie_node() {
            memset(&child[0], -1, ALPHABET * sizeof(int));
        }
    };

    static const int ROOT = 0;

    vector<trie_node> nodes = {trie_node()};

    array_trie(int total_length = -1) {
        if (total_length >= 0)
            nodes.reserve(total_length + 1);
    }

    int get_or_create_child(int node, int c) {
        if (nodes[node].child[c] < 0) {
            nodes[node].child[c] = int(nodes.size());
            nodes.emplace_back();
        }

        return nodes[node].child[c];
    }

    int build(const string &word, int delta) {
        int node = ROOT;

        for (char ch : word) {
            nodes[node].starting_with += delta;
            node = get_or_create_child(node, ch - MIN_CHAR);
        }

        nodes[node].starting_with += delta;
        return node;
    }

    int add(const string &word) {
        int node = build(word, +1);
        nodes[node].words_here++;
        return node;
    }

    int erase(const string &word) {
        int node = build(word, -1);
        nodes[node].words_here--;
        return node;
    }

    int find(const string &str) const {
        int node = ROOT;

        for (char ch : str) {
            node = nodes[node].child[ch - MIN_CHAR];

            if (node < 0)
                break;
        }

        return node;
    }

    // Given a string, how many words in the trie are prefixes of the string?
    int count_prefixes(const string &str, bool include_full) const {
        int node = ROOT, count = 0;

        for (char ch : str) {
            count += nodes[node].words_here;
            node = nodes[node].child[ch - MIN_CHAR];

            if (node < 0)
                break;
        }

        if (include_full && node >= 0)
            count += nodes[node].words_here;

        return count;
    }

    // Given a string, how many words in the trie start with the given string?
    int count_starting_with(const string &str, bool include_full) const {
        int node = find(str);

        if (node < 0)
            return 0;

        return nodes[node].starting_with - (include_full ? 0 : nodes[node].words_here);
    }
};


int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N;
    cin >> N;
    array_trie trie;
    vector<string> strings(N);

    for (int i = 0; i < N; i++) {
        cin >> strings[i];
        cout << trie.count_prefixes(strings[i], true) << ' ' << trie.count_starting_with(strings[i], true) << '\n';
        trie.add(strings[i]);

        if (i >= N / 2)
            trie.erase(strings[i - N / 2]);
    }
}
