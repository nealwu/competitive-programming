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
        int words = 0;

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

    int add(const string &word) {
        int node = ROOT;

        for (char ch : word)
            node = get_or_create_child(node, ch - MIN_CHAR);

        nodes[node].words++;
        return node;
    }

    // Given a string, how many words in the trie are prefixes of the string?
    int count_prefixes(const string &str, bool include_full) {
        int node = ROOT, count = 0;

        for (char ch : str) {
            count += nodes[node].words;
            node = nodes[node].child[ch - MIN_CHAR];

            if (node < 0)
                break;
        }

        if (include_full && node >= 0)
            count += nodes[node].words;

        return count;
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

    for (int i = 0; i < N; i++) {
        string str;
        cin >> str;
        cout << trie.count_prefixes(str, true) << (i < N - 1 ? ' ' : '\n');
        trie.add(str);
    }
}
