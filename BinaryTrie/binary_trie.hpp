#ifndef TSUKASA_DIARY_BINARY_TRIE_HPP
#define TSUKASA_DIARY_BINARY_TRIE_HPP

#include <vector>
#include <algorithm>

template<typename T> class DoubleBinaryTrie {
	struct Node {
		size_t ch[4];
		Node() { std::fill(ch, ch+4, 0); }
	};
	
	std::vector<Node> node;
	std::vector<std::vector<T>> item;
	
	size_t expand(size_t i, size_t k) {
		node[i].ch[k] = node.size();
		node.push_back(Node());
		item.push_back(std::vector<T>());
		return node[i].ch[k];
	}
	
	size_t get_next(size_t& a, size_t& b) {
		size_t k = ((a & 1) << 1) | (b & 1);
		a >>= 1;
		b >>= 1;
		return k;
	}
	
public:
	DoubleBinaryTrie() : node(1, Node()), item(1, std::vector<T>()) {}
	
	void add(size_t a, size_t b, T x) {
		size_t i = 0;
		while (a > 0 or b > 0) {
			size_t k = get_next(a, b);
			if (node[i].ch[k] == 0) i = expand(i, k);
			else i = node[i].ch[k];
		}
		item[i].push_back(x);
	}
	
	std::vector<T> find(size_t a, size_t b) {
		size_t i = 0;
		while (a > 0 or b > 0) {
			size_t k = get_next(a, b);
			if (node[i].ch[k] == 0) return std::vector<T>();
			else i = node[i].ch[k];
		}
		return item[i];
	}
	
	bool erase(size_t a, size_t b) {
		size_t i = 0;
		while (a > 0 or b > 0) {
			size_t k = get_next(a, b);
			if (node[i].ch[k] == 0) return false;
			else i = node[i].ch[k];
		}
		item[i] = std::vector<T>();
		return true;
	}
};

#endif
