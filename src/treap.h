#pragma once

#include <unordered_map>

// ==============================
// 模块 1：Treap（BST + 大根堆）
// ==============================
class Treap {
private:
    struct Node {
        int key;
        double priority;
        Node* left;
        Node* right;
        Node(int k, double p) : key(k), priority(p), left(nullptr), right(nullptr) {}
    };
    Node* root;

    Node* rotateRight(Node* x);
    Node* rotateLeft(Node* x);
    Node* insert(Node* t, Node* node);
    bool search(Node* t, int key, int& depth) const;
    void collectDepths(Node* t, int d, std::unordered_map<int, int>& out) const;
    int height(Node* t) const;
    void destroy(Node* t);

public:
    Treap();
    ~Treap();
    Treap(const Treap&) = delete;
    Treap& operator=(const Treap&) = delete;

    void insert(int key, double priority);
    bool search(int key, int& depth) const;
    std::unordered_map<int, int> allDepths() const;
    int height() const;
};
