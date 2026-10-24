#include "treap.h"

#include <algorithm>

using namespace std;

Treap::Treap() : root(nullptr) {}
Treap::~Treap() { destroy(root); }

Treap::Node* Treap::rotateRight(Node* x) {
    Node* y = x->left;
    x->left = y->right;
    y->right = x;
    return y;
}

Treap::Node* Treap::rotateLeft(Node* x) {
    Node* y = x->right;
    x->right = y->left;
    y->left = x;
    return y;
}

Treap::Node* Treap::insert(Node* t, Node* node) {
    if (!t) return node;
    if (node->key < t->key) {
        t->left = insert(t->left, node);
        if (t->left->priority > t->priority) t = rotateRight(t);
    }
    else if (node->key > t->key) {
        t->right = insert(t->right, node);
        if (t->right->priority > t->priority) t = rotateLeft(t);
    }
    return t;
}

bool Treap::search(Node* t, int key, int& depth) const {
    depth = 0;
    while (t) {
        depth++;
        if (key == t->key) return true;
        t = (key < t->key) ? t->left : t->right;
    }
    return false;
}

void Treap::collectDepths(Node* t, int d, unordered_map<int, int>& out) const {
    if (!t) return;
    out[t->key] = d;
    collectDepths(t->left, d + 1, out);
    collectDepths(t->right, d + 1, out);
}

int Treap::height(Node* t) const {
    if (!t) return 0;
    return 1 + max(height(t->left), height(t->right));
}

void Treap::destroy(Node* t) {
    if (!t) return;
    destroy(t->left);
    destroy(t->right);
    delete t;
}

void Treap::insert(int key, double priority) { root = insert(root, new Node(key, priority)); }
bool Treap::search(int key, int& depth) const { return search(root, key, depth); }
unordered_map<int, int> Treap::allDepths() const {
    unordered_map<int, int> m;
    collectDepths(root, 1, m);
    return m;
}
int Treap::height() const { return height(root); }
