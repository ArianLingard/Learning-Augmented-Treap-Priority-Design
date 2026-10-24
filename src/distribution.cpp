#include "distribution.h"

#include <algorithm>
#include <cmath>
#include <numeric>

using namespace std;

vector<double> genDistribution(const string& name, double param, int n) {
    vector<double> p(n);
    if (name == "Uniform") {
        for (auto& x : p) x = 1.0 / n;
    }
    else if (name == "Zipf") {
        double alpha = param, sum = 0.0;
        for (int i = 1; i <= n; i++) sum += 1.0 / pow((double)i, alpha);
        for (int i = 0; i < n; i++) p[i] = 1.0 / pow((double)(i + 1), alpha) / sum;
    }
    else if (name == "Adversarial") {
        // 对抗分布：p_i = 2(n-i+1)/(n(n+1)), i=1..n
        for (int i = 0; i < n; i++) p[i] = 2.0 * (double)(n - i) / ((double)n * (n + 1));
    }
    else if (name == "Hotspot") {
        // 极端热点：p_1=0.999，其余均分
        p[0] = 0.999;
        for (int i = 1; i < n; i++) p[i] = (1.0 - 0.999) / (double)(n - 1);
    }
    return p;
}

vector<double> assignWeights(const vector<double>& p, int n, mt19937_64& rng) {
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), rng);
    vector<double> w(n);
    for (int i = 0; i < n; i++) w[perm[i]] = p[i];
    return w;
}

vector<int> genAccessSequence(const vector<double>& w, int m, mt19937_64& rng) {
    discrete_distribution<int> dist(w.begin(), w.end());
    vector<int> seq(m);
    for (auto& x : seq) x = dist(rng);
    return seq;
}
