#include "priority.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <utility>

using namespace std;

KeyInfo computeKeyInfo(const vector<double>& w, mt19937_64& rng) {
    int n = (int)w.size();
    vector<int> tie(n);
    iota(tie.begin(), tie.end(), 0);
    shuffle(tie.begin(), tie.end(), rng);       // 随机打破同权重并列
    vector<int> orderDesc(n), orderAsc(n);
    iota(orderDesc.begin(), orderDesc.end(), 0);
    iota(orderAsc.begin(), orderAsc.end(), 0);
    sort(orderDesc.begin(), orderDesc.end(), [&](int a, int b) {
        if (w[a] != w[b]) return w[a] > w[b];   // 降序（热->冷）
        return tie[a] < tie[b];
    });
    sort(orderAsc.begin(), orderAsc.end(), [&](int a, int b) {
        if (w[a] != w[b]) return w[a] < w[b];   // 升序（冷->热）
        return tie[a] < tie[b];
    });
    vector<int> rankDesc(n);
    for (int i = 0; i < n; i++) rankDesc[orderDesc[i]] = i + 1;

    double total = accumulate(w.begin(), w.end(), 0.0);
    vector<double> F(n, 0.0);
    double prefix = 0.0;
    for (int x : orderAsc) { F[x] = prefix / total; prefix += w[x]; }  // 严格更轻质量占比

    return { n, w, rankDesc, F };
}

double basePriority(const string& method, const KeyInfo& ki, int x) {
    double w = ki.w[x];
    if (method == "Chen") {
        return -floor(log2(log2(1.0 / w)));          // w==1 时 log2(0)=-inf -> +inf
    }
    if (method == "f1") return -log2(log2(1.0 / w) + EPS);
    if (method == "f2") {
        double K = log2(log2((double)ki.n));
        return K - (double)ki.rankDesc[x] * K / (double)ki.n;
    }
    if (method == "f3") {
        double p = (double)(ki.n - ki.rankDesc[x] + 1) / (double)ki.n;
        return -log2(log2(1.0 / p) + EPS);
    }
    if (method == "f4") {
        double F = ki.F[x];
        return -log2(log2(1.0 / (F + 1.0 / ki.n)) + EPS);
    }
    return 0.0;
}

vector<double> buildPriorities(const string& method, const KeyInfo& ki, PerturbMode mode, mt19937_64& rng) {
    int n = ki.n;
    vector<double> pri(n);
    if (method == "Random") {
        uniform_real_distribution<double> u(0.0, 1.0);
        for (int i = 0; i < n; i++) pri[i] = u(rng);
        return pri;
    }
    vector<double> base(n);
    for (int i = 0; i < n; i++) base[i] = basePriority(method, ki, i);

    if (method == "Chen") {   // Chen 自带层间距 1，扰动 U(0,1)
        uniform_real_distribution<double> u(0.0, 1.0);
        for (int i = 0; i < n; i++) pri[i] = base[i] + u(rng);
        return pri;
    }
    // f1/f2/f3/f4：Q(f)=floor(f)，扰动加在桶外
    if (mode == NONE) {
        for (int i = 0; i < n; i++) pri[i] = base[i];
        return pri;
    }
    if (mode == FIXED) {
        uniform_real_distribution<double> u(0.0, FIXED_PERTURB);
        for (int i = 0; i < n; i++) pri[i] = floor(base[i]) + u(rng);
        return pri;
    }
    // ADAPTIVE：delta_max,k = min(0.49, 0.49*gap_k)，gap_k 为当前桶与下一个不同桶的间距。
    // 只要 delta < 桶间距就不会跨桶（保留学习到的排序）；floor 分桶的 gap 恒为 1，故 delta_max=0.49。
    vector<int> bucket(n);
    unordered_map<int, int> cnt;
    for (int i = 0; i < n; i++) { bucket[i] = (int)floor(base[i]); cnt[bucket[i]]++; }
    vector<int> uniq;
    for (auto& kv : cnt) uniq.push_back(kv.first);
    sort(uniq.begin(), uniq.end());
    unordered_map<int, double> gap;
    for (size_t t = 0; t < uniq.size(); t++)
        gap[uniq[t]] = (t + 1 < uniq.size()) ? (double)(uniq[t + 1] - uniq[t]) : 1.0;
    for (int i = 0; i < n; i++) {
        double dmax = min(0.49, 0.49 * gap[bucket[i]]);
        uniform_real_distribution<double> u(0.0, dmax);
        pri[i] = (double)bucket[i] + u(rng);
    }
    return pri;
}

vector<double> buildPrioritiesPhi(const string& method, const KeyInfo& ki, double phi, mt19937_64& rng) {
    int n = ki.n;
    vector<double> pri(n);
    uniform_real_distribution<double> u(0.0, phi);
    for (int i = 0; i < n; i++) pri[i] = floor(basePriority(method, ki, i)) + u(rng);
    return pri;
}

static long long countInversions(vector<double>& a, vector<double>& tmp, int l, int r) {
    if (r - l <= 1) return 0;
    int m = (l + r) / 2;
    long long inv = countInversions(a, tmp, l, m) + countInversions(a, tmp, m, r);
    int i = l, j = m, k = l;
    while (i < m && j < r) {
        if (a[i] <= a[j]) tmp[k++] = a[i++];
        else { tmp[k++] = a[j++]; inv += (long long)(m - i); }
    }
    while (i < m) tmp[k++] = a[i++];
    while (j < r) tmp[k++] = a[j++];
    for (int t = l; t < r; t++) a[t] = tmp[t];
    return inv;
}

double kendallTau(const vector<double>& a, const vector<double>& b) {
    int n = (int)a.size();
    vector<pair<double, double>> pairs(n);
    for (int i = 0; i < n; i++) pairs[i] = { a[i], b[i] };
    sort(pairs.begin(), pairs.end());
    vector<double> bb(n), tmp(n);
    for (int i = 0; i < n; i++) bb[i] = pairs[i].second;
    long long inv = countInversions(bb, tmp, 0, n);
    long long total = (long long)n * (n - 1) / 2;
    return 1.0 - 2.0 * (double)inv / (double)total;
}

double crossBucketInversionRatio(const string& method, const KeyInfo& ki, const vector<double>& pri) {
    int n = ki.n;
    vector<int> bucket(n);
    for (int i = 0; i < n; i++) bucket[i] = (int)floor(basePriority(method, ki, i));
    vector<int> idx(n);
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int a, int b) {
        if (bucket[a] != bucket[b]) return bucket[a] < bucket[b];
        return a < b;
    });
    vector<double> sortedPri(n), tmp(n);
    for (int i = 0; i < n; i++) sortedPri[i] = pri[idx[i]];
    vector<double> a = sortedPri;
    long long totalInv = countInversions(a, tmp, 0, n);
    long long withinInv = 0;
    for (int i = 0; i < n; ) {
        int j = i;
        while (j < n && bucket[idx[j]] == bucket[idx[i]]) j++;
        int sz = j - i;
        if (sz > 1) {
            vector<double> sub(sz), stmp(sz);
            for (int k = 0; k < sz; k++) sub[k] = pri[idx[i + k]];
            withinInv += countInversions(sub, stmp, 0, sz);
        }
        i = j;
    }
    long long crossInv = totalInv - withinInv;
    unordered_map<int, int> cnt;
    for (int i = 0; i < n; i++) cnt[bucket[i]]++;
    long long same = 0;
    for (auto& kv : cnt) { long long k = kv.second; same += k * (k - 1) / 2; }
    long long totalCross = (long long)n * (n - 1) / 2 - same;
    if (totalCross == 0) return 0.0;
    return (double)crossInv / (double)totalCross;
}
