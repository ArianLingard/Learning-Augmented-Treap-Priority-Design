#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <string>
#include <fstream>
#include <chrono>
#include <map>
using namespace std;

// =====================================================================
// 全局统一设置（见 实验设计.md "全局统一设置"）
//   - 键序随机打乱；每个配置 30 次随机种子，报 均值±标准差
//   - 分桶 Q(f)=floor(f)，固定扰动 delta ~ U(0,0.25)
//   - 学习增强方法同组共用同一个预测频率 w_x
// =====================================================================

const int   M             = 100000;    // 访问序列长度（E1/E2/E4/E5/E6）
const int   SEEDS         = 30;        // 每个配置的随机种子数
const double EPS          = 1e-9;      // 平滑双对数的正则常数（f1/f3/f4）
const double FIXED_PERTURB= 0.25;      // 固定扰动上界 U(0,0.25)

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

    Node* rotateRight(Node* x) {
        Node* y = x->left;
        x->left = y->right;
        y->right = x;
        return y;
    }

    Node* rotateLeft(Node* x) {
        Node* y = x->right;
        x->right = y->left;
        y->left = x;
        return y;
    }

    Node* insert(Node* t, Node* node) {
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

    bool search(Node* t, int key, int& depth) const {
        depth = 0;
        while (t) {
            depth++;
            if (key == t->key) return true;
            t = (key < t->key) ? t->left : t->right;
        }
        return false;
    }

    void collectDepths(Node* t, int d, unordered_map<int, int>& out) const {
        if (!t) return;
        out[t->key] = d;
        collectDepths(t->left, d + 1, out);
        collectDepths(t->right, d + 1, out);
    }

    int height(Node* t) const {
        if (!t) return 0;
        return 1 + max(height(t->left), height(t->right));
    }

    void destroy(Node* t) {
        if (!t) return;
        destroy(t->left);
        destroy(t->right);
        delete t;
    }

public:
    Treap() : root(nullptr) {}
    ~Treap() { destroy(root); }
    Treap(const Treap&) = delete;
    Treap& operator=(const Treap&) = delete;

    void insert(int key, double priority) { root = insert(root, new Node(key, priority)); }
    bool search(int key, int& depth) const { return search(root, key, depth); }
    unordered_map<int, int> allDepths() const {
        unordered_map<int, int> m;
        collectDepths(root, 1, m);
        return m;
    }
    int height() const { return height(root); }
};

// ==============================
// 统计工具
// ==============================
struct Stat { double mean; double sd; };
Stat stat(const vector<double>& v) {
    double m = accumulate(v.begin(), v.end(), 0.0) / (double)v.size();
    double var = 0.0;
    for (double x : v) var += (x - m) * (x - m);
    var /= (double)(v.size() - 1);   // 样本标准差
    return { m, sqrt(var) };
}

// ====================================================
// 模块 2：访问分布（键序随机打乱 -> 权重按随机排列赋给键）
// ====================================================
// 返回长度为 n 的概率向量 p（按 rank 从热到冷排列：p[0] 最热）
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

// 键序随机打乱：把 rank 概率 p 随机赋给键 0..n-1
vector<double> assignWeights(const vector<double>& p, int n, mt19937_64& rng) {
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), rng);
    vector<double> w(n);
    for (int i = 0; i < n; i++) w[perm[i]] = p[i];
    return w;
}

// 依据权重向量 w 采样长度为 m 的访问序列（同一组实验内所有方法共用）
vector<int> genAccessSequence(const vector<double>& w, int m, mt19937_64& rng) {
    discrete_distribution<int> dist(w.begin(), w.end());
    vector<int> seq(m);
    for (auto& x : seq) x = dist(rng);
    return seq;
}

// ====================================================
// 模块 3：预测器 Oracle（完美频率 / 量纲漂移 / KL 噪声）
// ====================================================
// 完美频率：w_x = p_x（直接使用）
// 缩放/幂律（E2）：w_x = c*p_x 或 (p_x)^beta（严格单调变换），
// 双对数定义域为 (0,1)：把可能越界的预测权重夹回 [1e-300, 1-1e-12]
static inline double clampWeight(double x) { return min(max(x, 1e-300), 1.0 - 1e-12); }

double klDivergence(const vector<double>& P, const vector<double>& Q) {
    double kl = 0.0;
    for (size_t i = 0; i < P.size(); i++) if (P[i] > 0) kl += P[i] * log2(P[i] / Q[i]);
    return kl;
}

// 逐键乘性对数正态噪声：Q ∝ P·exp(σ z)，z~N(0,1)，二分搜索 σ 使 D_KL(P||Q) 逼近 target。
// 乘性误差既扰动幅度也扰动相对排序（比纯温度缩放更贴近真实预测器）。
vector<double> makeNoisy(const vector<double>& P, double targetKL, mt19937_64& rng) {
    int n = (int)P.size();
    normal_distribution<double> nd(0.0, 1.0);
    vector<double> z(n);
    for (auto& x : z) x = nd(rng);
    auto qForSigma = [&](double sigma) {
        vector<double> q(n);
        double s = 0.0;
        for (int i = 0; i < n; i++) { q[i] = P[i] * exp(sigma * z[i]); s += q[i]; }
        for (auto& x : q) x /= s;
        return q;
    };
    double lo = 0.0, hi = 20.0;   // KL 随 sigma 单调递增：sigma=0 完美(0)，sigma 大则噪声主导
    for (int it = 0; it < 200; it++) {
        double mid = (lo + hi) / 2.0;
        double kl = klDivergence(P, qForSigma(mid));
        if (kl < targetKL) lo = mid; else hi = mid;
    }
    return qForSigma((lo + hi) / 2.0);
}

// ====================================================
// 模块 4：priority 映射函数
// ====================================================
struct KeyInfo {
    int n;
    vector<double> w;        // 每个键的预测权重
    vector<int>    rankDesc; // 降序排名（1 = 最热）
    vector<double> F;        // 归一化 CDF（严格更轻键的质量占比，尺度不变）
};

// 由预测权重计算降序排名与平滑秩 CDF（同权重用随机排列打破，避免均匀分布退化为链）
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

// 基础映射 h（分桶/扰动之前的确定性部分）
//   Chen : -floor(log2(log2(1/w)))                 （w->1 时 -> +inf）
//   f1   : -log2(log2(1/w) + eps)                  （平滑双对数，右端有界）
//   f2   : K - R*K/n, K=log2(log2 n)               （线性秩映射，丢弃幅度）
//   f3   : -log2(log2(1/p) + eps), p=(n-R+1)/n     （秩双对数）
//   f4   : -log2(log2(1/(F+1/n)) + eps)            （CDF 平滑秩双对数）
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

enum PerturbMode { NONE, FIXED, ADAPTIVE };

// 由基础映射 + 分桶 + 扰动构造最终 priority 向量
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

// 固定幅度扰动（E6 扫描用）：Q(f)=floor(f)，δ~U(0,phi)。phi 可跨过桶宽 1 以暴露跨桶逆序。
vector<double> buildPrioritiesPhi(const string& method, const KeyInfo& ki, double phi, mt19937_64& rng) {
    int n = ki.n;
    vector<double> pri(n);
    uniform_real_distribution<double> u(0.0, phi);
    for (int i = 0; i < n; i++) pri[i] = floor(basePriority(method, ki, i)) + u(rng);
    return pri;
}

// ====================================================
// 模块 5：评估与统计
// ====================================================

// Kendall tau：O(n log n)，tau=1 表示两个 priority 序完全一致（拓扑不变）
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

// 跨桶逆序率：只统计“不同桶”之间的 key 对里，优先级序与权重序相反的比例。
// 桶 = floor(basePriority)。它直接度量扰动把学习到的排序破坏了多少（扰动幅度 >= 桶间距时才会 >0）。
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

struct EvalResult { long long comparisons; double avgDepth; int maxDepth; };
EvalResult evaluate(int n, const vector<int>& keysOrder, const vector<double>& pri, const vector<int>& seq) {
    Treap t;
    for (int k : keysOrder) t.insert(k, pri[k]);
    long long comps = 0;
    for (int q : seq) { int d; t.search(q, d); comps += d; }
    auto depths = t.allDepths();
    double avg = 0.0;
    for (auto& kv : depths) avg += kv.second;
    avg /= (double)n;
    return { comps, avg, t.height() };
}

double entropyBits(const vector<double>& p) {
    double h = 0.0;
    for (double x : p) if (x > 0) h -= x * log2(x);
    return h;
}

// ====================================================
// 模块 6：实验设计（E1 - E7）
// ====================================================

// E1：常规分布下基础性能对比
void runE1(ofstream& out) {
    vector<pair<string, double>> dists = {
        {"Zipf", 1.0}, {"Zipf", 1.25}, {"Zipf", 2.0}, {"Adversarial", 0.0}, {"Uniform", 0.0}
    };
    vector<int> ns = { 1000, 5000 };
    vector<string> methods = { "Random", "Chen", "f1", "f3", "f4" };
    out << "distribution,param,n,method,mean_comparisons,std_comparisons,mean_avg_depth,std_avg_depth,mean_max_depth,std_max_depth\n";
    out << setprecision(6);
    for (auto& d : dists) {
        for (int n : ns) {
            vector<vector<double>> comps(methods.size()), avgd(methods.size()), maxd(methods.size());
            vector<double> p = genDistribution(d.first, d.second, n);
            for (int seed = 0; seed < SEEDS; seed++) {
                mt19937_64 rng(seed + 1);
                vector<double> w = assignWeights(p, n, rng);
                KeyInfo ki = computeKeyInfo(w, rng);
                vector<int> seq = genAccessSequence(w, M, rng);
                vector<int> keysOrder(n); iota(keysOrder.begin(), keysOrder.end(), 0);
                shuffle(keysOrder.begin(), keysOrder.end(), rng);
                for (size_t mi = 0; mi < methods.size(); mi++) {
                    vector<double> pri = buildPriorities(methods[mi], ki, FIXED, rng);
                    auto res = evaluate(n, keysOrder, pri, seq);
                    comps[mi].push_back((double)res.comparisons);
                    avgd[mi].push_back(res.avgDepth);
                    maxd[mi].push_back((double)res.maxDepth);
                }
            }
            for (size_t mi = 0; mi < methods.size(); mi++) {
                auto cs = stat(comps[mi]), as = stat(avgd[mi]), ms = stat(maxd[mi]);
                string param = (d.first == "Zipf") ? to_string(d.second) : "";
                out << d.first << "," << param << "," << n << "," << methods[mi] << ","
                    << cs.mean << "," << cs.sd << "," << as.mean << "," << as.sd << ","
                    << ms.mean << "," << ms.sd << "\n";
            }
        }
    }
}

// E2：绝对量纲漂移的鲁棒性（缩放 / 幂律）
void runE2(ofstream& out) {
    int n = 5000;
    vector<pair<string, double>> transforms = { {"scale", 0.1}, {"scale", 10.0}, {"power", 0.5}, {"power", 2.0} };
    vector<string> methods = { "Chen", "f1", "f3", "f4" };
    out << "transform,param,method,mean_comparisons,std_comparisons,mean_kendall_tau,std_kendall_tau\n";
    out << setprecision(6);
    vector<double> p = genDistribution("Zipf", 1.0, n);
    for (auto& tr : transforms) {
        for (auto& method : methods) {
            vector<double> comps, taus;
            bool isChen = (method == "Chen");
            for (int seed = 0; seed < SEEDS; seed++) {
                mt19937_64 rng(seed + 100);
                vector<double> w = assignWeights(p, n, rng);   // 真实分布仍是 Zipf(1)
                KeyInfo kiPerfect = computeKeyInfo(w, rng);
                vector<int> seq = genAccessSequence(w, M, rng);
                vector<int> keysOrder(n); iota(keysOrder.begin(), keysOrder.end(), 0);
                shuffle(keysOrder.begin(), keysOrder.end(), rng);
                // 变换作用于“完美频率”（每个键的真实权重 w），保持键序对应
                vector<double> wTr(n);
                for (int k = 0; k < n; k++)
                    wTr[k] = clampWeight((tr.first == "power") ? pow(w[k], tr.second) : tr.second * w[k]);
                KeyInfo kiTr = computeKeyInfo(wTr, rng);
                // 固定扰动（基线/变换共用同一 delta，隔离变换本身的拓扑影响）
                uniform_real_distribution<double> u(0.0, isChen ? 1.0 : FIXED_PERTURB);
                vector<double> delta(n);
                for (int k = 0; k < n; k++) delta[k] = u(rng);
                // 完整 priority（base + 固定 delta），用于 Kendall tau 与建树评估
                vector<double> priP(n), priT(n);
                for (int k = 0; k < n; k++) {
                    double bp = basePriority(method, kiPerfect, k);
                    double bt = basePriority(method, kiTr, k);
                    priP[k] = (isChen ? bp : floor(bp)) + delta[k];
                    priT[k] = (isChen ? bt : floor(bt)) + delta[k];
                }
                taus.push_back(kendallTau(priP, priT));
                auto res = evaluate(n, keysOrder, priT, seq);
                comps.push_back((double)res.comparisons);
            }
            auto cs = stat(comps), ts = stat(taus);
            out << tr.first << "," << tr.second << "," << method << ","
                << cs.mean << "," << cs.sd << "," << ts.mean << "," << ts.sd << "\n";
        }
    }
}

// E3：右端奇点与数值溢出压力测试
void runE3(ofstream& out) {
    int n = 1000;
    vector<double> testWs = { 0.9999, 0.99999, 0.999999, 1.0 };
    vector<string> methods = { "Chen", "f1", "f3", "f4" };
    out << "w,method,max_priority,is_inf,is_nan,build_time_us\n";
    for (double w0 : testWs) {
        for (auto& method : methods) {
            // 构造极端热点：预测权重 键0 = w0，其余均分 (1-0.999)/(n-1)
            vector<double> w(n, (1.0 - 0.999) / (double)(n - 1));
            w[0] = w0;
            mt19937_64 rng(12345);   // 仅用于同权重随机并列打破，max priority 与其无关
            KeyInfo ki = computeKeyInfo(w, rng);
            double maxp = -INFINITY; bool inf = false, nan = false;
            for (int i = 0; i < n; i++) {
                double b = basePriority(method, ki, i);
                if (isinf(b)) inf = true;
                if (isnan(b)) nan = true;
                if (b > maxp) maxp = b;
            }
            // 建树耗时（10 种子平均）
            double tsum = 0.0;
            for (int seed = 0; seed < SEEDS; seed++) {
                mt19937_64 rng(seed + 1000);
                vector<double> pri = buildPriorities(method, ki, FIXED, rng);
                vector<int> keysOrder(n); iota(keysOrder.begin(), keysOrder.end(), 0);
                shuffle(keysOrder.begin(), keysOrder.end(), rng);
                auto t0 = chrono::high_resolution_clock::now();
                Treap tr;
                for (int k : keysOrder) tr.insert(k, pri[k]);
                auto t1 = chrono::high_resolution_clock::now();
                tsum += (double)chrono::duration_cast<chrono::microseconds>(t1 - t0).count();
            }
            out << setprecision(9) << w0 << "," << method << "," << maxp << ","
                << (inf ? 1 : 0) << "," << (nan ? 1 : 0) << "," << setprecision(3) << (tsum / SEEDS) << "\n";
        }
    }
}

// E4：理论静态最优代价分析
void runE4(ofstream& out) {
    vector<pair<string, double>> dists = {
        {"Zipf", 1.0}, {"Zipf", 1.25}, {"Zipf", 2.0}, {"Adversarial", 0.0}, {"Uniform", 0.0}
    };
    vector<int> ns = { 1000, 5000 };
    vector<string> methods = { "Chen", "f1", "f2", "f3", "f4" };
    out << "distribution,param,n,method,mean_rho,std_rho\n";
    out << setprecision(6);
    for (auto& d : dists) {
        for (int n : ns) {
            vector<double> pn = genDistribution(d.first, d.second, n);
            double denom = M * entropyBits(pn);   // 香农熵下界：sum f_x log2(m/f_x)
            vector<vector<double>> rhos(methods.size());
            for (int seed = 0; seed < SEEDS; seed++) {
                mt19937_64 rng(seed + 200);
                vector<double> w = assignWeights(pn, n, rng);
                KeyInfo ki = computeKeyInfo(w, rng);
                vector<int> seq = genAccessSequence(w, M, rng);
                vector<int> keysOrder(n); iota(keysOrder.begin(), keysOrder.end(), 0);
                shuffle(keysOrder.begin(), keysOrder.end(), rng);
                for (size_t mi = 0; mi < methods.size(); mi++) {
                    vector<double> pri = buildPriorities(methods[mi], ki, FIXED, rng);
                    auto res = evaluate(n, keysOrder, pri, seq);
                    rhos[mi].push_back((double)res.comparisons / denom);
                }
            }
            for (size_t mi = 0; mi < methods.size(); mi++) {
                auto rs = stat(rhos[mi]);
                string param = (d.first == "Zipf") ? to_string(d.second) : "";
                out << d.first << "," << param << "," << n << "," << methods[mi] << ","
                    << rs.mean << "," << rs.sd << "\n";
            }
        }
    }
}

// E5：非完美预测器（KL 噪声）鲁棒性
void runE5(ofstream& out) {
    int n = 5000;
    vector<double> kls = { 0.1, 0.3, 0.5, 0.7 };
    vector<string> methods = { "Random", "Chen", "f1", "f3", "f4" };
    out << "kl,method,mean_comparisons,std_comparisons,mean_relative_degradation,std_relative_degradation\n";
    out << setprecision(6);
    vector<double> p = genDistribution("Zipf", 1.0, n);
    // 成对比较：每个 seed 的访问序列/键序固定，仅噪声预测分布 Q 变化
    vector<vector<double>> base(methods.size(), vector<double>(SEEDS, 0.0));
    vector<vector<vector<double>>> comps(kls.size(),
        vector<vector<double>>(methods.size(), vector<double>(SEEDS, 0.0)));
    for (int seed = 0; seed < SEEDS; seed++) {
        mt19937_64 rng(seed + 300);
        vector<double> w = assignWeights(p, n, rng);
        KeyInfo kiPerfect = computeKeyInfo(w, rng);
        vector<int> seq = genAccessSequence(w, M, rng);
        vector<int> keysOrder(n); iota(keysOrder.begin(), keysOrder.end(), 0);
        shuffle(keysOrder.begin(), keysOrder.end(), rng);
        // KL=0 完美预测基线
        for (size_t mi = 0; mi < methods.size(); mi++) {
            vector<double> pri = buildPriorities(methods[mi], kiPerfect, FIXED, rng);
            auto res = evaluate(n, keysOrder, pri, seq);
            base[mi][seed] = (double)res.comparisons;
        }
        // 各 KL 噪声预测
        for (size_t ki = 0; ki < kls.size(); ki++) {
            mt19937_64 noiseRng((uint64_t)seed * 1000 + (uint64_t)(kls[ki] * 100));
            vector<double> q = makeNoisy(p, kls[ki], noiseRng);
            vector<double> wPred(n);
            for (int k = 0; k < n; k++) wPred[k] = q[kiPerfect.rankDesc[k] - 1];
            KeyInfo kiQ = computeKeyInfo(wPred, rng);
            for (size_t mi = 0; mi < methods.size(); mi++) {
                vector<double> pri = buildPriorities(methods[mi], kiQ, FIXED, rng);
                auto res = evaluate(n, keysOrder, pri, seq);
                comps[ki][mi][seed] = (double)res.comparisons;
            }
        }
    }
    for (size_t ki = 0; ki < kls.size(); ki++) {
        for (size_t mi = 0; mi < methods.size(); mi++) {
            vector<double> csV(SEEDS), dgV(SEEDS);
            for (int seed = 0; seed < SEEDS; seed++) {
                csV[seed] = comps[ki][mi][seed];
                dgV[seed] = comps[ki][mi][seed] / base[mi][seed];
            }
            auto cs = stat(csV), dg = stat(dgV);
            out << kls[ki] << "," << methods[mi] << "," << cs.mean << "," << cs.sd << ","
                << dg.mean << "," << dg.sd << "\n";
        }
    }
}

// E6：自适应扰动 vs 固定幅度扫描（跨桶逆序率 + 查询代价）
void runE6(ofstream& out) {
    int n = 50000;
    vector<string> methods = { "f3", "f4" };
    vector<pair<string, double>> fixedPhis = {
        {"0.01", 0.01}, {"0.25", 0.25}, {"0.49", 0.49},
        {"0.75", 0.75}, {"0.99", 0.99}, {"1.25", 1.25}, {"1.75", 1.75}, {"2.0", 2.0}
    };
    struct Cfg { string method, label; bool adaptive; double phi; };
    vector<Cfg> configs;
    for (auto& m : methods) {
        for (auto& fp : fixedPhis) configs.push_back({ m, fp.first, false, fp.second });
        configs.push_back({ m, "adaptive", true, 0.0 });
    }
    out << "method,perturb,mean_inversion_ratio,std_inversion_ratio,mean_comparisons,std_comparisons,mean_max_depth,std_max_depth\n";
    out << setprecision(6);
    vector<double> p = genDistribution("Zipf", 2.0, n);
    int C = (int)configs.size();
    vector<vector<double>> invs(C), comps(C), maxds(C);
    for (int seed = 0; seed < SEEDS; seed++) {
        mt19937_64 rng(seed + 400);
        vector<double> w = assignWeights(p, n, rng);
        KeyInfo ki = computeKeyInfo(w, rng);
        vector<int> seq = genAccessSequence(w, M, rng);
        vector<int> keysOrder(n); iota(keysOrder.begin(), keysOrder.end(), 0);
        shuffle(keysOrder.begin(), keysOrder.end(), rng);
        for (int ci = 0; ci < C; ci++) {
            vector<double> pri = configs[ci].adaptive
                ? buildPriorities(configs[ci].method, ki, ADAPTIVE, rng)
                : buildPrioritiesPhi(configs[ci].method, ki, configs[ci].phi, rng);
            invs[ci].push_back(crossBucketInversionRatio(configs[ci].method, ki, pri));
            auto res = evaluate(n, keysOrder, pri, seq);
            comps[ci].push_back((double)res.comparisons);
            maxds[ci].push_back((double)res.maxDepth);
        }
    }
    for (int ci = 0; ci < C; ci++) {
        auto iv = stat(invs[ci]), cp = stat(comps[ci]), md = stat(maxds[ci]);
        out << configs[ci].method << "," << configs[ci].label << ","
            << iv.mean << "," << iv.sd << "," << cp.mean << "," << cp.sd << ","
            << md.mean << "," << md.sd << "\n";
    }
}

// E7：分桶 + 扰动必要性消融（对抗键序：key 0 最热，插入序 = 键序，去掉随机化以暴露退化）
void runE7(ofstream& out) {
    int n = 5000;
    const int M7 = 20000;   // 查询序列长度（消融实验：代价比值与 m 无关，取较小 m 控制运行时间）
    vector<string> settings = { "no_bucket_no_perturb", "bucket_no_perturb", "bucket_perturb", "bucket_adaptive" };
    out << "setting,mean_height,std_height,mean_comparisons,std_comparisons\n";
    out << setprecision(6);
    vector<double> p = genDistribution("Zipf", 1.0, n);
    vector<vector<double>> hts(settings.size()), comps(settings.size());
    for (int seed = 0; seed < SEEDS; seed++) {
        mt19937_64 rng(seed + 500);
        // 对抗键序：key i 直接拿 rank i 的权重（key 0 最热）；插入序 = 键序 0..n-1。
        vector<double> w = p;   // p 已按 rank 降序
        KeyInfo ki = computeKeyInfo(w, rng);
        vector<int> seq = genAccessSequence(w, M7, rng);
        vector<int> keysOrder(n); iota(keysOrder.begin(), keysOrder.end(), 0);
        vector<PerturbMode> modes = { NONE, NONE, FIXED, ADAPTIVE };
        for (size_t si = 0; si < settings.size(); si++) {
            vector<double> pri(n);
            if (settings[si] == "no_bucket_no_perturb") {
                for (int k = 0; k < n; k++) pri[k] = basePriority("f3", ki, k);   // 连续、不分桶不扰动
            } else {
                pri = buildPriorities("f3", ki, modes[si], rng);
            }
            Treap t;
            for (int k : keysOrder) t.insert(k, pri[k]);
            hts[si].push_back((double)t.height());
            long long cc = 0;
            for (int q : seq) { int d; t.search(q, d); cc += d; }
            comps[si].push_back((double)cc);
        }
    }
    for (size_t si = 0; si < settings.size(); si++) {
        auto hs = stat(hts[si]), cs = stat(comps[si]);
        out << settings[si] << "," << hs.mean << "," << hs.sd << "," << cs.mean << "," << cs.sd << "\n";
    }
}

// 简单计时器：返回从起点 t0 到现在的秒数（高精度时钟，微秒级）
static double elapsedSeconds(chrono::high_resolution_clock::time_point t0) {
    return chrono::duration_cast<chrono::microseconds>(
        chrono::high_resolution_clock::now() - t0).count() / 1e6;
}

int main(int argc, char** argv) {
    vector<string> run;
    for (int i = 1; i < argc; i++) run.push_back(argv[i]);
    bool all = run.empty();
    auto shouldRun = [&](const string& id) {
        return all || find(run.begin(), run.end(), id) != run.end();
    };

    cout << "Learning-Augmented Treap experiments (E1-E7)\n";
    auto tAll = chrono::high_resolution_clock::now();

    if (shouldRun("e1")) {
        auto t0 = chrono::high_resolution_clock::now();
        ofstream out("data/e1_basic.csv");
        runE1(out); out.close();
        cout << "  wrote data/e1_basic.csv  [" << fixed << setprecision(2) << elapsedSeconds(t0) << " s]\n";
    }
    if (shouldRun("e2")) {
        auto t0 = chrono::high_resolution_clock::now();
        ofstream out("data/e2_drift.csv");
        runE2(out); out.close();
        cout << "  wrote data/e2_drift.csv  [" << fixed << setprecision(2) << elapsedSeconds(t0) << " s]\n";
    }
    if (shouldRun("e3")) {
        auto t0 = chrono::high_resolution_clock::now();
        ofstream out("data/e3_singularity.csv");
        runE3(out); out.close();
        cout << "  wrote data/e3_singularity.csv  [" << fixed << setprecision(2) << elapsedSeconds(t0) << " s]\n";
    }
    if (shouldRun("e4")) {
        auto t0 = chrono::high_resolution_clock::now();
        ofstream out("data/e4_static_optimal.csv");
        runE4(out); out.close();
        cout << "  wrote data/e4_static_optimal.csv  [" << fixed << setprecision(2) << elapsedSeconds(t0) << " s]\n";
    }
    if (shouldRun("e5")) {
        auto t0 = chrono::high_resolution_clock::now();
        ofstream out("data/e5_kl_noise.csv");
        runE5(out); out.close();
        cout << "  wrote data/e5_kl_noise.csv  [" << fixed << setprecision(2) << elapsedSeconds(t0) << " s]\n";
    }
    if (shouldRun("e6")) {
        auto t0 = chrono::high_resolution_clock::now();
        ofstream out("data/e6_adaptive_perturb.csv");
        runE6(out); out.close();
        cout << "  wrote data/e6_adaptive_perturb.csv  [" << fixed << setprecision(2) << elapsedSeconds(t0) << " s]\n";
    }
    if (shouldRun("e7")) {
        auto t0 = chrono::high_resolution_clock::now();
        ofstream out("data/e7_bucketing.csv");
        runE7(out); out.close();
        cout << "  wrote data/e7_bucketing.csv  [" << fixed << setprecision(2) << elapsedSeconds(t0) << " s]\n";
    }
    cout << "done. total " << fixed << setprecision(2) << elapsedSeconds(tAll) << " s\n";
    return 0;
}
