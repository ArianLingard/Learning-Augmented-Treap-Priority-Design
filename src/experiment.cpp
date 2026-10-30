#include "experiment.h"

#include "distribution.h"
#include "oracle.h"
#include "priority.h"
#include "treap.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace std;

namespace {

// =====================================================================
// 全局统一设置（见 实验设计.md "全局统一设置"）
//   - 键序随机打乱；每个配置 30 次随机种子，报 均值±标准差
//   - 分桶 Q(f)=floor(f)，固定扰动 delta ~ U(0,0.25)
//   - 学习增强方法同组共用同一个预测频率 w_x
// =====================================================================

const int M     = 100000;   // 访问序列长度（E1/E2/E4/E5/E6）
const int SEEDS = 30;       // 每个配置的随机种子数

// 统计工具
struct Stat { double mean; double sd; };
Stat stat(const vector<double>& v) {
    double m = accumulate(v.begin(), v.end(), 0.0) / (double)v.size();
    double var = 0.0;
    for (double x : v) var += (x - m) * (x - m);
    var /= (double)(v.size() - 1);   // 样本标准差
    return { m, sqrt(var) };
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

// 最优静态 BST 代价：Knuth O(n²) 动态规划。
// 输入每个键的访问频率 f[k]（k=0..n-1，已按键值升序），返回
//   OPT = min_{BST} Σ_k f[k] · depth(k)，
// 其中 depth(k) 为节点访问次数（根=1），与 Treap::search 统计的比较次数口径一致。
// 因此 ρ = 实测比较次数 / OPT 恒 >= 1，是真正意义上的静态最优差距。
double optimalBSTCost(const vector<long long>& f) {
    int n = (int)f.size();
    if (n == 0) return 0.0;
    vector<double> w(n + 1);                       // 1-indexed 权重
    for (int k = 0; k < n; k++) w[k + 1] = (double)f[k];
    vector<double> sum(n + 1, 0.0);                // 前缀和 W(i,j)=sum[j]-sum[i-1]
    for (int i = 1; i <= n; i++) sum[i] = sum[i - 1] + w[i];

    vector<vector<double>> e(n + 2, vector<double>(n + 1, 0.0));  // 空区间 e[i][i-1]=0
    vector<vector<int>> root(n + 2, vector<int>(n + 1, 0));
    for (int i = 1; i <= n; i++) { e[i][i] = w[i]; root[i][i] = i; }

    for (int len = 2; len <= n; len++) {
        for (int i = 1; i + len - 1 <= n; i++) {
            int j = i + len - 1;
            double W = sum[j] - sum[i - 1];
            // Knuth 单调性：root[i][j-1] <= root[i][j] <= root[i+1][j]
            int lo = root[i][j - 1];
            int hi = root[i + 1][j];
            double best = numeric_limits<double>::infinity();
            int bestRoot = lo;
            for (int r = lo; r <= hi; r++) {
                double t = e[i][r - 1] + e[r + 1][j];
                if (t < best) { best = t; bestRoot = r; }
            }
            e[i][j] = best + W;
            root[i][j] = bestRoot;
        }
    }
    return e[1][n];
}

} // namespace

// E1：常规分布下基础性能对比
void runE1(ofstream& out) {
    vector<pair<string, double>> dists = {
        {"Zipf", 1.0}, {"Zipf", 1.25}, {"Zipf", 2.0}, {"Adversarial", 0.0}, {"Uniform", 0.0}
    };
    vector<int> ns = { 1000, 5000 };
    vector<string> methods = { "Random", "Chen", "f1", "f2", "f3", "f4" };
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
            vector<vector<double>> rhos(methods.size());
            for (int seed = 0; seed < SEEDS; seed++) {
                mt19937_64 rng(seed + 200);
                vector<double> w = assignWeights(pn, n, rng);
                KeyInfo ki = computeKeyInfo(w, rng);
                vector<int> seq = genAccessSequence(w, M, rng);
                vector<int> keysOrder(n); iota(keysOrder.begin(), keysOrder.end(), 0);
                shuffle(keysOrder.begin(), keysOrder.end(), rng);
                // 真实访问频率 -> 最优静态 BST 代价（所有方法共用同一 seq，故 OPT 每 seed 只算一次）
                vector<long long> f(n, 0);
                for (int q : seq) f[q]++;
                double opt = optimalBSTCost(f);
                for (size_t mi = 0; mi < methods.size(); mi++) {
                    vector<double> pri = buildPriorities(methods[mi], ki, FIXED, rng);
                    auto res = evaluate(n, keysOrder, pri, seq);
                    rhos[mi].push_back((double)res.comparisons / opt);
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
