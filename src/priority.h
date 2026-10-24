#pragma once

#include <random>
#include <string>
#include <vector>

// ====================================================
// 模块 4：priority 映射函数
// ====================================================

// Softplus 阻尼系数 D(x)=ln(1+e^{γx})-ln2（仅 f1 使用）。
//   下界 γ>ln(2e-1)≈1.49（文档记为 ≈1.55）保证 D(1)>1、消除右端奇点；
//   上界 γ·e^γ/(1+e^γ)<1/ln2（即 γ≲1.71）保证 f 随权重单调不减；取 γ=1.6。
constexpr double GAMMA         = 1.6;
constexpr double EPS           = 1e-9;   // 平滑双对数的正则常数（f3/f4）
constexpr double FIXED_PERTURB = 0.25;   // 固定扰动上界 U(0,0.25)

struct KeyInfo {
    int n;
    std::vector<double> w;        // 每个键的预测权重
    std::vector<int>    rankDesc; // 降序排名（1 = 最热）
    std::vector<double> F;        // 归一化 CDF（严格更轻键的质量占比，尺度不变）
};

// 由预测权重计算降序排名与平滑秩 CDF（同权重用随机排列打破，避免均匀分布退化为链）
KeyInfo computeKeyInfo(const std::vector<double>& w, std::mt19937_64& rng);

// 基础映射 h（分桶/扰动之前的确定性部分），D(x)=ln(1+e^{γx})-ln2 为 Softplus 平滑阻尼项
//   Chen : -floor(log2(log2(1/w)))                 （w->1 时 -> +inf）
//   f1   : -log2(log2(1/w) + D(w))                 （平滑双对数，右端有界）
//   f2   : K - floor(R*K/n), K=log2(log2 n)        （线性秩映射，丢弃幅度）
//   f3   : -log2(log2(1/p) + eps), p=(n-R+1)/n     （秩双对数）
//   f4   : -log2(log2(1/(F+1/n)) + eps)            （CDF 平滑秩双对数）
double basePriority(const std::string& method, const KeyInfo& ki, int x);

enum PerturbMode { NONE, FIXED, ADAPTIVE };

// 由基础映射 + 分桶 + 扰动构造最终 priority 向量
std::vector<double> buildPriorities(const std::string& method, const KeyInfo& ki, PerturbMode mode, std::mt19937_64& rng);

// 固定幅度扰动（E6 扫描用）：Q(f)=floor(f)，δ~U(0,phi)。phi 可跨过桶宽 1 以暴露跨桶逆序。
std::vector<double> buildPrioritiesPhi(const std::string& method, const KeyInfo& ki, double phi, std::mt19937_64& rng);

// Kendall tau：O(n log n)，tau=1 表示两个 priority 序完全一致（拓扑不变）
double kendallTau(const std::vector<double>& a, const std::vector<double>& b);

// 跨桶逆序率：只统计“不同桶”之间的 key 对里，优先级序与权重序相反的比例。
// 桶 = floor(basePriority)。它直接度量扰动把学习到的排序破坏了多少（扰动幅度 >= 桶间距时才会 >0）。
double crossBucketInversionRatio(const std::string& method, const KeyInfo& ki, const std::vector<double>& pri);
