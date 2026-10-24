#pragma once

#include <algorithm>
#include <random>
#include <vector>

// ====================================================
// 模块 3：预测器 Oracle（完美频率 / 量纲漂移 / KL 噪声）
// ====================================================
// 完美频率：w_x = p_x（直接使用）
// 缩放/幂律（E2）：w_x = c*p_x 或 (p_x)^beta（严格单调变换），
// 双对数定义域为 (0,1)：把可能越界的预测权重夹回 [1e-300, 1-1e-12]
inline double clampWeight(double x) { return std::min(std::max(x, 1e-300), 1.0 - 1e-12); }

double klDivergence(const std::vector<double>& P, const std::vector<double>& Q);

// 逐键乘性对数正态噪声：Q ∝ P·exp(σ z)，z~N(0,1)，二分搜索 σ 使 D_KL(P||Q) 逼近 target。
// 乘性误差既扰动幅度也扰动相对排序（比纯温度缩放更贴近真实预测器）。
std::vector<double> makeNoisy(const std::vector<double>& P, double targetKL, std::mt19937_64& rng);
