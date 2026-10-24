#pragma once

#include <random>
#include <string>
#include <vector>

// ====================================================
// 模块 2：访问分布（键序随机打乱 -> 权重按随机排列赋给键）
// ====================================================
// 返回长度为 n 的概率向量 p（按 rank 从热到冷排列：p[0] 最热）
std::vector<double> genDistribution(const std::string& name, double param, int n);

// 键序随机打乱：把 rank 概率 p 随机赋给键 0..n-1
std::vector<double> assignWeights(const std::vector<double>& p, int n, std::mt19937_64& rng);

// 依据权重向量 w 采样长度为 m 的访问序列（同一组实验内所有方法共用）
std::vector<int> genAccessSequence(const std::vector<double>& w, int m, std::mt19937_64& rng);
