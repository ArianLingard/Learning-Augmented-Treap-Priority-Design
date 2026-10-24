#pragma once

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

// ====================================================
// 统一模板：打开文件 / 记录计时 / 写文件 / 记录运行时间
//   runExperiment(id, path, fn, log)：创建 ofstream -> 计时执行 fn(out) -> 关闭，
//   在控制台打印耗时，并把一行 "id,path,seconds" 写入 log。
//   所有实验（E1-E7）共用这套打开/计时/写文件逻辑。
// ====================================================

// 简单计时器：返回从起点 t0 到现在的秒数（高精度时钟，微秒级）
inline double elapsedSeconds(std::chrono::high_resolution_clock::time_point t0) {
    return std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::high_resolution_clock::now() - t0).count() / 1e6;
}

template <typename Fn>
void runExperiment(const std::string& id, const std::string& path, Fn&& fn, std::ofstream& log) {
    auto t0 = std::chrono::high_resolution_clock::now();
    std::ofstream out(path.c_str());
    if (!out.is_open()) {
        std::cerr << "  [error] cannot open " << path << "\n";
        return;
    }
    fn(out);
    out.close();
    double secs = elapsedSeconds(t0);
    std::cout << "  wrote " << path << "  ["
              << std::fixed << std::setprecision(2) << secs << " s]\n";
    log << id << "," << path << "," << std::fixed << std::setprecision(2) << secs << "\n";
}
