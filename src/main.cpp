#include "experiment.h"
#include "runner.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
    vector<string> run;
    for (int i = 1; i < argc; i++) run.push_back(argv[i]);
    bool all = run.empty();
    auto shouldRun = [&](const string& id) {
        return all || find(run.begin(), run.end(), id) != run.end();
    };

    cout << "Learning-Augmented Treap experiments (E1-E7)\n";
    auto tAll = chrono::high_resolution_clock::now();

    // 每次运行生成 data/run_times.csv，记录每个实验的耗时（experiment,path,seconds）
    ofstream log("data/run_times.csv");
    log << "experiment,path,seconds\n";

    if (shouldRun("e1")) runExperiment("e1", "data/e1_basic.csv", runE1, log);
    if (shouldRun("e2")) runExperiment("e2", "data/e2_drift.csv", runE2, log);
    if (shouldRun("e3")) runExperiment("e3", "data/e3_singularity.csv", runE3, log);
    if (shouldRun("e4")) runExperiment("e4", "data/e4_static_optimal.csv", runE4, log);
    if (shouldRun("e5")) runExperiment("e5", "data/e5_kl_noise.csv", runE5, log);
    if (shouldRun("e6")) runExperiment("e6", "data/e6_adaptive_perturb.csv", runE6, log);
    if (shouldRun("e7")) runExperiment("e7", "data/e7_bucketing.csv", runE7, log);

    double total = elapsedSeconds(tAll);
    cout << "done. total " << fixed << setprecision(2) << total << " s\n";
    log << "total,," << fixed << setprecision(2) << total << "\n";
    log.close();
    return 0;
}
