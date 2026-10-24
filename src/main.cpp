#include "experiment.h"
#include "runner.h"

#include <algorithm>
#include <chrono>
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

    if (shouldRun("e1")) runExperiment("data/e1_basic.csv", runE1);
    if (shouldRun("e2")) runExperiment("data/e2_drift.csv", runE2);
    if (shouldRun("e3")) runExperiment("data/e3_singularity.csv", runE3);
    if (shouldRun("e4")) runExperiment("data/e4_static_optimal.csv", runE4);
    if (shouldRun("e5")) runExperiment("data/e5_kl_noise.csv", runE5);
    if (shouldRun("e6")) runExperiment("data/e6_adaptive_perturb.csv", runE6);
    if (shouldRun("e7")) runExperiment("data/e7_bucketing.csv", runE7);

    cout << "done. total " << fixed << setprecision(2) << elapsedSeconds(tAll) << " s\n";
    return 0;
}
