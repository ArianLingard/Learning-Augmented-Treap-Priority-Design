#include "oracle.h"

#include <cmath>

using namespace std;

double klDivergence(const vector<double>& P, const vector<double>& Q) {
    double kl = 0.0;
    for (size_t i = 0; i < P.size(); i++) if (P[i] > 0) kl += P[i] * log2(P[i] / Q[i]);
    return kl;
}

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
