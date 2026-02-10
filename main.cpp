#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <algorithm>

using namespace std;

struct DataPoint { double x, t, u; };
struct Alg2Exception {};
struct Alg3Exception {};

map<string, vector<DataPoint>> cache;

void load_file(string name) {
    if (cache.count(name)) return;
    ifstream file(name);
    if (!file.is_open()) return;
    double cx, ct, cu;
    while (file >> cx >> ct >> cu) cache[name].push_back({ cx, ct, cu });
    sort(cache[name].begin(), cache[name].end(), [](auto& a, auto& b) { return a.x < b.x; });
}

void get_val(double x, double& T, double& U) {
    string file = abs(x) <= 1.0 ? "dat_X_1_1.dat" : (x < -1.0 ? "dat_X_00_1.dat" : "dat_X_1_00.dat");
    double tx = abs(x) <= 1.0 ? x : 1.0 / x;

    load_file(file);
    auto& data = cache[file];
    if (data.empty()) throw Alg3Exception();

    auto it = lower_bound(data.begin(), data.end(), tx, [](auto& p, double val) { return p.x < val; });

    if (it != data.end() && abs(it->x - tx) < 1e-9) {
        T = it->t; U = it->u;
    }
    else if (it != data.begin() && it != data.end()) {
        auto prev = it - 1;
        double ratio = (tx - prev->x) / (it->x - prev->x);
        T = prev->t + (it->t - prev->t) * ratio;
        U = prev->u + (it->u - prev->u) * ratio;
    }
    else {
        T = (it == data.begin()) ? data.front().t : data.back().t;
        U = (it == data.begin()) ? data.front().u : data.back().u;
    }
}

double Srz(double x, double y, double z) {
    double tx, ux, ty, uy, tz, uz;
    get_val(x, tx, ux); get_val(y, ty, uy); get_val(z, tz, uz);
    return (x > y) ? tx + uz - ty : ty + uy - uz;
}

double Gold1(double x, double y) {
    if (y == 0) return 0.1;
    if (x > y && abs(y) > 0.1) return x / y;
    if (x <= y && abs(x) > 0.1) return y / x;
    if (x < y && abs(x) > 0.1) return 0.15;
    return 0;
}

double Grs1(double x, double y) {
    auto Glr1 = [](double a, double b) { return abs(a) < 1 ? a : b; };
    return 0.14 * Srz(x + y, Gold1(x, y), Glr1(x, x * y)) +
        1.83 * Srz(x - y, Gold1(y, x / 5.0), Glr1(4.0 * x, x * y)) +
        0.83 * Srz(x, Glr1(y, x / 4.0), Gold1(4.0 * y, y));
}

double fun_alg3(double x, double y, double z) {
    return 1.3498 * z + 2.2362 * y - 2.348 * x * y;
}

double fun_alg2(double x, double y, double z) {
    try {
        return x * Grs1(x, y) + y * Grs1(y, z) + z * Grs1(z, x);
    }
    catch (Alg3Exception&) { return fun_alg3(x, y, z); }
}

double Gold(double x, double y) {
    if ((x > y && y == 0) || (x < y && x == 0)) throw Alg2Exception();
    return (x > y) ? x / y : y / x;
}

double Glr(double x, double y) {
    if (abs(x) < 1) return x;
    if (abs(y) < 1) return y;
    double d = x * x + y * y - 4;
    if (d >= 0 && sqrt(d) > 0.1) return y / sqrt(d);
    throw Alg2Exception();
}

double Grs(double x, double y) {
    return 0.1389 * Srz(x + y, Gold(x, y), Glr(x, x * y)) +
        1.8389 * Srz(x - y, Gold(y, x / 5.0), Glr(5.0 * x, x * y)) +
        0.83 * Srz(x - 0.9, Glr(y, x / 5.0), Gold(5.0 * y, y));
}

int main() {
    double x, y, z;
    if (!(cin >> x >> y >> z)) return 1;
    try {
        cout << (x * x * Grs(y, z) + y * y * Grs(x, z) + 0.33 * x * y * Grs(x, z)) << endl;
    }
    catch (Alg2Exception&) {
        cout << fun_alg2(x, y, z) << endl;
    }
    catch (Alg3Exception&) {
        cout << fun_alg3(x, y, z) << endl;
    }
    return 0;
}