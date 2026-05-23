#ifndef FIT_CPU_H
#define FIT_CPU_H

#include <vector>


// 6次関数のパラメータを保持する構造体
struct Pol6 {
    double a[7];
};

// ベースライン（直線）のフィッティング
void GetBaseline(const std::vector<double>& y, double& baseline);

// 6次関数を評価する関数
double evaluatePolynomial(const Pol6& poly, double x);

// 6次関数の最小二乗法によるフィッティング（xデータは自動生成）
void FitPol6(const std::vector<double>& y, Pol6& poly);

// フィッティングされた6次関数の最大値を探索
void GetPeak(const Pol6& poly, double x_min, double x_max, double& peak, int n_point);

#endif
