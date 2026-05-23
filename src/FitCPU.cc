#include "FitCPU.hh"
#include <cmath>
#include <algorithm>

const int ORDER = 6;

// Get baseline values
void GetBaseline(const std::vector<double>& y, double& baseline) {
  int    n = y.size();
  double sumY = 0.0;
  int    count = 0;

  // 先頭10サンプルのy値を使用して平均値を計算
  for (int i = 0; i < 10 && i < n; ++i) {
    sumY += y[i];
    count++;
  }

  // 最後の10サンプルのy値を使用して平均値を計算
  for (int i = std::max(0, n - 10); i < n; ++i) {
    sumY += y[i];
    count++;
  }

  // ベースラインの平均値を slope に格納
  baseline = sumY / count;
}

// 6次関数の評価
double evaluatePolynomial(const Pol6& poly, double x) {
    double result = 0.0;
    double xi = 1.0;
    for (int i = 0; i <= ORDER; ++i) {
        result += poly.a[i] * xi;
        xi *= x;
    }
    return result;
}

// 6次関数の最小二乗法によるフィッティング（xデータは自動生成）
void FitPol6(const std::vector<double>& y, Pol6& poly) {
    int n = y.size();
    double A[ORDER + 1][ORDER + 1] = {0.0};
    double B[ORDER + 1] = {0.0};

    // 行列AとベクトルBの計算
    for (int i = 0; i < n; ++i) {
        double xi = 1.0;
        for (int j = 0; j <= ORDER; ++j) {
            B[j] += y[i] * xi;
            double xj = xi;
            for (int k = 0; k <= ORDER; ++k) {
                A[j][k] += xj;
                xj *= i;  // xデータはインデックス値として扱う
            }
            xi *= i;
        }
    }

    // ガウスの消去法で連立方程式を解く
    for (int i = 0; i <= ORDER; ++i) {
        for (int j = i + 1; j <= ORDER; ++j) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k <= ORDER; ++k) {
                A[j][k] -= factor * A[i][k];
            }
            B[j] -= factor * B[i];
        }
    }
    for (int i = ORDER; i >= 0; --i) {
        poly.a[i] = B[i];
        for (int j = i + 1; j <= ORDER; ++j) {
            poly.a[i] -= A[i][j] * poly.a[j];
        }
        poly.a[i] /= A[i][i];
    }
}

// フィッティングされた6次関数の最大値を探索
void GetPeak(const Pol6& poly, double x_min, double x_max, double& peak, int n_point) {
    peak = -1e20;

    double dx = (x_max - x_min) / n_point;

    for (int i = 0; i < n_point; ++i) {
        double x = x_min + i * dx;
        double y = evaluatePolynomial(poly, x);

        if (y > peak) {
            peak = y;
        }
    }
}
