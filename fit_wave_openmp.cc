// #include <iostream>
// #include <vector>
// #include <omp.h>
// #include "TFile.h"
// #include "TTree.h"
// #include "TGraph.h"
// #include "TF1.h"
// #include "TMath.h"

// // フィッティング関数
// double fitFunction(double *x, double *par) {
//     double t      = x[0];
//     double A      = par[0];
//     double t0     = par[1];
//     double tau    = par[2];
//     double sigma  = par[3];
//     double offset = par[4];

//     double exponent = sigma * sigma / (2 * tau * tau) - (t - t0) / tau;
//     double erfc_arg = (sigma / (TMath::Sqrt(2) * tau)) - ((t - t0) / (sigma * TMath::Sqrt(2)));

//     return offset + A * 0.5 * TMath::Exp(exponent) * TMath::Erfc(erfc_arg);
// }

// // 各イベントのデータをフィットする関数
// void fitEvent(const std::vector<double>& waveform, int waveformSize, double* chi2_ndf) {
//     double par[2] = {1.0, 0.1};

//     double chi2 = 0.0;
//     #pragma omp parallel for reduction(+:chi2)
//     for (int i = 0; i < waveformSize; ++i) {
//         double x = static_cast<double>(i);
//         double y = waveform[i];
//         double y_fit = fitFunction(&x, par);
//         chi2 += (y - y_fit) * (y - y_fit);
//     }

//     *chi2_ndf = chi2 / (waveformSize - 2);
// }

// int main() {
//     TFile file("run01040_Rayraw.root");
//     TTree *tree = (TTree*)file.Get("tree");

//     std::vector<std::vector<double>> *waveforms = nullptr;
//     tree->SetBranchAddress("waveform", &waveforms);

//     int n_event = tree->GetEntries();
//     int n_loop  = 10;
//     std::vector<double> chi2_ndf(n_loop);

//     double start = omp_get_wtime();

//     // OpenMPによる並列化
//     #pragma omp parallel for
//     for (int i = 0; i < n_loop; ++i) {
//         tree->GetEntry(i);
//         int n_ch = waveforms->size();

//         for (int j = 0; j < n_ch; ++j) {
//             const std::vector<double>& waveform = waveforms->at(j);
//             int n_sample = waveform.size();

//             TGraph *gr = new TGraph(n_sample);
//             for (int k = 0; k < n_sample; ++k) {
//                 gr->SetPoint(k, k, waveform[k]);
//             }

//             gr->SetTitle(Form("run01040_seg%d_waveform", j));
//             // gr->GetXaxis()->SetTitle("sample #");
//             // gr->GetYaxis()->SetTitle("ADC [ch]");
//             gr->SetMarkerStyle(20);
//             gr->Draw("AP");

//             // フィッティング関数の設定
//             Int_t  x_min         = 25;
//             Int_t  x_max         = 40;
//             double A_ini        = 30.0;
//             double t0_ini       = 30.0;
//             double tau_ini      = 10.0;
//             double sigma_ini    = 2.0;
//             double baseline_ini = 512;

//             TF1 *fitFunc = new TF1("fitFunc", fitFunction, x_min, x_max, 5);
//             fitFunc->SetParameters(A_ini, t0_ini, tau_ini, sigma_ini, baseline_ini);

//             // フィッティングを実行
//             gr->Fit(fitFunc);
//             fitFunc->Draw("same");
//         }
//     }

//     double end = omp_get_wtime();

//     std::cout << "Processing time: " << end - start << std::endl;
//     // // 結果表示
//     // for (int i = 0; i < n_loop; ++i) {
//     //     std::cout << "Event " << i << " chi2/ndf: " << chi2_ndf[i] << std::endl;
//     // }

//     return 0;
// }

#include <iostream>
#include <vector>
#include <omp.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMath.h>
#include <ctime>

// フィッティング関数
double GEXP_w_offset(double *x, double *par) {
    double t      = x[0];
    double A      = par[0]; // Amplitude
    double t0     = par[1]; // Rising time
    double tau    = par[2]; // Decay time constant
    double sigma  = par[3]; // std. dev. of gaussian
    double offset = par[4]; // baseline

    double exponent = sigma * sigma / (2 * tau * tau) - (t - t0) / tau;
    double erfc_arg = (sigma / (TMath::Sqrt(2) * tau)) - ((t - t0) / (sigma * TMath::Sqrt(2)));

    return offset + A * 0.5 * TMath::Exp(exponent) * TMath::Erfc(erfc_arg);
}

// 各イベントのデータをフィットする関数
void fitEvent(const std::vector<double>& waveform, int waveformSize, double* chi2_ndf) {
    // フィットパラメータの初期値
    double par[5] = {1.0, 30.0, 10.0, 2.0, 0.1}; // A, t0, tau, sigma, offset

    // chi2計算
    double chi2 = 0.0;
    for (int i = 0; i < waveformSize; ++i) {
        double x = static_cast<double>(i); // x軸の値
        double y = waveform[i];             // 実際の波形データ
        double y_fit = GEXP_w_offset(&x, par);
        chi2 += (y - y_fit) * (y - y_fit); // chi2を計算
    }

    // 最小二乗法によるフィッティングの結果を保存
    *chi2_ndf = chi2 / (waveformSize - 5); // chi2/ndfの計算
}

int main() {
    TFile file("run01040_Rayraw.root");
    TTree *tree = (TTree*)file.Get("tree");

    std::vector<std::vector<double>> *waveforms = nullptr; // 2次元ベクタで波形を格納
    tree->SetBranchAddress("waveform", &waveforms);

    int n_event = tree->GetEntries();

    // chi2/ndfの結果を保存する配列
    std::vector<double> chi2_ndf(n_event);

    // 計測開始
    clock_t start = clock();

    // OpenMPを使用して各イベントのフィッティングを並列処理
    #pragma omp parallel for
    for (int i = 0; i < n_event; ++i) {
        tree->GetEntry(i);
        int n_ch = waveforms->size();

        // 各チャネルの波形データをフィッティング
        for (int j = 0; j < n_ch; ++j) {
            const std::vector<double>& waveform = waveforms->at(j);
            int n_sample = waveform.size();

            // chi2/ndfを保存するための変数
            double chi2_ndf_value;

            // fitEvent関数を呼び出し
            fitEvent(waveform, n_sample, &chi2_ndf_value);

            // chi2/ndfの結果を保存
            chi2_ndf[i] += chi2_ndf_value;
        }
    }

    // 計測終了
    clock_t end = clock();
    double elapsed_time = static_cast<double>(end - start) / CLOCKS_PER_SEC;

    // 結果表示
    for (int i = 0; i < n_event; ++i) {
        std::cout << "Event " << i << " chi2/ndf: " << chi2_ndf[i] << std::endl;
    }

    std::cout << "Total elapsed time: " << elapsed_time << " seconds." << std::endl;

    return 0;
}
