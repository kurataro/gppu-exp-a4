#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TF1.h>
#include <vector>
#include <TMath.h>

double GEXP_w_offset(double *x, double *par) {
    double t = x[0];
    double A = par[0];     // 振幅
    double t0 = par[1];    // 開始時刻
    double tau = par[2];   // 減衰時間定数
    double sigma = par[3]; // ガウスの標準偏差
    double offset = par[4];

    double exponent = sigma * sigma / (2 * tau * tau) - (t - t0) / tau;
    double erfc_arg = (sigma / (TMath::Sqrt(2) * tau)) - ((t - t0) / (sigma * TMath::Sqrt(2)));

    return offset + A * 0.5 * TMath::Exp(exponent) * TMath::Erfc(erfc_arg);
}

double two_GEXP_w_offset(double *x, double *par) {
    double t1 = x[0];
    double A1 = par[0];     // 振幅
    double t0_1 = par[1];    // 開始時刻
    double tau1 = par[2];   // 減衰時間定数
    double sigma1 = par[3]; // ガウスの標準偏差
    double t2 = x[1];
    double A2 = par[4];     // 振幅
    double t0_2 = par[5];    // 開始時刻
    double tau2 = par[6];   // 減衰時間定数
    double sigma2 = par[7]; // ガウスの標準偏差
    double offset = par[8];

    double exponent1 = sigma1 * sigma1 / (2 * tau1 * tau1) - (t1 - t0_1) / tau1;
    double erfc_arg1 = (sigma1 / (TMath::Sqrt(2) * tau1)) - ((t1 - t0_1) / (sigma1 * TMath::Sqrt(2)));
    double exponent2 = sigma2 * sigma2 / (2 * tau2 * tau2) - (t2 - t0_2) / tau2;
    double erfc_arg2 = (sigma2 / (TMath::Sqrt(2) * tau2)) - ((t2 - t0_2) / (sigma2 * TMath::Sqrt(2)));

    double GEXP1 = A1 * 0.5 * TMath::Exp(exponent1) * TMath::Erfc(erfc_arg1);
    double GEXP2 = A2 * 0.5 * TMath::Exp(exponent2) * TMath::Erfc(erfc_arg2);
    return offset + GEXP1 + GEXP2;
}

double asym_gaus_w_offset(double *x, double *par){
  double A = par[0];
  double mu = par[1];
  double sigma_left = par[2];
  double sigma_right = par[3];
  double offset = par[4];

  double asym_gaus;
  if(x[0] < mu){
    asym_gaus = A * TMath::Exp(-0.5 * TMath::Power((x[0] -mu)/ sigma_left, 2));
  }else {
    asym_gaus = A * TMath::Exp(-0.5 * TMath::Power((x[0] - mu)/ sigma_right,2));
  }
  return offset + asym_gaus;
}

void old_fit_waveform() {
    // ROOTファイルを開く
    TFile file("run01040_Rayraw.root");
    TTree *tree = (TTree*)file.Get("tree");

    // 2次元vectorを格納するポインタを用意してブランチをセット
    std::vector<std::vector<double>> *waveforms = nullptr;
    tree->SetBranchAddress("waveform", &waveforms);

    // 1イベント分のデータを取得（例: 最初のイベント）
    // tree->GetEntry(0);
    tree->GetEntry(10000);

    // 取得したデータから特定の波形（例: 1つ目の波形）を選択
    const std::vector<double> &waveform = waveforms->at(27); // 1つ目の波形を使用
    int nPoints = waveform.size();

    // TGraphを作成し、波形データをセット
    TGraph *graph = new TGraph(nPoints);
    for (int i = 0; i < nPoints; ++i) {
        graph->SetPoint(i, i, waveform[i]); // x軸はインデックス、y軸は波形の値
    }

    // グラフを描画
    graph->SetTitle("Waveform from 2D Vector with Fitting;Time index;Amplitude");
    graph->SetMarkerStyle(20);
    graph->Draw("AP");


    // フィット関数を定義（例: ガウシアン）
    // TF1 *fitFunc = new TF1("fitFunc", "gaus", 0, nPoints);
    Int_t x_min = 25;
    Int_t x_max = 40;
    // TF1 *fitFunc = new TF1("fitFunc", asym_gaus_w_offset, x_min, x_max, 5);
    // fitFunc->SetParameters(30, 25.0, 1.0, 2.0, 512);

    TF1 *fitFunc = new TF1("fitFunc", GEXP_w_offset, x_min, x_max, 5);
    fitFunc->SetParameters(30, 30.0, 10.0, 2.0, 512);
    // fitFunc->SetParameter(0, 550);
    // fitFunc->SetParameter(1, 30);
    // fitFunc->SetParameter(2, 5);
    // fitFunc->SetParameter(3, 2);

    // fitFunc->SetRange(25, 40);
    // fitFunc->SetRange(25, 50);

    // フィットを実行
    graph->Fit(fitFunc);
    fitFunc->Draw("same"); // フィット結果を重ねて表示

}
