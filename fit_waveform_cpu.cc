// #include <TFile.h>
// #include <TTree.h>
// #include <TGraph.h>
// #include <TF1.h>
// #include <vector>
// #include <TMath.h>

#include <iostream>
#include <vector>
#include <cmath>

// double GEXP_w_offset(double *x, double *par) {
double fitFunction(double *x, double *par) {
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

// // フィッティング関数（例: 単純な指数関数）
// double fitFunction(double x, double *par) {
//     return par[0] * exp(-par[1] * x);
// }

// 各イベントのデータをフィットする関数
void fitEvent(const std::vector<double>& waveform, int waveformSize, double* chi2_ndf) {
    // フィットパラメータの初期値（例として2つのパラメータを持つ）
    double par[2] = {1.0, 0.1};

    // フィッティングの反復処理（ここでは単純にchi2計算のみを例示）
    double chi2 = 0.0;
    for (int i = 0; i < waveformSize; ++i) {
        double x = static_cast<double>(i); // x軸の値（例: インデックスをそのまま使用）
        double y = waveform[i];             // 実際の波形データ
        double y_fit = fitFunction(x, par);
        chi2 += (y - y_fit) * (y - y_fit);  // chi2を計算
    }

    // 最小二乗法によるフィッティングの結果を保存
    *chi2_ndf = chi2 / (waveformSize - 2); // chi2/ndfの計算
}

int main() {
    int numEvents = 10;         // イベント数
    int waveformSize = 80;      // 波形のサイズ
    std::vector<double> waveforms(numEvents * waveformSize, 0.1); // ダミー波形データ

    // chi2/ndfの結果を保存する配列
    std::vector<double> chi2_ndf(numEvents);

    // 各イベントごとにフィッティングを実行
    for (int eventIdx = 0; eventIdx < numEvents; ++eventIdx) {
        // イベントごとの波形データを取得
        std::vector<double> waveform(waveforms.begin() + eventIdx * waveformSize,
                                     waveforms.begin() + (eventIdx + 1) * waveformSize);

        // chi2/ndfの計算
        fitEvent(waveform, waveformSize, &chi2_ndf[eventIdx]);
    }

    // 結果表示
    for (int i = 0; i < 10; ++i) {
        std::cout << "Event " << i << " chi2/ndf: " << chi2_ndf[i] << std::endl;
    }
    return 0;
}

// void fit_waveform_ROOT() {

    TFile file("run01040_Rayraw.root");
    TTree *tree = (TTree*)file.Get("tree");

    std::vector<std::vector<double>> *waveforms = nullptr; // 2-dim vector for store the tree waveform
    tree->SetBranchAddress("waveform", &waveforms);

    int n_event = tree->GetEntries();
    int n_loop  = n_event;

    for(int i=0; i<n_loop; ++i){
      tree->GetEntry(i);
      int n_ch = waveforms->size();
      std::cout << "n_ch: " << n_ch << std::endl;
      for(int j=0; j<n_ch; ++j){
	const std::vector<double> &waveform = waveforms->at(j);
	int n_sample = waveform.size();

	TGraph *gr = new TGraph(n_sample);
	for(int k=0; k<n_sample; ++k){
	  gr->SetPoint(k, k, waveform[k]);
	}

	gr->SetTitle(Form("run01040_seg%d_waveform", j));
	gr->GetXaxis()->SetTitle("sample #");
	gr->GetYaxis()->SetTitle("ADC [ch]");
	gr->SetMarkerStyle(20);
	gr->Draw("AP");

	// Define fitting function
	Int_t  x_min         = 25;
	Int_t  x_max         = 40;
	double A_ini        = 30.0;
	double t0_ini       = 30.0;
	double tau_ini      = 10.0;
	double sigma_ini    = 2.0;
	double baseline_ini = 512;

	TF1 *fitFunc = new TF1("fitFunc", GEXP_w_offset, x_min, x_max, 5);
	fitFunc->SetParameters(A_ini, t0_ini, tau_ini, sigma_ini, baseline_ini);

	// Do fitting
	gr->Fit(fitFunc);

	//Draw fitting result
	fitFunc->Draw("same");
      }
    }
// }
