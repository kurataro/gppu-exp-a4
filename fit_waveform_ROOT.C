#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TF1.h>
#include <vector>
#include <TMath.h>
#include <TStopwatch.h>

// double GEXP_w_offset(double *x, double *par) {
//     double t      = x[0];
//     double A      = par[0]; // Amplitude
//     double t0     = par[1]; // Rising time
//     double tau    = par[2]; // Decay time constant
//     double sigma  = par[3]; // std. dev. of gaussian
//     double offset = par[4]; // baseline

//     double exponent = sigma * sigma / (2 * tau * tau) - (t - t0) / tau;
//     double erfc_arg = (sigma / (TMath::Sqrt(2) * tau)) - ((t - t0) / (sigma * TMath::Sqrt(2)));

//     return offset + A * 0.5 * TMath::Exp(exponent) * TMath::Erfc(erfc_arg);
// }

void fit_waveform_ROOT() {

    TFile file("run01040_Rayraw.root");
    TTree *tree = (TTree*)file.Get("tree");

    std::vector<std::vector<double>> *waveforms = nullptr; // 2-dim vector for store the tree waveform
    tree->SetBranchAddress("waveform", &waveforms);

    int n_event = tree->GetEntries();
    // int n_loop  = 1;
    int n_loop  = n_event;


    std::vector<std::vector<double>> chi2_ndf_all;

    TStopwatch stopwatch;
    stopwatch.Start();

    for(int i=0; i<n_loop; ++i){
      if(i%10000 == 0)
	std::cout << "loop: " << i << std::endl;

      // tree->GetEntry(i+1001);
      tree->GetEntry(i);
      int n_ch = waveforms->size();
      // int n_ch = 1;
      // std::cout << "n_ch: " << n_ch << std::endl;

      std::vector<double> chi2_ndf; //for each channel

      for(int j=0; j<n_ch; ++j){
	const std::vector<double> &waveform = waveforms->at(j);
	int n_sample = waveform.size();

	TGraph *gr = new TGraph(n_sample);
	for(int k=0; k<n_sample; ++k){
	  gr->SetPoint(k, k, waveform[k]);
	}

	// gr->SetTitle(Form("run01040_seg%d_waveform", j));
	// gr->GetXaxis()->SetTitle("sample #");
	// gr->GetYaxis()->SetTitle("ADC [ch]");
	// gr->SetMarkerStyle(20);
	// gr->Draw("AP");

	// Define fitting function
	// Int_t  x_min         = 25;
	// Int_t  x_max         = 40;
	// double A_ini        = 30.0;
	// double t0_ini       = 30.0;
	// double tau_ini      = 10.0;
	// double sigma_ini    = 2.0;
	// double baseline_ini = 512;

	Int_t x_min1 = 60;
	Int_t x_max1 = 80;

	// TF1 *fitpol1 = new TF1("fitpol1", "pol1" , x_min1, x_max1);
	// fitFunc->SetParameters(A_ini, t0_ini, tau_ini, sigma_ini, baseline_ini);
	// gr->Fit(fitpol1);
	gr->Fit("pol1", "Q", "", x_min1, x_max1);
	TF1 *fitpol1 = gr->GetFunction("pol1"); // "pol1" のフィット関数オブジェクト
	double offset = fitpol1->GetParameter(0);
	// double offset = fitpol1->GetParameter(0);
	// double offset = 511;
	// std::cout << "offset: " << offset << std::endl;

	TGraph *gr2 = new TGraph(n_sample);
	for(int k=0; k<n_sample; ++k){
	  // gr2->SetPoint(k, k, waveform[k] - offset);
	  gr2->SetPoint(k, k, waveform[k] - offset);
	}
	// gr2->SetTitle(Form("run01040_seg%d_waveform2", j));
	// gr2->GetXaxis()->SetTitle("sample #");
	// gr2->GetYaxis()->SetTitle("ADC [ch]");
	// gr2->SetMarkerStyle(20);
	// gr2->Draw("AP");

	Int_t x_min2 = 25;
	Int_t x_max2 = 40;

	// TF1 *fitpol4 = new TF1("fitpol4", "pol4" , x_min2, x_max2);
	// fitFunc->SetParameters(A_ini, t0_ini, tau_ini, sigma_ini, baseline_ini);
	// Do fitting
	// gr2->Fit(fitpol4);
	gr2->Fit("pol6", "Q", "", x_min2, x_max2);
	// TF1 *fitpol2 = gr->GetFunction("pol4"); // "pol1" のフィット関数オブジェクト
	// fitpol2->Draw("same");

	//Draw fitting result
	// fitpol4->Draw("same");

	// TF1 *fitFunc = new TF1("fitFunc", GEXP_w_offset, x_min, x_max, 5);
	// fitFunc->SetParameters(A_ini, t0_ini, tau_ini, sigma_ini, baseline_ini);

	// // Do fitting
	// gr->Fit(fitFunc);

	// double chi2 = fitFunc->GetChisquare();
	// int    ndf  = fitFunc->GetNDF();
	// chi2_ndf.push_back(chi2 / ndf);

	// //Draw fitting result
	// fitFunc->Draw("same");
      }
      // chi2_ndf_all.push_back(chi2_ndf);
    }

    stopwatch.Stop();

    std::cout << "n_fit:" << n_loop * 32 << ", Processing Time)  Real: " << stopwatch.RealTime() << "[s], CPU: " << stopwatch.CpuTime() << "[s]" << std::endl;
}
