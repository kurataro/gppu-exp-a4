#include <iostream>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <omp.h>

#include "FitCPU.hh"

#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define MAGENTA "\033[35m"
#define CYAN    "\033[36m"
#define RESET   "\033[0m"


int main(int argc, char* argv[]){

  if(argc == 1)
    std::cout << "#of event is not specified " << MAGENTA << "-> All events will be analyzed!" << RESET << std::endl;

  TString ifname = "/group/had/sks/Users/rkurata/rootfile/rayraw/20241031_gaincheck_BIAS0_recalc/run01040_Rayraw.root";
  TFile *ifile = TFile::Open(ifname);
  if(!ifile || ifile->IsZombie()){
    std::cerr << RED << "#E Cannot open the file: " << RESET << ifname << std::endl;
    return 1;
  }
  std::cout << CYAN << "#D File opened: " << RESET << ifname << std::endl;

  TString itname = "tree";
  TTree *itree = static_cast<TTree*>(ifile->Get(itname));
  if(!itree){
    std::cerr << RED << "#E Cannot find the tree: " << RESET << itname << std::endl;
    return 1;
  }
  std::cout << CYAN << "#D Tree loaded: " << RESET << itname << std::endl;

  std::vector<std::vector<Double_t>> *waveforms = nullptr;
  itree->SetBranchAddress("waveform", &waveforms);

  // Determine #of event fot analyzing
  int n_loop;
  if(argc == 2){
    n_loop = atoi(argv[1]);
  }else{
    n_loop = itree->GetEntries();
  }
  static const int n_event = n_loop;

  // Vectors for store the fitting results
  std::vector<double> Contbaseline;
  std::vector<double> Contpeak;
  // std::vector<double> Contchi2_ndf;

  Contbaseline.reserve(n_event * waveforms->size());
  Contpeak.reserve(n_event * waveforms->size());

  TStopwatch stopwatch;
  stopwatch.Start();
  for(int ev=0; ev<n_event; ++ev){
    if(ev%10000 == 0)
    // if(ev%100 == 0)
      std::cout << GREEN << "#D [for loop] event#: " << RESET << ev << std::endl;

    itree->GetEntry(ev);

    std::vector<double> local_baseline;
    std::vector<double> local_peak;

#pragma omp parallel for
    for(const auto& waveform : *waveforms){
    // for(size_t w =0; w<waveforms->size(); ++w){
    //   const auto& waveform = (*waveforms)[w];
      // std::vector<double> waveform = (*waveforms)[w];

      double baseline;
      GetBaseline(waveform, baseline);
      local_baseline.push_back(baseline);
      // Contbaseline.push_back(baseline);

      std::vector<double> waveform_crct(waveform.size());

      #pragma omp parallel for
	for(size_t i=0; i<waveform.size(); ++i){
	  waveform_crct[i] = waveform[i] - baseline;
	}

      Pol6 poly;
      FitPol6(waveform_crct, poly);

      double peak;
      double x_min=25;
      double x_max=40;
      int    n_point=(x_max - x_min +1)*100;

      GetPeak(poly, x_min, x_max, peak, n_point);

      // Contpeak.push_back(peak);
      local_peak.push_back(peak);
      // Contchi2_ndef.push_back(chi2_nddef);

      // std::vector<double> fitResult = least_square(waveform); // Individual Fitting results
      // baseline.push_back(fitResult[0]);
      // peak.push_back(fitResult[1]);
      // peak_baseline.push_back(fitResult[2]);
      // chi2_ndef.push_back(fitResult[3]);
    }

    #pragma omp critical
    {
        Contbaseline.insert(Contbaseline.end(), local_baseline.begin(), local_baseline.end());
        Contpeak.insert(Contpeak.end(), local_peak.begin(), local_peak.end());
    }

  }
  stopwatch.Stop();
  // auto end = std::chrono::high_resolution_clock::now();
  std::cout << CYAN << "#D for loop done" << RESET << std::endl;

  TString ofname = "fit_cpu.root";
  TFile *ofile = TFile::Open(ofname, "RECREATE");
  TTree *otree = new TTree("tree", "Tree for Fitting Results");
  std::cout << CYAN << "#D output file & tree created: " << RESET << ofname << std::endl;

  otree->Branch("baseline",      &Contbaseline);
  otree->Branch("peak",          &Contpeak);
  // otree->Branch("chi2/ndf",      &Contchi2_ndf);

  otree->Fill();
  otree->Write();
  ofile->Close();
  std::cout << CYAN << "#D output file closed: " << RESET << ofname << std::endl;

  ifile->Close();
  std::cout << CYAN << "#D input file closed: " << RESET << ifname << std::endl;

  std::cout << MAGENTA << "Total processing time for fit (Real): " << stopwatch.RealTime() << "s" << RESET << std::endl;
  std::cout << MAGENTA << "Total processing time for fit (CPU): " << stopwatch.CpuTime() << "s" << RESET << std::endl;

  return 0;
}
