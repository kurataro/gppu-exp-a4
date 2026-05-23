#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <fstream>
#include <iostream>

// int main() {
int create_txtdata() {
    TFile file("run01040_Rayraw.root");
    TTree *tree = (TTree*)file.Get("tree");

    std::vector<std::vector<double>> *waveforms = nullptr;
    tree->SetBranchAddress("waveform", &waveforms);

    int n_event = tree->GetEntries();
    int n_loop  = 100;

    // 出力用のテキストファイル
    std::ofstream outFile("run01040_waveform_data.txt");
    if (!outFile) {
        std::cerr << "Failed to open output file!" << std::endl;
        return -1;
    }

    // データのエクスポート
    for (int i = 0; i < n_loop; ++i) {
        tree->GetEntry(i);

        int n_ch = waveforms->size();
        for (int j = 0; j < n_ch; ++j) {
            const std::vector<double> &waveform = waveforms->at(j);
            outFile << i << " " << j;  // イベント番号とチャンネル番号を書き出し

            for (const auto &sample : waveform) {
                outFile << " " << sample;  // サンプルデータを書き出し
            }
            outFile << "\n";  // 行の終わり
        }
    }

    outFile.close();
    file.Close();

    std::cout << "Data export completed!" << std::endl;
    return 0;
}
