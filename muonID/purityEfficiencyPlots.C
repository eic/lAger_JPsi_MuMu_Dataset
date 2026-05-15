#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <numeric>
#include <TFile.h>
#include <TH1.h>
#include <TH3.h>
#include <TBranch.h>
#include <TH2.h>
#include <TTree.h>
#include <TChain.h>
#include <TCut.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TRandom.h>
#include <TEventList.h>
#include <TMultiLayerPerceptron.h>
#include <TComplex.h>
#include <TVirtualGeoPainter.h>
#include <TFile.h>
#include <TSystem.h>
#include <TClassTree.h>
#include <TPaveLabel.h>
#include <TCanvas.h>
#include <TGClient.h>
#include <RQ_OBJECT.h>
#include <TApplication.h>
#include <TRint.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TF1.h>
#include <TF2.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <Riostream.h>
#include <TObjString.h>
#include <TChain.h>
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TLatex.h"
#include "TPie.h"

#include "../ePICStyle.C"


void purityEfficiencyPlots()
{
    //gROOT->SetBatch(kTRUE);
    gROOT->ProcessLine("SetePICStyle()");
    gStyle->SetOptStat(0);
    gStyle->SetCanvasPreferGL(true); // enables transparency

    std::vector<int> sixColourScheme = {
    TColor::GetColor("#5790fc"),     // blue
    TColor::GetColor("#e42536"),     // red
    TColor::GetColor("#7021dd"),     // violet
    TColor::GetColor("#964a8b"),     // grape
    TColor::GetColor("#f89c20"),     // yellow
    TColor::GetColor("#9c9ca1"),     // grey
    };

    std::vector<double> purity;
    std::vector<double> efficiency;
    std::vector<double> likelihoodThresholds;
    double JPsiContent = 0.0;
    double DISContent = 0.0;

    TString infileJPsi = "outputs/SimCampaign_26020_JPsiMuMu_10ifb_10x130ep_likelihoodPlots.root";
    TFile *inFileJPsi = new TFile(infileJPsi);

    TString infileDIS = "outputs/DIS_Q2_1_100_10x130ep_likelihoodPlots.root";
    TFile *inFileDIS = new TFile(infileDIS);

    TH1D* h_JPsiLikelihood = (TH1D*)inFileJPsi->Get("h_muonLikelihood");
    TH1D* h_DISLikelihood = (TH1D*)inFileDIS->Get("h_nonMuonLikelihood");



    for (int i = h_JPsiLikelihood->GetNbinsX(); i >= 1; i--) {
        JPsiContent += h_JPsiLikelihood->GetBinContent(i);
        DISContent += h_DISLikelihood->GetBinContent(i);
        if (JPsiContent + DISContent > 0) purity.push_back(JPsiContent / (JPsiContent + DISContent));
        else purity.push_back(0.0);
        efficiency.push_back(JPsiContent / h_JPsiLikelihood->Integral());
        likelihoodThresholds.push_back(h_JPsiLikelihood->GetBinCenter(i));

        std::cout << "Likelihood Threshold: " << h_JPsiLikelihood->GetBinCenter(i) << " | Purity: " << purity.back() << " | Efficiency: " << efficiency.back() << std::endl;
    }

    TGraph* g_purity = new TGraph(purity.size(), &likelihoodThresholds[0], &purity[0]);
    g_purity->SetName("g_purity");
    g_purity->SetTitle("Purity vs Likelihood Threshold;Likelihood Threshold;Purity");

    TGraph* g_efficiency = new TGraph(efficiency.size(), &likelihoodThresholds[0], &efficiency[0]);
    g_efficiency->SetName("g_efficiency");
    g_efficiency->SetTitle("Efficiency vs Likelihood Threshold;Likelihood Threshold;Efficiency");

    TLegend* legend_muPi = new TLegend(0.2, 0.6, 0.4, 0.8);
    legend_muPi->AddEntry(h_JPsiLikelihood, "JPsi->MuMu", "P");
    legend_muPi->AddEntry(h_DISLikelihood, "DIS Background", "P");
    legend_muPi->SetBorderSize(0);
    legend_muPi->SetFillStyle(0);
    legend_muPi->SetTextSize(0.04);

    TCanvas* c_likelihood = new TCanvas("c_likelihood", "c_likelihood", 800, 600);
    h_JPsiLikelihood->SetMarkerStyle(20);
    h_JPsiLikelihood->SetMarkerColor(sixColourScheme[0]);
    h_JPsiLikelihood->SetMarkerSize(1.75);
    h_JPsiLikelihood->Draw("P");
    h_DISLikelihood->SetMarkerStyle(21);
    h_DISLikelihood->SetMarkerColor(sixColourScheme[1]);
    h_DISLikelihood->SetMarkerSize(1.75);
    h_DISLikelihood->Draw("P SAME");
    legend_muPi->Draw();
    c_likelihood->Update();

    TLegend* legend_PE = new TLegend(0.2, 0.2, 0.4, 0.4);
    legend_PE->AddEntry(g_purity, "Purity", "P");
    legend_PE->AddEntry(g_efficiency, "Efficiency", "P");
    legend_PE->SetBorderSize(0);
    legend_PE->SetFillStyle(0);
    legend_PE->SetTextSize(0.04);

    TCanvas *c_purityEfficiency = new TCanvas("c_purityEfficiency", "c_purityEfficiency", 800, 600);
    g_purity->SetMarkerStyle(20);
    g_purity->SetMarkerColor(sixColourScheme[2]);
    g_purity->SetMarkerSize(1.75);
    g_purity->GetYaxis()->SetRangeUser(0.0, 1.05);
    g_purity->GetYaxis()->SetTitle("");
    g_purity->Draw("AP");
    g_efficiency->SetMarkerStyle(20);
    g_efficiency->SetMarkerColor(sixColourScheme[3]);
    g_efficiency->SetMarkerSize(1.75);
    g_efficiency->Draw("P SAME");
    legend_PE->Draw();
    c_purityEfficiency->Update();
}