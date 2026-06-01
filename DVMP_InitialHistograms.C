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
#include <TPaveText.h>
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
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TLatex.h>
#include <Math/Boost.h>
#include <TLegend.h>
#include <TLine.h>

#include "ePICStyle.C"

void DVMP_InitialHistograms()
{
    //gROOT->SetBatch(kTRUE);
    gROOT->ProcessLine("SetePICStyle()");
    gStyle->SetOptStat(0);

    gStyle->SetCanvasPreferGL(true);


    std::vector<int> sixColourScheme = {
    TColor::GetColor("#7021dd"),     // violet
    TColor::GetColor("#964a8b"),     // grape
    TColor::GetColor("#e42536"),     // red
    TColor::GetColor("#f89c20"),     // yellow
    TColor::GetColor("#5790fc"),     // blue
    TColor::GetColor("#9c9ca1"),     // grey
    };

    TString infileN = "outputs/DVMP_JPsi_AnalysisOutput_10ifb_10x130ep.root";
    TFile *inFile = new TFile(infileN);

    std::string outFileDir = "outputs/InitialHistgrams/";
    std::string outfilename = outFileDir + "DVMP_JPsi_InitialHists_10ifb_10x130ep.root";

    // Set output file for the histograms
    TFile *ofile = TFile::Open(outfilename.c_str(),"RECREATE");

    // Load Eta and Momentum Histograms

    TH2D *electronEtaMom = (TH2D*) inFile->Get("electrons/electronEtaMom");
    TH2D *matchedElectronEtaMom = (TH2D*) inFile->Get("electrons/matchedElectronEtaMom");

    TH2D *protonEtaMom = (TH2D*) inFile->Get("protons/protonEtaMom");
    TH2D *matchedProtonEtaMom = (TH2D*) inFile->Get("protons/matchedProtonEtaMom");

    TH2D *muonEtaMom = (TH2D*) inFile->Get("muons/muonEtaMom");
    TH2D *matchedMuonEtaMom = (TH2D*) inFile->Get("muons/matchedMuonEtaMom");

    TH2D *JPsiEtaMom = (TH2D*) inFile->Get("JPsi/JPsiEtaMom");
    TH2D *matchedJPsiEtaMom = (TH2D*) inFile->Get("JPsi/matchedJPsiEtaMom");

    // Define canvases and legends for the 2D histograms

    TCanvas *C_electronEtaMom = new TCanvas("C_electronEtaMom", "Electron eta vs Momentum", 1900, 1000);
    electronEtaMom->GetXaxis()->SetRangeUser(-5, 1);
    electronEtaMom->GetYaxis()->SetRangeUser(-5, 20);
    electronEtaMom->Draw("COLZ");
    C_electronEtaMom->Update();

    TPaletteAxis *palette_electronEtaMom = (TPaletteAxis*)electronEtaMom->GetListOfFunctions()->FindObject("palette");

    if (palette_electronEtaMom) {
        // SetPosition(x1, y1, x2, y2)
        palette_electronEtaMom->SetX1NDC(0.840); // Move slightly right of the frame
        palette_electronEtaMom->SetX2NDC(0.875); // Width of the palette
        palette_electronEtaMom->SetY1NDC(0.155); // Align with bottom margin
        palette_electronEtaMom->SetY2NDC(0.955); // Align with top margin
        
        // Optional: Offset the axis labels so they don't overlap the box
        palette_electronEtaMom->SetLabelOffset(0.01);
        
        C_electronEtaMom->Modified();
        C_electronEtaMom->Update();
    }

    C_electronEtaMom->SaveAs("outputs/InitialHistgrams/10x130ep_electronEtaMom.pdf");

    TCanvas *C_matchedElectronEtaMom = new TCanvas("C_matchedElectronEtaMom", "Matched Electron eta vs Momentum", 1900, 1000);
    matchedElectronEtaMom->GetXaxis()->SetRangeUser(-5, 1);
    matchedElectronEtaMom->GetYaxis()->SetRangeUser(-5, 20);
    matchedElectronEtaMom->Draw("COLZ");
    C_matchedElectronEtaMom->Update();

    TPaletteAxis *palette_matchedElectronEtaMom = (TPaletteAxis*)matchedElectronEtaMom->GetListOfFunctions()->FindObject("palette");

    if (palette_matchedElectronEtaMom) {
        // SetPosition(x1, y1, x2, y2)
        palette_matchedElectronEtaMom->SetX1NDC(0.840); // Move slightly right of the frame
        palette_matchedElectronEtaMom->SetX2NDC(0.875); // Width of the palette
        palette_matchedElectronEtaMom->SetY1NDC(0.155); // Align with bottom margin
        palette_matchedElectronEtaMom->SetY2NDC(0.955); // Align with top margin
        
        // Optional: Offset the axis labels so they don't overlap the box
        palette_matchedElectronEtaMom->SetLabelOffset(0.01);
        
        C_matchedElectronEtaMom->Modified();
        C_matchedElectronEtaMom->Update();
    }

    C_matchedElectronEtaMom->SaveAs("outputs/InitialHistgrams/10x130ep_matchedElectronEtaMom.pdf");

    TCanvas *C_protonEtaMom = new TCanvas("C_protonEtaMom", "Proton eta vs Momentum", 1900, 1000);
    protonEtaMom->GetXaxis()->SetRangeUser(3, 10);
    protonEtaMom->GetYaxis()->SetRangeUser(100, 140);
    protonEtaMom->Draw("COLZ");
    C_protonEtaMom->SetLogz();
    C_protonEtaMom->Update();

    TPaletteAxis *palette_protonEtaMom = (TPaletteAxis*)protonEtaMom->GetListOfFunctions()->FindObject("palette");

    if (palette_protonEtaMom) {
        // SetPosition(x1, y1, x2, y2)
        palette_protonEtaMom->SetX1NDC(0.840); // Move slightly right of the frame
        palette_protonEtaMom->SetX2NDC(0.875); // Width of the palette
        palette_protonEtaMom->SetY1NDC(0.155); // Align with bottom margin
        palette_protonEtaMom->SetY2NDC(0.955); // Align with top margin
        
        // Optional: Offset the axis labels so they don't overlap the box
        palette_protonEtaMom->SetLabelOffset(0.01);
        
        C_protonEtaMom->Modified();
        C_protonEtaMom->Update();
    }

    C_protonEtaMom->SaveAs("outputs/InitialHistgrams/10x130ep_protonEtaMom.pdf");


    TCanvas *C_matchedProtonEtaMom = new TCanvas("C_matchedProtonEtaMom", "Matched Proton eta vs Momentum", 1900, 1000);
    matchedProtonEtaMom->GetXaxis()->SetRangeUser(3, 10);
    matchedProtonEtaMom->GetYaxis()->SetRangeUser(100, 140);
    matchedProtonEtaMom->Draw("COLZ");
    C_matchedProtonEtaMom->SetLogz();
    C_matchedProtonEtaMom->Update();

    TPaletteAxis *palette_matchedProtonEtaMom = (TPaletteAxis*)matchedProtonEtaMom->GetListOfFunctions()->FindObject("palette");

    if (palette_matchedProtonEtaMom) {
        // SetPosition(x1, y1, x2, y2)
        palette_matchedProtonEtaMom->SetX1NDC(0.840); // Move slightly right of the frame
        palette_matchedProtonEtaMom->SetX2NDC(0.875); // Width of the palette
        palette_matchedProtonEtaMom->SetY1NDC(0.155); // Align with bottom margin
        palette_matchedProtonEtaMom->SetY2NDC(0.955); // Align with top margin
        
        // Optional: Offset the axis labels so they don't overlap the box
        palette_matchedProtonEtaMom->SetLabelOffset(0.01);
        
        C_matchedProtonEtaMom->Modified();
        C_matchedProtonEtaMom->Update();
    }

    C_matchedProtonEtaMom->SaveAs("outputs/InitialHistgrams/10x130ep_matchedProtonEtaMom.pdf");

    TCanvas *C_muonEtaMom = new TCanvas("C_muonEtaMom", "Muon eta vs Momentum", 1900, 1000);
    muonEtaMom->GetXaxis()->SetRangeUser(-5, 5);
    muonEtaMom->GetYaxis()->SetRangeUser(0, 20);
    muonEtaMom->Draw("COLZ");
    C_muonEtaMom->Update();

    TPaletteAxis *palette_muonEtaMom = (TPaletteAxis*)muonEtaMom->GetListOfFunctions()->FindObject("palette");

    if (palette_muonEtaMom) {
        // SetPosition(x1, y1, x2, y2)
        palette_muonEtaMom->SetX1NDC(0.840); // Move slightly right of the frame
        palette_muonEtaMom->SetX2NDC(0.875); // Width of the palette
        palette_muonEtaMom->SetY1NDC(0.155); // Align with bottom margin
        palette_muonEtaMom->SetY2NDC(0.955); // Align with top margin
        
        // Optional: Offset the axis labels so they don't overlap the box
        palette_muonEtaMom->SetLabelOffset(0.01);
        
        C_muonEtaMom->Modified();
        C_muonEtaMom->Update();
    }

    C_muonEtaMom->SaveAs("outputs/InitialHistgrams/10x130ep_muonEtaMom.pdf");


    TCanvas *C_matchedMuonEtaMom = new TCanvas("C_matchedMuonEtaMom", "Matched Muon eta vs Momentum", 1900, 1000);
    matchedMuonEtaMom->GetXaxis()->SetRangeUser(-5,5);
    matchedMuonEtaMom->GetYaxis()->SetRangeUser(0, 20);
    matchedMuonEtaMom->Draw("COLZ");
    C_matchedMuonEtaMom->Update();

    TPaletteAxis *palette_matchedMuonEtaMom = (TPaletteAxis*)matchedMuonEtaMom->GetListOfFunctions()->FindObject("palette");

    if (palette_matchedMuonEtaMom) {
        // SetPosition(x1, y1, x2, y2)
        palette_matchedMuonEtaMom->SetX1NDC(0.840); // Move slightly right of the frame
        palette_matchedMuonEtaMom->SetX2NDC(0.875); // Width of the palette
        palette_matchedMuonEtaMom->SetY1NDC(0.155); // Align with bottom margin
        palette_matchedMuonEtaMom->SetY2NDC(0.955); // Align with top margin
        
        // Optional: Offset the axis labels so they don't overlap the box
        palette_matchedMuonEtaMom->SetLabelOffset(0.01);
        
        C_matchedMuonEtaMom->Modified();
        C_matchedMuonEtaMom->Update();
    }

    C_matchedMuonEtaMom->SaveAs("outputs/InitialHistgrams/10x130ep_matchedMuonEtaMom.pdf");

}