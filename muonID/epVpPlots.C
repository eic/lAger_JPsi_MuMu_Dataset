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

void epVpPlots()
{
    //gROOT->SetBatch(kTRUE);
    //gROOT->ProcessLine("SetePICStyle()");
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

    TString infileN = "outputs/mu-pi_20-100GeV_epPlots.root";
    TFile *inFile = new TFile(infileN);

    std::string outfilename = "outputs/mu-pi_20-100GeV_EpVp.root";

    // Ecal ep Plots by particle
    TH2D *pionEpTrueEcal = (TH2D*) inFile->Get("trueEcalEp/pionEpTrueEcal");
    TH2D *muonEpTrueEcal = (TH2D*) inFile->Get("trueEcalEp/muonEpTrueEcal");

    TH2D *pionEpRecEcal = (TH2D*) inFile->Get("recoEcalEp/pionEpRecEcal");
    TH2D *muonEpRecEcal = (TH2D*) inFile->Get("recoEcalEp/muonEpRecEcal");

    // Hcal ep Plots by particle
    TH2D *pionEpTrueHcal = (TH2D*) inFile->Get("trueHcalEp/pionEpTrueHcal");
    TH2D *muonEpTrueHcal = (TH2D*) inFile->Get("trueHcalEp/muonEpTrueHcal");

    TH2D *pionEpRecHcal = (TH2D*) inFile->Get("recoHcalEp/pionEpRecHcal");
    TH2D *muonEpRecHcal = (TH2D*) inFile->Get("recoHcalEp/muonEpRecHcal");


    TH1D *pionEcalProbabilities[pionEpRecEcal->GetNbinsX()];
    TH1D *muonEcalProbabilities[muonEpRecEcal->GetNbinsX()];

    TH1D *pionHcalProbabilities[pionEpRecHcal->GetNbinsX()];
    TH1D *muonHcalProbabilities[muonEpRecHcal->GetNbinsX()];
    

    for (int iE = 1; iE <= pionEpRecEcal->GetNbinsX(); iE++) {
        pionEcalProbabilities[iE-1] = pionEpRecEcal->ProjectionY(Form("pionEcalProbabilities_bin%d", iE), iE, iE);
        muonEcalProbabilities[iE-1] = muonEpRecEcal->ProjectionY(Form("muonEcalProbabilities_bin%d", iE), iE, iE);
        
        // CRITICAL: Detach from the input file so they stay in memory
        pionEcalProbabilities[iE-1]->SetDirectory(0);
        muonEcalProbabilities[iE-1]->SetDirectory(0);

        if (pionEcalProbabilities[iE-1]->Integral() != 0) pionEcalProbabilities[iE-1]->Scale(1.0/pionEcalProbabilities[iE-1]->Integral());
        if (muonEcalProbabilities[iE-1]->Integral() != 0) muonEcalProbabilities[iE-1]->Scale(1.0/muonEcalProbabilities[iE-1]->Integral());
    }

    for (int iH = 1; iH <= pionEpRecHcal->GetNbinsX(); iH++) {
        pionHcalProbabilities[iH-1] = pionEpRecHcal->ProjectionY(Form("pionHcalProbabilities_bin%d", iH), iH, iH);
        muonHcalProbabilities[iH-1] = muonEpRecHcal->ProjectionY(Form("muonHcalProbabilities_bin%d", iH), iH, iH);

        // CRITICAL: Detach from the input file so they stay in memory
        pionHcalProbabilities[iH-1]->SetDirectory(0);
        muonHcalProbabilities[iH-1]->SetDirectory(0);

        if (pionHcalProbabilities[iH-1]->Integral() != 0) pionHcalProbabilities[iH-1]->Scale(1.0/pionHcalProbabilities[iH-1]->Integral());
        if (muonHcalProbabilities[iH-1]->Integral() != 0) muonHcalProbabilities[iH-1]->Scale(1.0/muonHcalProbabilities[iH-1]->Integral());
    }

    TCanvas *c_EcalProbabilities = new TCanvas("c_EcalProbabilities", "Ecal Probabilities", 1200, 800);
    c_EcalProbabilities->Divide(pionEpRecEcal->GetNbinsX()/4, 4);
    for (int i = 0; i < pionEpRecEcal->GetNbinsX(); i++) {
        c_EcalProbabilities->cd(i+1);
        pionEcalProbabilities[i]->SetLineColor(sixColourScheme[0]);
        muonEcalProbabilities[i]->SetLineColor(sixColourScheme[1]);
        pionEcalProbabilities[i]->SetMarkerColor(sixColourScheme[0]);
        muonEcalProbabilities[i]->SetMarkerColor(sixColourScheme[1]);
        pionEcalProbabilities[i]->SetMarkerStyle(20);
        muonEcalProbabilities[i]->SetMarkerStyle(21);
        pionEcalProbabilities[i]->Draw("HIST P");
        muonEcalProbabilities[i]->Draw("HIST P SAME");
    }
    c_EcalProbabilities->Update();

    TCanvas *c_HcalProbabilities = new TCanvas("c_HcalProbabilities", "Hcal Probabilities", 1200, 800);
    c_HcalProbabilities->Divide(pionEpRecHcal->GetNbinsX()/4, 4);
    for (int i = 0; i < pionEpRecHcal->GetNbinsX(); i++) {
        c_HcalProbabilities->cd(i+1);
        pionHcalProbabilities[i]->SetLineColor(sixColourScheme[0]);
        muonHcalProbabilities[i]->SetLineColor(sixColourScheme[1]);
        pionHcalProbabilities[i]->SetMarkerColor(sixColourScheme[0]);
        muonHcalProbabilities[i]->SetMarkerColor(sixColourScheme[1]);
        pionHcalProbabilities[i]->SetMarkerStyle(20);
        muonHcalProbabilities[i]->SetMarkerStyle(21);
        pionHcalProbabilities[i]->Draw("HIST P");
        muonHcalProbabilities[i]->Draw("HIST P SAME");
    }
    c_HcalProbabilities->Update();


    // Set output file for the histograms
    TFile *ofile = TFile::Open(outfilename.c_str(),"RECREATE");


    ofile->cd();
    ofile->mkdir("eCalProbabilities");
    ofile->cd("eCalProbabilities");
    for (int i = 0; i < pionEpRecEcal->GetNbinsX(); i++) {
        pionEcalProbabilities[i]->Write();
        muonEcalProbabilities[i]->Write();
    }
    ofile->cd("..");    
    ofile->mkdir("hCalProbabilities");
    ofile->cd("hCalProbabilities");
    for (int i = 0; i < pionEpRecHcal->GetNbinsX(); i++) {
        pionHcalProbabilities[i]->Write();
        muonHcalProbabilities[i]->Write();  
    }
    ofile->cd("..");
    ofile->Close();
}