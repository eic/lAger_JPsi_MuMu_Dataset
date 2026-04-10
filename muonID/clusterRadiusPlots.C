#include <iomanip>
#include <cstdlib>
#include <cmath>
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
#include <TLegend.h>
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
#include <TColor.h>

#include "/shared/physics/physdata/nuclear/ePIC_EIC/gbx505/ePIC_JPsi/ePICStyle.C"

void clusterRadiusPlots()
{
    //gROOT->SetBatch(kTRUE);
    gROOT->ProcessLine("SetePICStyle()");

    // New Colour Scheme palette
    std::vector<int> sixColourScheme = {
        TColor::GetColor("#7021dd"),     // violet
        TColor::GetColor("#964a8b"),     // grape
        TColor::GetColor("#e42536"),     // red
        TColor::GetColor("#f89c20"),     // yellow
        TColor::GetColor("#5790fc"),     // blue
        TColor::GetColor("#9c9ca1"),     // grey
    };

    std::string fileNameList[2] = {
                                  "outputs/muCalPlotsOutput.root",
                                  "outputs/piCalPlotsOutput.root",                                  
                                };

    std::string ParticleTypeList[2] = {"muons", "pions"};
    int numberOfFiles = 2;

    TH1D *h_EcalBarrelImagingHitsRadius[2];
    TH1D *h_EcalBarrelScFiHitsRadius[2];
    TH1D *h_EcalEndcapPHitsRadius[2];
    TH1D *h_EcalEndcapNHitsRadius[2];
    TH1D *h_HcalBarrelHitsRadius[2];
    TH1D *h_HcalEndcapPHitsRadius[2];
    TH1D *h_HcalEndcapNHitsRadius[2];
    TH1D *h_HcalLFHitsRadius[2];

    TLegend* legend = new TLegend(0.55, 0.7, 0.8, 0.9);
    legend->SetHeader("Particle", "C");

    for (int fileNumber= 0; fileNumber< numberOfFiles; fileNumber++){
        TFile *infile = TFile::Open(fileNameList[fileNumber].c_str(),"READ");
        if (!infile || infile->IsZombie()) {
            std::cerr << "Error opening file: " << fileNameList[fileNumber] << std::endl;
            continue;
        }

        h_EcalBarrelImagingHitsRadius[fileNumber] = (TH1D*)infile->Get("radiusPlots/h_EcalBarrelImagingHitsRadius");
        h_EcalBarrelScFiHitsRadius[fileNumber] = (TH1D*)infile->Get("radiusPlots/h_EcalBarrelScFiHitsRadius");
        h_EcalEndcapPHitsRadius[fileNumber] = (TH1D*)infile->Get("radiusPlots/h_EcalEndcapPHitsRadius");
        h_EcalEndcapNHitsRadius[fileNumber] = (TH1D*)infile->Get("radiusPlots/h_EcalEndcapNHitsRadius");
        h_HcalBarrelHitsRadius[fileNumber] = (TH1D*)infile->Get("radiusPlots/h_HcalBarrelHitsRadius");
        h_HcalEndcapPHitsRadius[fileNumber] = (TH1D*)infile->Get("radiusPlots/h_HcalEndcapPHitsRadius");
        h_HcalEndcapNHitsRadius[fileNumber] = (TH1D*)infile->Get("radiusPlots/h_HcalEndcapNHitsRadius");
        h_HcalLFHitsRadius[fileNumber] = (TH1D*)infile->Get("radiusPlots/h_HcalLFHitsRadius");
        legend->AddEntry(h_EcalBarrelImagingHitsRadius[fileNumber], ParticleTypeList[fileNumber].c_str(), "p");

    }

    TCanvas *c_EcalBarrelImagingHitsRadius = new TCanvas("c_EcalBarrelImagingHitsRadius", "Ecal Barrel Imaging Hits Radius", 800, 600);
    h_EcalBarrelImagingHitsRadius[0]->SetLineColor(sixColourScheme[2]);
    h_EcalBarrelImagingHitsRadius[0]->SetMarkerColor(sixColourScheme[2]);
    h_EcalBarrelImagingHitsRadius[0]->SetMarkerStyle(20);
    h_EcalBarrelImagingHitsRadius[0]->Draw("HIST PE");
    h_EcalBarrelImagingHitsRadius[1]->SetLineColor(sixColourScheme[4]);
    h_EcalBarrelImagingHitsRadius[1]->SetMarkerColor(sixColourScheme[4]);
    h_EcalBarrelImagingHitsRadius[1]->SetMarkerStyle(22);
    h_EcalBarrelImagingHitsRadius[1]->Draw("HIST PE SAME");
    legend->Draw();
    c_EcalBarrelImagingHitsRadius->SetLogy();
    c_EcalBarrelImagingHitsRadius->Draw();
    
    TCanvas *c_EcalBarrelScFiHitsRadius = new TCanvas("c_EcalBarrelScFiHitsRadius", "Ecal Barrel ScFi Hits Radius", 800, 600);
    h_EcalBarrelScFiHitsRadius[0]->SetLineColor(sixColourScheme[2]);
    h_EcalBarrelScFiHitsRadius[0]->SetMarkerColor(sixColourScheme[2]);
    h_EcalBarrelScFiHitsRadius[0]->SetMarkerStyle(20);
    h_EcalBarrelScFiHitsRadius[0]->Draw("HIST PE");
    h_EcalBarrelScFiHitsRadius[1]->SetLineColor(sixColourScheme[4]);
    h_EcalBarrelScFiHitsRadius[1]->SetMarkerColor(sixColourScheme[4]);
    h_EcalBarrelScFiHitsRadius[1]->SetMarkerStyle(22);
    h_EcalBarrelScFiHitsRadius[1]->Draw("HIST PE SAME");
    legend->Draw();
    c_EcalBarrelScFiHitsRadius->SetLogy();
    c_EcalBarrelScFiHitsRadius->Draw();


    TCanvas *c_EcalEndcapPHitsRadius = new TCanvas("c_EcalEndcapPHitsRadius", "Ecal Endcap P Hits Radius", 800, 600);
    h_EcalEndcapPHitsRadius[0]->SetLineColor(sixColourScheme[2]);
    h_EcalEndcapPHitsRadius[0]->SetMarkerColor(sixColourScheme[2]);
    h_EcalEndcapPHitsRadius[0]->SetMarkerStyle(20);
    h_EcalEndcapPHitsRadius[0]->Draw("HIST PE");
    h_EcalEndcapPHitsRadius[1]->SetLineColor(sixColourScheme[4]);
    h_EcalEndcapPHitsRadius[1]->SetMarkerColor(sixColourScheme[4]);
    h_EcalEndcapPHitsRadius[1]->SetMarkerStyle(22);
    h_EcalEndcapPHitsRadius[1]->Draw("HIST PE SAME");
    legend->Draw();
    c_EcalEndcapPHitsRadius->SetLogy();
    c_EcalEndcapPHitsRadius->Draw();

    TCanvas *c_EcalEndcapNHitsRadius = new TCanvas("c_EcalEndcapNHitsRadius", "Ecal Endcap N Hits Radius", 800, 600);
    h_EcalEndcapNHitsRadius[0]->SetLineColor(sixColourScheme[2]);
    h_EcalEndcapNHitsRadius[0]->SetMarkerColor(sixColourScheme[2]);
    h_EcalEndcapNHitsRadius[0]->SetMarkerStyle(20);
    h_EcalEndcapNHitsRadius[0]->Draw("HIST PE");
    h_EcalEndcapNHitsRadius[1]->SetLineColor(sixColourScheme[4]);
    h_EcalEndcapNHitsRadius[1]->SetMarkerColor(sixColourScheme[4]);
    h_EcalEndcapNHitsRadius[1]->SetMarkerStyle(22);
    h_EcalEndcapNHitsRadius[1]->Draw("HIST PE SAME");
    legend->Draw();
    c_EcalEndcapNHitsRadius->SetLogy();
    c_EcalEndcapNHitsRadius->Draw();

    TCanvas *c_HcalBarrelHitsRadius = new TCanvas("c_HcalBarrelHitsRadius", "Hcal Barrel Hits Radius", 800, 600);
    h_HcalBarrelHitsRadius[0]->SetLineColor(sixColourScheme[2]);
    h_HcalBarrelHitsRadius[0]->SetMarkerColor(sixColourScheme[2]);
    h_HcalBarrelHitsRadius[0]->SetMarkerStyle(20);
    h_HcalBarrelHitsRadius[0]->Draw("HIST PE");
    h_HcalBarrelHitsRadius[1]->SetLineColor(sixColourScheme[4]);
    h_HcalBarrelHitsRadius[1]->SetMarkerColor(sixColourScheme[4]);
    h_HcalBarrelHitsRadius[1]->SetMarkerStyle(22);
    h_HcalBarrelHitsRadius[1]->Draw("HIST PE SAME");
    legend->Draw();
    c_HcalBarrelHitsRadius->SetLogy();
    c_HcalBarrelHitsRadius->Draw();

    TCanvas *c_HcalEndcapPHitsRadius = new TCanvas("c_HcalEndcapPHitsRadius", "Hcal Endcap P Hits Radius", 800, 600);
    h_HcalEndcapPHitsRadius[1]->SetLineColor(sixColourScheme[4]);
    h_HcalEndcapPHitsRadius[1]->SetMarkerColor(sixColourScheme[4]);
    h_HcalEndcapPHitsRadius[1]->SetMarkerStyle(22);
    h_HcalEndcapPHitsRadius[1]->Draw("HIST PE");
    h_HcalEndcapPHitsRadius[0]->SetLineColor(sixColourScheme[2]);
    h_HcalEndcapPHitsRadius[0]->SetMarkerColor(sixColourScheme[2]);
    h_HcalEndcapPHitsRadius[0]->SetMarkerStyle(20);
    h_HcalEndcapPHitsRadius[0]->Draw("HIST PE SAME");  
    legend->Draw();
    c_HcalEndcapPHitsRadius->SetLogy();
    c_HcalEndcapPHitsRadius->Draw();

    TCanvas *c_HcalEndcapNHitsRadius = new TCanvas("c_HcalEndcapNHitsRadius", "Hcal Endcap N Hits Radius", 800, 600);
    h_HcalEndcapNHitsRadius[0]->SetLineColor(sixColourScheme[2]);
    h_HcalEndcapNHitsRadius[0]->SetMarkerColor(sixColourScheme[2]);
    h_HcalEndcapNHitsRadius[0]->SetMarkerStyle(20);
    h_HcalEndcapNHitsRadius[0]->Draw("HIST PE");
    h_HcalEndcapNHitsRadius[1]->SetLineColor(sixColourScheme[4]);
    h_HcalEndcapNHitsRadius[1]->SetMarkerColor(sixColourScheme[4]);
    h_HcalEndcapNHitsRadius[1]->SetMarkerStyle(22);
    h_HcalEndcapNHitsRadius[1]->Draw("HIST PE SAME");
    legend->Draw();
    c_HcalEndcapNHitsRadius->SetLogy();
    c_HcalEndcapNHitsRadius->Draw();

    TCanvas *c_HcalLFHitsRadius = new TCanvas("c_LFHcalHitsRadius", "LFHcal Hits Radius", 800, 600);
    h_HcalLFHitsRadius[0]->SetLineColor(sixColourScheme[2]);
    h_HcalLFHitsRadius[0]->SetMarkerColor(sixColourScheme[2]);
    h_HcalLFHitsRadius[0]->SetMarkerStyle(20);
    h_HcalLFHitsRadius[0]->Draw("HIST PE");
    h_HcalLFHitsRadius[1]->SetLineColor(sixColourScheme[4]);
    h_HcalLFHitsRadius[1]->SetMarkerColor(sixColourScheme[4]);
    h_HcalLFHitsRadius[1]->SetMarkerStyle(22);
    h_HcalLFHitsRadius[1]->Draw("HIST PE SAME");
    legend->Draw();
    c_HcalLFHitsRadius->SetLogy();
    c_HcalLFHitsRadius->Draw();
}