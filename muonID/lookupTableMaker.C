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

void lookupTableMaker()
{
    gROOT->SetBatch(kTRUE);

    TString infileN = "outputs/mu-pi_20-100GeV_epPlots.root";
    TFile *inFile = new TFile(infileN);

    std::string outfilename = "mu-pi_20-100GeV_Ep_lookupTable.root";

    // Ecal ep Plots by particle
    TH2D *pionEpTrueEcal = (TH2D*) inFile->Get("trueEcalEp/pionEpTrueEcal");
    TH2D *muonEpTrueEcal = (TH2D*) inFile->Get("trueEcalEp/muonEpTrueEcal");

    // Hcal ep Plots by particle
    TH2D *pionEpTrueHcal = (TH2D*) inFile->Get("trueHcalEp/pionEpTrueHcal");
    TH2D *muonEpTrueHcal = (TH2D*) inFile->Get("trueHcalEp/muonEpTrueHcal");


    TH1D *pionEcalProbabilities[pionEpTrueEcal->GetNbinsX()];
    TH1D *muonEcalProbabilities[muonEpTrueEcal->GetNbinsX()];

    TH1D *pionHcalProbabilities[pionEpTrueHcal->GetNbinsX()];
    TH1D *muonHcalProbabilities[muonEpTrueHcal->GetNbinsX()];


    std::vector<std::vector<double>> pionEcalLookupTable(pionEpTrueEcal->GetNbinsX() + 1, std::vector<double>(pionEpTrueEcal->GetNbinsY() + 1, 0.0));
    std::vector<std::vector<double>> muonEcalLookupTable(muonEpTrueEcal->GetNbinsX() + 1, std::vector<double>(muonEpTrueEcal->GetNbinsY() + 1, 0.0));
    std::vector<std::vector<double>> pionHcalLookupTable(pionEpTrueHcal->GetNbinsX() + 1, std::vector<double>(pionEpTrueHcal->GetNbinsY() + 1, 0.0));
    std::vector<std::vector<double>> muonHcalLookupTable(muonEpTrueHcal->GetNbinsX() + 1, std::vector<double>(muonEpTrueHcal->GetNbinsY() + 1, 0.0));

    std::vector<double> pValues(pionEpTrueEcal->GetNbinsX() + 1, 0.0);
    std::vector<double> EpValues(pionEpTrueEcal->GetNbinsY() + 1, 0.0);

    for (int iE = 1; iE < pionEpTrueEcal->GetNbinsX(); iE++) {
        pValues.at(iE) = pionEpTrueEcal->GetXaxis()->GetBinCenter(iE);
        pionEcalProbabilities[iE] = pionEpTrueEcal->ProjectionY(Form("pionEcalProbabilities_bin%d", iE), iE, iE);
        for (int jE = 0; jE < pionEpTrueEcal->GetNbinsY(); jE++) {
            double truep = pionEpTrueEcal->GetXaxis()->GetBinCenter(iE);
            double trueEp = pionEpTrueEcal->GetYaxis()->GetBinCenter(jE);
            double count = pionEpTrueEcal->GetBinContent(iE, jE);
            EpValues.at(jE) = trueEp;
            if (pionEcalProbabilities[iE]->Integral() == 0) continue;
            if (count > 0) {
                double probability = count / pionEcalProbabilities[iE]->Integral();
                int recBinX = pionEpTrueEcal->GetXaxis()->FindBin(truep);
                int recBinY = pionEpTrueEcal->GetYaxis()->FindBin(trueEp);
                if (recBinX > 0 && recBinY > 0) {
                    pionEcalLookupTable[recBinX][recBinY] = probability;
                }
            }
        }
    }

    for (int iE = 1; iE < muonEpTrueEcal->GetNbinsX(); iE++) {
        muonEcalProbabilities[iE] = muonEpTrueEcal->ProjectionY(Form("muonEcalProbabilities_bin%d", iE), iE, iE);
        for (int jE = 0; jE < muonEpTrueEcal->GetNbinsY(); jE++) {
            if (muonEcalProbabilities[iE]->Integral() == 0) continue;
            double truep = muonEpTrueEcal->GetXaxis()->GetBinCenter(iE);
            double trueEp = muonEpTrueEcal->GetYaxis()->GetBinCenter(jE);
            double count = muonEpTrueEcal->GetBinContent(iE, jE);
            if (count > 0) {
                double probability = count / muonEcalProbabilities[iE]->Integral();
                int recBinX = muonEpTrueEcal->GetXaxis()->FindBin(truep);
                int recBinY = muonEpTrueEcal->GetYaxis()->FindBin(trueEp);
                if (recBinX > 0 && recBinY > 0) {
                    muonEcalLookupTable[recBinX][recBinY] = probability;
                }
            }
        }
    }

    for (int iH = 1; iH < muonEpTrueHcal->GetNbinsX(); iH++) {
        pionHcalProbabilities[iH] = pionEpTrueHcal->ProjectionY(Form("pionHcalProbabilities_bin%d", iH), iH, iH);
        for (int jH = 0; jH < muonEpTrueHcal->GetNbinsY(); jH++) {
            if (pionHcalProbabilities[iH]->Integral() == 0) continue;
            double truep = pionEpTrueHcal->GetXaxis()->GetBinCenter(iH);
            double trueEp = pionEpTrueHcal->GetYaxis()->GetBinCenter(jH);
            double count = pionEpTrueHcal->GetBinContent(iH, jH);
            if (count > 0) {
                double probability = count / pionHcalProbabilities[iH]->Integral();
                int recBinX = pionEpTrueHcal->GetXaxis()->FindBin(truep);
                int recBinY = pionEpTrueHcal->GetYaxis()->FindBin(trueEp);
                if (recBinX > 0 && recBinY > 0) {
                    pionHcalLookupTable[recBinX][recBinY] = probability;
                }
            }
        }
    }

    //std::cout << "Number of bins in p: " << muonEpTrueHcal->GetNbinsX() << ", Number of bins in Ep: " << muonEpTrueHcal->GetNbinsY() << std::endl;

    for (int iH = 1; iH <= muonEpTrueHcal->GetNbinsX(); iH++) {
        //std::cout << "Processing muon Hcal bin " << iH + 1 << ", with momentum " << muonEpTrueHcal->GetXaxis()->GetBinCenter(iH) << std::endl;
        muonHcalProbabilities[iH] = muonEpTrueHcal->ProjectionY(Form("muonHcalProbabilities_bin%d", iH), iH, iH);
        for (int jH = 0; jH <= muonEpTrueHcal->GetNbinsY(); jH++) {
            if (muonHcalProbabilities[iH]->Integral() == 0) continue;
            double truep = muonEpTrueHcal->GetXaxis()->GetBinCenter(iH);
            double trueEp = muonEpTrueHcal->GetYaxis()->GetBinCenter(jH);
            double count = muonEpTrueHcal->GetBinContent(iH, jH);
            if (count > 0) {
                double probability = count / muonHcalProbabilities[iH]->Integral();
                int recBinX = muonEpTrueHcal->GetXaxis()->FindBin(truep);
                int recBinY = muonEpTrueHcal->GetYaxis()->FindBin(trueEp);
                if (recBinX > 0 && recBinY > 0) {
                    muonHcalLookupTable[recBinX][recBinY] = probability;
                    //std::cout << "Muon Hcal Lookup - True p: " << truep << ", True Ep: " << trueEp << ", Probability: " << probability << std::endl;
                }
            }
        }
    }
    

    // Set output file for the lookup tables
    TFile *ofile = TFile::Open(outfilename.c_str(),"RECREATE");

    double p_val;
    std::vector<double> ep_array;
    std::vector<double> prob_array;

    TTree *pionEcalLookupTree = new TTree("pionEcalLookupTree", "Pion EcalLookup Table Tree");
    pionEcalLookupTree->Branch("pionp", &p_val, "p_val/D");
    pionEcalLookupTree->Branch("trueEp", &ep_array);
    pionEcalLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        ep_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < EpValues.size(); j++) {
            ep_array.push_back(EpValues.at(j));
            prob_array.push_back(pionEcalLookupTable[i][j]);
        }
        pionEcalLookupTree->Fill();
    }

    pionEcalLookupTree->Write();


    TTree *muonEcalLookupTree = new TTree("muonEcalLookupTree", "Muon EcalLookup Table Tree");
    muonEcalLookupTree->Branch("muonp", &p_val, "p_val/D");
    muonEcalLookupTree->Branch("trueEp", &ep_array);
    muonEcalLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        ep_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < EpValues.size(); j++) {
            ep_array.push_back(EpValues.at(j));
            prob_array.push_back(muonEcalLookupTable[i][j]);
        }
        muonEcalLookupTree->Fill();
    }

    muonEcalLookupTree->Write();

    TTree *pionHcalLookupTree = new TTree("pionHcalLookupTree", "Pion HcalLookup Table Tree");
    pionHcalLookupTree->Branch("pionp", &p_val, "p_val/D");
    pionHcalLookupTree->Branch("trueEp", &ep_array);
    pionHcalLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        ep_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < EpValues.size(); j++) {
            ep_array.push_back(EpValues.at(j));
            prob_array.push_back(pionHcalLookupTable[i][j]);
        }
        pionHcalLookupTree->Fill();
    }

    pionHcalLookupTree->Write();


    TTree *muonHcalLookupTree = new TTree("muonHcalLookupTree", "Muon HcalLookup Table Tree");
    muonHcalLookupTree->Branch("muonp", &p_val, "p_val/D");
    muonHcalLookupTree->Branch("trueEp", &ep_array);
    muonHcalLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        ep_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < EpValues.size(); j++) {
            ep_array.push_back(EpValues.at(j));
            prob_array.push_back(muonHcalLookupTable[i][j]);
        }
        muonHcalLookupTree->Fill();
    }

    muonHcalLookupTree->Write();


    ofile->cd();
    ofile->Close();
}