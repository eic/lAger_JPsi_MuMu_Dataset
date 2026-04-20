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

void likelihood()
{
    //gROOT->SetBatch(kTRUE);
    //gROOT->ProcessLine("SetePICStyle()");
    gStyle->SetOptStat(0);
    gStyle->SetCanvasPreferGL(true); // enables transparency

    std::vector<int> sixColourScheme = {
    TColor::GetColor("#5790fc"),     // blue
    TColor::GetColor("#e42536"),     // red
    TColor::GetColor("#964a8b"),     // grape
    TColor::GetColor("#f89c20"),     // yellow
    TColor::GetColor("#7021dd"),     // violet
    TColor::GetColor("#9c9ca1"),     // grey
    };

    TString lookupN = "mu-pi_20-100GeV_Ep_lookupTable.root";
    TFile *lookupFile = new TFile(lookupN);

    //TString infile="../eicReconOutput/SimCampaign_JPsiMuMu_10ifb_10x130ep_Pruned.root";
    TString infile="reconOut/pi_100GeV_reconOut.root";
    //TString infile="../dis_background/DIS_Q2_1_10_10x130ep_Pruned.root";

    double EndcapNHcal_Factor = 1.0/6.0;

    std::string outfilename = "outputs/pi_20-100GeV_likelihoodPlots.root";
    //std::string outfilename = "outputs/SimCampaign_JPsiMuMu_10ifb_10x130ep_likelihoodPlots.root";
    //std::string outfilename = "outputs/DIS_Q2_1_10_10x130ep_likelihoodPlots.root";

    TTree *muonEcalLookupTree = (TTree*) lookupFile->Get("muonEcalLookupTree");

    std::vector<double> *Ep_E_ptr = nullptr;
    std::vector<double> *prob_E_ptr = nullptr;
    double lookup_E_p;

    muonEcalLookupTree->SetBranchAddress("trueEp", &Ep_E_ptr);
    muonEcalLookupTree->SetBranchAddress("probability", &prob_E_ptr);
    muonEcalLookupTree->SetBranchAddress("muonp", &lookup_E_p);

    TTree *muonHcalLookupTree = (TTree*) lookupFile->Get("muonHcalLookupTree");

    std::vector<double> *Ep_H_ptr = nullptr;
    std::vector<double> *prob_H_ptr = nullptr;
    double lookup_H_p;

    muonHcalLookupTree->SetBranchAddress("trueEp", &Ep_H_ptr);
    muonHcalLookupTree->SetBranchAddress("probability", &prob_H_ptr);
    muonHcalLookupTree->SetBranchAddress("muonp", &lookup_H_p);

    double pMin = muonEcalLookupTree->GetMinimum("muonp");
    double pMax = muonEcalLookupTree->GetMaximum("muonp");
    int nRows = muonEcalLookupTree->GetEntries();

    // Set up input file chain
    TChain *mychain = new TChain("events");
    mychain->Add(infile);

    // Initialize reader
    TTreeReader particleSource_Reader(mychain);

    // Get Particle Information
    TTreeReaderArray<int> partGenStat(particleSource_Reader, "MCParticles.generatorStatus");
    TTreeReaderArray<double> partMomX(particleSource_Reader, "MCParticles.momentum.x");
    TTreeReaderArray<double> partMomY(particleSource_Reader, "MCParticles.momentum.y");
    TTreeReaderArray<double> partMomZ(particleSource_Reader, "MCParticles.momentum.z");
    TTreeReaderArray<int> partPdg(particleSource_Reader, "MCParticles.PDG");
    TTreeReaderArray<double> partMass(particleSource_Reader, "MCParticles.mass");
    TTreeReaderArray<float> partCharge(particleSource_Reader, "MCParticles.charge");
    TTreeReaderArray<unsigned int> partParb(particleSource_Reader, "MCParticles.parents_begin");
    TTreeReaderArray<unsigned int> partPare(particleSource_Reader, "MCParticles.parents_end");
    TTreeReaderArray<int> partParI(particleSource_Reader, "_MCParticles_parents.index");

    // Get Reconstructed Track Information
    TTreeReaderArray<float> recoMomX(particleSource_Reader, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float> recoMomY(particleSource_Reader, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float> recoMomZ(particleSource_Reader, "ReconstructedChargedParticles.momentum.z");
    TTreeReaderArray<int> trackPDG(particleSource_Reader, "ReconstructedChargedParticles.PDG");
    TTreeReaderArray<float> trackMass(particleSource_Reader, "ReconstructedChargedParticles.mass");
    TTreeReaderArray<float> trackCharge(particleSource_Reader, "ReconstructedChargedParticles.charge");
    TTreeReaderArray<float> trackEng(particleSource_Reader, "ReconstructedChargedParticles.energy");

    // Get Truth Seeded Reconstructed Track Information
    TTreeReaderArray<float> trecoMomX(particleSource_Reader, "ReconstructedTruthSeededChargedParticles.momentum.x");
    TTreeReaderArray<float> trecoMomY(particleSource_Reader, "ReconstructedTruthSeededChargedParticles.momentum.y");
    TTreeReaderArray<float> trecoMomZ(particleSource_Reader, "ReconstructedTruthSeededChargedParticles.momentum.z");
    TTreeReaderArray<int> ttrackPDG(particleSource_Reader, "ReconstructedTruthSeededChargedParticles.PDG");
    TTreeReaderArray<float> ttrackMass(particleSource_Reader, "ReconstructedTruthSeededChargedParticles.mass");
    TTreeReaderArray<float> ttrackCharge(particleSource_Reader, "ReconstructedTruthSeededChargedParticles.charge");
    TTreeReaderArray<float> ttrackEng(particleSource_Reader, "ReconstructedTruthSeededChargedParticles.energy");

    // Get Associations Between MCParticles and ReconstructedChargedParticles
    TTreeReaderArray<int> recoAssoc(particleSource_Reader, "_ReconstructedChargedParticleAssociations_rec.index");
    TTreeReaderArray<int> simuAssoc(particleSource_Reader, "_ReconstructedChargedParticleAssociations_sim.index");

    // Get B0 Information
    TTreeReaderArray<int> recoAssocB0(particleSource_Reader, "_B0ECalClusterAssociations_rec.index");
    TTreeReaderArray<int> simuAssocB0(particleSource_Reader, "_B0ECalClusterAssociations_sim.index");
    TTreeReaderArray<float> B0Eng(particleSource_Reader, "B0ECalClusters.energy");
    TTreeReaderArray<float> B0z(particleSource_Reader, "B0ECalClusters.position.z");

    // Get Forward Detector Information
    TTreeReaderArray<float> RPEng(particleSource_Reader, "ForwardRomanPotRecParticles.energy");
    TTreeReaderArray<int> RPpdg(particleSource_Reader, "ForwardRomanPotRecParticles.PDG");
    TTreeReaderArray<float> RPMomX(particleSource_Reader, "ForwardRomanPotRecParticles.momentum.x");
    TTreeReaderArray<float> RPMomY(particleSource_Reader, "ForwardRomanPotRecParticles.momentum.y");
    TTreeReaderArray<float> RPMomZ(particleSource_Reader, "ForwardRomanPotRecParticles.momentum.z");

    TTreeReaderArray<float> OffMEng(particleSource_Reader, "ForwardOffMRecParticles.energy");

    // Ecal Information
    TTreeReaderArray<int> simuAssocEcalBarrel(particleSource_Reader, "_EcalBarrelClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocEcalBarrel(particleSource_Reader, "_EcalBarrelClusterAssociations_rec.index");
    TTreeReaderArray<float> EcalBarrelEng(particleSource_Reader, "EcalBarrelClusters.energy");
    /*
    TTreeReaderArray<int> simuAssocEcalBarrelImg(particleSource_Reader, "EcalBarrelImagingClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocEcalBarrelImg(particleSource_Reader, "EcalBarrelImagingClusterAssociations_rec.index");
    TTreeReaderArray<float> EcalBarrelImgEng(particleSource_Reader, "EcalBarrelImagingClusters.energy");
    TTreeReaderArray<int> simuAssocEcalBarrelScFi(particleSource_Reader, "EcalBarrelScFiClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocEcalBarrelScFi(particleSource_Reader, "EcalBarrelScFiClusterAssociations_rec.index");
    TTreeReaderArray<float> EcalBarrelScFiEng(particleSource_Reader, "EcalBarrelScFiClusters.energy");
    */

    TTreeReaderArray<int> simuAssocEcalEndcapP(particleSource_Reader, "_EcalEndcapPClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocEcalEndcapP(particleSource_Reader, "_EcalEndcapPClusterAssociations_rec.index");    
    TTreeReaderArray<float> EcalEndcapPEng(particleSource_Reader, "EcalEndcapPClusters.energy");

    TTreeReaderArray<int> simuAssocEcalEndcapN(particleSource_Reader, "_EcalEndcapNClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocEcalEndcapN(particleSource_Reader, "_EcalEndcapNClusterAssociations_rec.index");
    TTreeReaderArray<float> EcalEndcapNEng(particleSource_Reader, "EcalEndcapNClusters.energy");

    // Hcal Information
    TTreeReaderArray<int> simuAssocHcalBarrel(particleSource_Reader, "_HcalBarrelClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocHcalBarrel(particleSource_Reader, "_HcalBarrelClusterAssociations_rec.index");
    TTreeReaderArray<float> HcalBarrelEng(particleSource_Reader, "HcalBarrelClusters.energy");

    TTreeReaderArray<int> simuAssocHcalEndcapP(particleSource_Reader, "_HcalEndcapPInsertClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocHcalEndcapP(particleSource_Reader, "_HcalEndcapPInsertClusterAssociations_rec.index");    
    TTreeReaderArray<float> HcalEndcapPEng(particleSource_Reader, "HcalEndcapPInsertClusters.energy");

    TTreeReaderArray<int> simuAssocLFHcal(particleSource_Reader, "_LFHCALClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocLFHcal(particleSource_Reader, "_LFHCALClusterAssociations_rec.index");    
    TTreeReaderArray<float> LFHcalEng(particleSource_Reader, "LFHCALClusters.energy");

    TTreeReaderArray<int> simuAssocHcalEndcapN(particleSource_Reader, "_HcalEndcapNClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocHcalEndcapN(particleSource_Reader, "_HcalEndcapNClusterAssociations_rec.index");
    TTreeReaderArray<float> HcalEndcapNEng(particleSource_Reader, "HcalEndcapNClusters.energy");

    // Histograms

    TH1D *h_muonLikelihood = new TH1D("h_muonLikelihood", "Muon Likelihood;Likelihood;Entries", 100.0, -60.0, 10.0);
    TH2D *h_muonLikelihoodVsp = new TH2D("h_muonLikelihoodVsp", "Muon Likelihood vs p;Momentum (GeV);Likelihood", 500, 0.0, 100.0, 100.0, -60.0, 10.0);

    TH1D *h_nonMuonLikelihood = new TH1D("h_nonMuonLikelihood", "Non-Muon Likelihood;Likelihood;Entries", 100.0, -60.0, 10.0);
    TH2D *h_nonMuonLikelihoodVsp = new TH2D("h_nonMuonLikelihoodVsp", "Non-Muon Likelihood vs p;Momentum (GeV);Likelihood", 500, 0.0, 100.0, 100.0, -60.0, 10.0);

    int eventID = 0;

    double ECalEnergy = 0.;
    double HCalEnergy = 0.;
    double Ep_Ecal = 0.;
    double Ep_Hcal = 0.;

    double particleMomentum = 0.;
    double particleEta = 0.;
    double truePDG = 0.;

    while(particleSource_Reader.Next()) // Loop over events
    {
        eventID++;

        if (eventID % 10 == 0)
        {
            fprintf (stderr, "%4.2f Percent\r ", eventID*100.0/mychain->GetEntries());
            fflush (stderr);
        }

        if (eventID > 10000) break; // Only run over first 10k events for testing purposes

        for (int i = 0; i < trackPDG.GetSize(); i++)
        {
            if (i >= simuAssoc.GetSize() || simuAssoc[i] >= partPdg.GetSize()) continue;

            ECalEnergy = 0.;
            HCalEnergy = 0.;

            TVector3 recoMom(recoMomX[i],recoMomY[i],recoMomZ[i]);
            ROOT::Math::PxPyPzEVector reco4Mom(recoMomX[i],recoMomY[i],recoMomZ[i], trackEng[i]);

            particleMomentum = recoMom.Mag();
            particleEta = recoMom.Eta();
            truePDG = TMath::Abs(partPdg[simuAssoc[i]]);

            if (TMath::Abs(particleEta) > 4) continue;
            if (TMath::Abs(particleEta) >= 1.0 && TMath::Abs(particleEta) <= 1.3) continue; // Remove overlap region between barrel and endcap calorimeters
            if (particleMomentum < 0.05 || particleMomentum > 100.0) continue;


            if (EcalBarrelEng.GetSize() == simuAssocEcalBarrel.GetSize())
            {
                for (int jB = 0; jB < simuAssocEcalBarrel.GetSize(); jB++) // Look for associations in the Ecal Barrel
                {
                    if (simuAssocEcalBarrel[jB] == simuAssoc[i])
                    {
                        ECalEnergy += EcalBarrelEng[jB];
                        
                    }
                }
            }
            /*
            if (EcalBarrelScFiEng.GetSize() == simuAssocEcalBarrelScFi.GetSize() && ECalEnergy == 0.0)
            {
                for (int jS = 0; jS < simuAssocEcalBarrelScFi.GetSize(); jS++) // Look for associations in the Ecal Barrel ScFi
                {
                    if (simuAssocEcalBarrelScFi[jS] == simuAssoc[i])
                    {
                        ECalEnergy += EcalBarrelScFiEng[jS];
                        
                    }
                }
            }
            if (EcalBarrelImgEng.GetSize() == simuAssocEcalBarrelImg.GetSize() && ECalEnergy == 0.0)
            {
                for (int jI = 0; jI < simuAssocEcalBarrelImg.GetSize(); jI++) // Look for associations in the Ecal Barrel Preshower
                {
                    if (simuAssocEcalBarrelImg[jI] == simuAssoc[i])
                    {
                        ECalEnergy += EcalBarrelImgEng[jI];
                        
                    }
                }
            }
            */
            if (EcalEndcapPEng.GetSize() == simuAssocEcalEndcapP.GetSize()) 
            {
                for (int jP = 0; jP < simuAssocEcalEndcapP.GetSize(); jP++) // Look for associations in the Ecal Endcap P
                {
                    if (simuAssocEcalEndcapP[jP] == simuAssoc[i])
                    {
                        ECalEnergy += EcalEndcapPEng[jP];
                        
                    }
                }
            }
            if (EcalEndcapNEng.GetSize() == simuAssocEcalEndcapN.GetSize())
            {
                for (int jN = 0; jN < simuAssocEcalEndcapN.GetSize(); jN++) // Look for associations in the Ecal Endcap N
                {
                    if (simuAssocEcalEndcapN[jN] == simuAssoc[i])
                    {
                        ECalEnergy += EcalEndcapNEng[jN];
                        
                    }
                }
            }   

            if (HcalBarrelEng.GetSize() == simuAssocHcalBarrel.GetSize())
            {
                for (int jB = 0; jB < simuAssocHcalBarrel.GetSize(); jB++) // Look for associations in the Ecal Barrel
                {
                    if (simuAssocHcalBarrel[jB] == simuAssoc[i])
                    {
                        HCalEnergy += HcalBarrelEng[jB];
                        
                    }
                }
            }
            
            if (HcalEndcapPEng.GetSize() == simuAssocHcalEndcapP.GetSize())
            {
                for (int jP = 0; jP < simuAssocHcalEndcapP.GetSize(); jP++) // Look for associations in the Ecal Endcap P
                {
                    if (simuAssocHcalEndcapP[jP] == simuAssoc[i])
                    {
                        HCalEnergy += HcalEndcapPEng[jP];
                        
                    }
                }
            }

            if (LFHcalEng.GetSize() == simuAssocLFHcal.GetSize())
            {
                for (int jL = 0; jL < simuAssocLFHcal.GetSize(); jL++) // Look for associations in the LFHcal
                {
                    if (simuAssocLFHcal[jL] == simuAssoc[i])
                    {
                        HCalEnergy += LFHcalEng[jL];
                        
                    }
                }
            }
            
            if (HcalEndcapNEng.GetSize() == simuAssocHcalEndcapN.GetSize())
            {
                for (int jN = 0; jN < simuAssocHcalEndcapN.GetSize(); jN++) // Look for associations in the Ecal Endcap N
                {
                    if (simuAssocHcalEndcapN[jN] == simuAssoc[i])
                    {
                        HCalEnergy += HcalEndcapNEng[jN] * EndcapNHcal_Factor;
                        
                    }
                }
            }

            if (ECalEnergy == 0.0 && HCalEnergy == 0.0) continue; // Skip if no calorimeter energy is associated

            Ep_Ecal = ECalEnergy / particleMomentum;
            Ep_Hcal = HCalEnergy / particleMomentum;

            double finalProbE = 1e-6;
            double finalProbH = 1e-6;

            double minDeltaP = 999.0;

            int bestRow = 0;


            for (int iLE = 0; iLE < nRows; iLE++) {
                muonEcalLookupTree->GetEntry(iLE);
                double deltaP = TMath::Abs(lookup_E_p - particleMomentum);
                
                if (deltaP < minDeltaP) {
                    minDeltaP = deltaP;
                    bestRow = iLE;
                } else {
                    break; 
                }
            }
            muonEcalLookupTree->GetEntry(bestRow);
            double minDeltaEp = 999.0;
            int bestCol = 0;

            for (size_t jLE = 0; jLE < Ep_E_ptr->size(); jLE++) {
                double deltaEp = TMath::Abs(Ep_E_ptr->at(jLE) - Ep_Ecal);
                
                if (deltaEp < minDeltaEp) {
                    minDeltaEp = deltaEp;
                    bestCol = jLE;
                } else {
                    break;
                }
            }

            finalProbE = prob_E_ptr->at(bestCol);

            if (finalProbE <= 0)
            {
                if (bestCol == 0) finalProbE = prob_E_ptr->at(bestCol+1);
                else if (bestCol == Ep_E_ptr->size() - 1) finalProbE = prob_E_ptr->at(bestCol-1);
                else finalProbE = (prob_E_ptr->at(bestCol-1) + prob_E_ptr->at(bestCol+1)) / 2.0;
            }

            //std::cout << "Final Ecal Probability: " << finalProbE << std::endl;
            
            minDeltaP = 999.0;

            bestRow = 0;

            for (int iLH = 0; iLH < nRows; iLH++) {
                muonHcalLookupTree->GetEntry(iLH);
                double deltaP = TMath::Abs(lookup_H_p - particleMomentum);
                
                if (deltaP < minDeltaP) {
                    minDeltaP = deltaP;
                    bestRow = iLH;
                } else {
                    break; 
                }
            }
            muonHcalLookupTree->GetEntry(bestRow);
            minDeltaEp = 999.0;
            bestCol = 0;

            for (size_t jLH = 0; jLH < Ep_H_ptr->size(); jLH++) {
                double deltaEp = TMath::Abs(Ep_H_ptr->at(jLH) - Ep_Hcal);
                
                if (deltaEp < minDeltaEp) {
                    minDeltaEp = deltaEp;
                    bestCol = jLH;
                } else {
                    break;
                }
            }

            finalProbH = prob_H_ptr->at(bestCol);

            if (finalProbH <= 0)
            {
                if (bestCol == 0) finalProbH = prob_H_ptr->at(bestCol+1);
                else if (bestCol == Ep_H_ptr->size() - 1) finalProbH = prob_H_ptr->at(bestCol-1);
                else finalProbH = (prob_H_ptr->at(bestCol-1) + prob_H_ptr->at(bestCol+1)) / 2.0;
            }

            //std::cout << "Final Hcal Probability: " << finalProbH << std::endl;

            if (ECalEnergy == 0.0) finalProbE = 1; // Ignore Ecal probability if no energy is deposited in Ecal

            if (finalProbE <= 0) finalProbE = 1e-6; // Avoid log(0)
            if (finalProbH <= 0) finalProbH = 1e-6; // Avoid log(0)

            double logL_muon = std::log(finalProbE) + std::log(finalProbH);

            if (truePDG == 13) // If it's a muon
            {
                h_muonLikelihood->Fill(logL_muon);
                h_muonLikelihoodVsp->Fill(particleMomentum, logL_muon);
                //std::cout << "Log-Likelihood: " << logL_muon << std::endl;
            }

            if (truePDG != 13) // If it's not a muon
            {
                h_nonMuonLikelihood->Fill(logL_muon);
                h_nonMuonLikelihoodVsp->Fill(particleMomentum, logL_muon);
                //std::cout << "Log-Likelihood: " << logL_muon << std::endl;
            }

            if (finalProbE < 1e-6 && finalProbH == 1e-6)
            {
                std::cout << "Event ID: " << eventID << std::endl;
                std::cout << "Particle Momentum: " << particleMomentum << std::endl;
                std::cout << "Particle Eta: " << particleEta << std::endl;
                std::cout << "True PDG: " << truePDG << std::endl;
                std::cout << "ECal Energy: " << ECalEnergy << std::endl;
                std::cout << "HCal Energy: " << HCalEnergy << std::endl;
                std::cout << "Ep ECal: " << Ep_Ecal << std::endl;
                std::cout << "Ep HCal: " << Ep_Hcal << std::endl;
                std::cout << "Final Ecal Probability: " << finalProbE << std::endl;
                std::cout << "Final Hcal Probability: " << finalProbH << std::endl;
                std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

            }

        }
    }

    TCanvas *c_muonLikelihood = new TCanvas("c_muonLikelihood", "Muon Likelihood", 800, 600);
    h_muonLikelihood->SetMarkerStyle(20);
    h_muonLikelihood->SetMarkerColor(sixColourScheme[0]);
    h_muonLikelihood->Draw("HIST P");
    h_nonMuonLikelihood->SetMarkerStyle(20);
    h_nonMuonLikelihood->SetMarkerColor(sixColourScheme[1]);
    h_nonMuonLikelihood->Draw("HIST P SAME");
    c_muonLikelihood->Update();

    TCanvas *c_muonLikelihoodVsp = new TCanvas("c_muonLikelihoodVsp", "Muon Likelihood vs p", 800, 600);
    h_muonLikelihoodVsp->Draw("COLZ");
    c_muonLikelihoodVsp->Update();

    TCanvas *c_nonMuonLikelihoodVsp = new TCanvas("c_nonMuonLikelihoodVsp", "Non-Muon Likelihood vs p", 800, 600);
    h_nonMuonLikelihoodVsp->Draw("COLZ");
    c_nonMuonLikelihoodVsp->Update();

    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    h_muonLikelihood->Write();
    h_nonMuonLikelihood->Write();
    h_muonLikelihoodVsp->Write();
    h_nonMuonLikelihoodVsp->Write();
    outFile->Close();

}