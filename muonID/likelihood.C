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

struct EpLookupTree {
    TTree* tree = nullptr;
    std::vector<double> *Ep_ptr = nullptr;
    std::vector<double> *prob_ptr = nullptr;
    double p;
};

struct SizeLookupTree {
    TTree* tree = nullptr;
    std::vector<double> *size_ptr = nullptr;
    std::vector<double> *prob_ptr = nullptr;
    double p;
};

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

    TString lookupN = "mu_20-100GeV_lookupTable.root";
    TFile *lookupFile = new TFile(lookupN);

    TString infile="../eicReconOutput/SimCampaign_26020_JPsiMuMu_10ifb_10x130ep_Pruned.root";
    //TString infile="reconOut/pi_20-100GeV_reconOut.root";
    //TString infile="../dis_background/DIS_Q2_1_100_10x130ep_Pruned.root";

    double EndcapNHcal_Factor = 1.0/6.0;

    std::string outfilename = "outputs/SimCampaign_26020_JPsiMuMu_10ifb_10x130ep_likelihoodPlots.root";
    //std::string outfilename = "outputs/pi_20-100GeV_likelihoodPlots.root";
    //std::string outfilename = "outputs/DIS_Q2_1_100_10x130ep_likelihoodPlots.root";


    std::string detectorNames[8] = {"EcalBarrelImaging", "EcalBarrelScFi", "EcalEndcapP", "EcalEndcapN", "HcalBarrel", "HcalEndcapP", "LFHcal", "HcalEndcapN"};

    EpLookupTree epTrees[8];

    for (int i = 0; i < 8; i++) {
        TString treeName = Form("%sEpLookupTree", detectorNames[i].c_str());
        epTrees[i].tree = (TTree*)lookupFile->Get(treeName);
        
        if (!epTrees[i].tree) {
            printf("Warning: %s not found!\n", treeName.Data());
            continue;
        }

        epTrees[i].tree->SetBranchAddress("trueEp", &epTrees[i].Ep_ptr);
        epTrees[i].tree->SetBranchAddress("probability", &epTrees[i].prob_ptr);
        epTrees[i].tree->SetBranchAddress("p", &epTrees[i].p);
    }

    SizeLookupTree sizeTrees[8];

    for (int i = 0; i < 8; i++) {
        TString treeName = Form("%sSizeLookupTree", detectorNames[i].c_str());
        sizeTrees[i].tree = (TTree*)lookupFile->Get(treeName);
        
        if (!sizeTrees[i].tree) {
            printf("Warning: %s not found!\n", treeName.Data());
            continue;
        }

        sizeTrees[i].tree->SetBranchAddress("size", &sizeTrees[i].size_ptr);
        sizeTrees[i].tree->SetBranchAddress("probability", &sizeTrees[i].prob_ptr);
        sizeTrees[i].tree->SetBranchAddress("p", &sizeTrees[i].p);
    }


    double pMin = epTrees[0].tree->GetMinimum("p");

    double pMax = epTrees[0].tree->GetMaximum("p");

    int nRows = epTrees[0].tree->GetEntries();

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
    
    TTreeReaderArray<int> simuAssocEcalBarrelImg(particleSource_Reader, "_EcalBarrelImagingClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocEcalBarrelImg(particleSource_Reader, "_EcalBarrelImagingClusterAssociations_rec.index");
    TTreeReaderArray<float> EcalBarrelImgEng(particleSource_Reader, "EcalBarrelImagingClusters.energy");
    TTreeReaderArray<int> simuAssocEcalBarrelScFi(particleSource_Reader, "_EcalBarrelScFiClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocEcalBarrelScFi(particleSource_Reader, "_EcalBarrelScFiClusterAssociations_rec.index");
    TTreeReaderArray<float> EcalBarrelScFiEng(particleSource_Reader, "EcalBarrelScFiClusters.energy");
    

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

    TH1D *h_muonLikelihood = new TH1D("h_muonLikelihood", "Muon Likelihood;Likelihood;Entries", 100.0, -30.0, 0.0);
    TH2D *h_muonLikelihoodVsp = new TH2D("h_muonLikelihoodVsp", "Muon Likelihood vs p;Momentum (GeV);Likelihood", 500, 0.0, 100.0, 100.0, -30.0, 0.0);
    TH2D *h_muonLikelihoodVsEta = new TH2D("h_muonLikelihoodVsEta", "Muon Likelihood vs Eta;Eta;Likelihood", 100, -4.0, 4.0, 100.0, -30.0, 0.0);

    TH1D *h_nonMuonLikelihood = new TH1D("h_nonMuonLikelihood", "Non-Muon Likelihood;Likelihood;Entries", 100.0, -30.0, 0.0);
    TH2D *h_nonMuonLikelihoodVsp = new TH2D("h_nonMuonLikelihoodVsp", "Non-Muon Likelihood vs p;Momentum (GeV);Likelihood", 500, 0.0, 100.0, 100.0, -30.0, 0.0);
    TH2D *h_nonMuonLikelihoodVsEta = new TH2D("h_nonMuonLikelihoodVsEta", "Non-Muon Likelihood vs Eta;Eta;Likelihood", 100, -4.0, 4.0, 100.0, -30.0, 0.0);
    
    int eventID = 0;

    double energy[8], size[8];

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

        if (eventID > 50000) break; // Only run over first 50k events for testing purposes

        for (int i = 0; i < trackPDG.GetSize(); i++)
        {
            if (i >= simuAssoc.GetSize() || simuAssoc[i] >= partPdg.GetSize()) continue;

            for (int j = 0; j < 8; j++)
            {
                energy[j] = 0.;
                size[j] = 0.;
            }

            TVector3 recoMom(recoMomX[i],recoMomY[i],recoMomZ[i]);

            particleMomentum = recoMom.Mag();
            particleEta = recoMom.Eta();
            truePDG = TMath::Abs(partPdg[simuAssoc[i]]);

            if (TMath::Abs(particleEta) > 4) continue;
            if (TMath::Abs(particleEta) >= 1.0 && TMath::Abs(particleEta) <= 1.3) continue; // Remove overlap region between barrel and endcap calorimeters
            if (particleMomentum < 0.05 || particleMomentum > 100.0) continue;


            if (EcalBarrelImgEng.GetSize() == simuAssocEcalBarrelImg.GetSize())
            {
                for (int jI = 0; jI < simuAssocEcalBarrelImg.GetSize(); jI++) // Look for associations in the Ecal Barrel Imaging
                {
                    if (simuAssocEcalBarrelImg[jI] == simuAssoc[i])
                    {
                        energy[0] += EcalBarrelImgEng[jI];
                        size[0] += 1.;
                    }
                }
            }
            if (EcalBarrelScFiEng.GetSize() == simuAssocEcalBarrelScFi.GetSize())
            {
                for (int jS = 0; jS < simuAssocEcalBarrelScFi.GetSize(); jS++) // Look for associations in the Ecal Barrel ScFi
                {
                    if (simuAssocEcalBarrelScFi[jS] == simuAssoc[i])
                    {
                        energy[1] += EcalBarrelScFiEng[jS];
                        size[1] += 1.;
                    }
                }
            }
            
            if (EcalEndcapPEng.GetSize() == simuAssocEcalEndcapP.GetSize()) 
            {
                for (int jP = 0; jP < simuAssocEcalEndcapP.GetSize(); jP++) // Look for associations in the Ecal Endcap P
                {
                    if (simuAssocEcalEndcapP[jP] == simuAssoc[i])
                    {
                        energy[2] += EcalEndcapPEng[jP];
                        size[2] += 1.;
                    }
                }
            }
            if (EcalEndcapNEng.GetSize() == simuAssocEcalEndcapN.GetSize())
            {
                for (int jN = 0; jN < simuAssocEcalEndcapN.GetSize(); jN++) // Look for associations in the Ecal Endcap N
                {
                    if (simuAssocEcalEndcapN[jN] == simuAssoc[i])
                    {
                        energy[3] += EcalEndcapNEng[jN];
                        size[3] += 1.;
                    }
                }
            }   

            if (HcalBarrelEng.GetSize() == simuAssocHcalBarrel.GetSize())
            {
                for (int jB = 0; jB < simuAssocHcalBarrel.GetSize(); jB++) // Look for associations in the HCal Barrel
                {
                    if (simuAssocHcalBarrel[jB] == simuAssoc[i])
                    {
                        energy[4] += HcalBarrelEng[jB];
                        size[4] += 1.;
                    }
                }
            }
            
            if (HcalEndcapPEng.GetSize() == simuAssocHcalEndcapP.GetSize())
            {
                for (int jP = 0; jP < simuAssocHcalEndcapP.GetSize(); jP++) // Look for associations in the HCal Endcap P
                {
                    if (simuAssocHcalEndcapP[jP] == simuAssoc[i])
                    {
                        energy[5] += HcalEndcapPEng[jP];
                        size[5] += 1.;
                    }
                }
            }

            if (LFHcalEng.GetSize() == simuAssocLFHcal.GetSize())
            {
                for (int jL = 0; jL < simuAssocLFHcal.GetSize(); jL++) // Look for associations in the LFHcal
                {
                    if (simuAssocLFHcal[jL] == simuAssoc[i])
                    {
                        energy[6] += LFHcalEng[jL];
                        size[6] += 1.;
                    }
                }
            }
            
            if (HcalEndcapNEng.GetSize() == simuAssocHcalEndcapN.GetSize())
            {
                for (int jN = 0; jN < simuAssocHcalEndcapN.GetSize(); jN++) // Look for associations in the HCal Endcap N
                {
                    if (simuAssocHcalEndcapN[jN] == simuAssoc[i])
                    {
                        energy[7] += HcalEndcapNEng[jN] * EndcapNHcal_Factor;
                        size[7] += 1.;
                    }
                }
            }

            if (energy[0] == 0.0 && energy[1] == 0.0 && energy[2] == 0.0 && energy[3] == 0.0 && energy[4] == 0.0 && energy[5] == 0.0 && energy[6] == 0.0 && energy[7] == 0.0) continue; // Skip if no calorimeter energy is associated


            double finalProbEp[8] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
            double finalProbSize[8] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};

            double minDeltaP = 999.0;

            int bestRow = 0;

            for (int detectorsI = 0; detectorsI < 8; detectorsI++)
            {
                if (energy[detectorsI] == 0.0) continue; // Skip if no energy is deposited in this detector
                for (int iEP = 0; iEP < nRows; iEP++) {
                    epTrees[detectorsI].tree->GetEntry(iEP);
                    double deltaP = TMath::Abs(epTrees[detectorsI].p - particleMomentum);
                    
                    if (deltaP < minDeltaP) {
                        minDeltaP = deltaP;
                        bestRow = iEP;
                    } else {
                        break; 
                    }
                }
                epTrees[detectorsI].tree->GetEntry(bestRow);
                double minDeltaEp = 999.0;
                int bestCol = 0;

                for (size_t jEP = 0; jEP < epTrees[detectorsI].Ep_ptr->size(); jEP++) {
                    double deltaEp = TMath::Abs(epTrees[detectorsI].Ep_ptr->at(jEP) - energy[detectorsI]/particleMomentum);
                    
                    if (deltaEp < minDeltaEp) {
                        minDeltaEp = deltaEp;
                        bestCol = jEP;
                    } else {
                        break;
                    }
                }

                finalProbEp[detectorsI] = epTrees[detectorsI].prob_ptr->at(bestCol);

                if (finalProbEp[detectorsI] <= 0)
                {
                    if (bestCol == 0) finalProbEp[detectorsI] = epTrees[detectorsI].prob_ptr->at(bestCol+1);
                    else if (bestCol == epTrees[detectorsI].Ep_ptr->size() - 1) finalProbEp[detectorsI] = epTrees[detectorsI].prob_ptr->at(bestCol-1);
                    else finalProbEp[detectorsI] = (epTrees[detectorsI].prob_ptr->at(bestCol-1) + epTrees[detectorsI].prob_ptr->at(bestCol+1)) / 2.0;
                }           

            
                minDeltaP = 999.0;

                bestRow = 0;

                for (int iSize = 0; iSize < nRows; iSize++) {
                    sizeTrees[detectorsI].tree->GetEntry(iSize);
                    double deltaP = TMath::Abs(sizeTrees[detectorsI].p - particleMomentum);
                    
                    if (deltaP < minDeltaP) {
                        minDeltaP = deltaP;
                        bestRow = iSize;
                    } else {
                        break; 
                    }
                }
                sizeTrees[detectorsI].tree->GetEntry(bestRow);
                minDeltaEp = 999.0;
                bestCol = 0;

                for (size_t jSize = 0; jSize < sizeTrees[detectorsI].size_ptr->size(); jSize++) {
                    double deltaEp = TMath::Abs(sizeTrees[detectorsI].size_ptr->at(jSize) - size[detectorsI]);
                    
                    if (deltaEp < minDeltaEp) {
                        minDeltaEp = deltaEp;
                        bestCol = jSize;
                    } else {
                        break;
                    }
                }

                finalProbSize[detectorsI] = sizeTrees[detectorsI].prob_ptr->at(bestCol);

                if (finalProbSize[detectorsI] <= 0)
                {
                    if (bestCol == 0) finalProbSize[detectorsI] = sizeTrees[detectorsI].prob_ptr->at(bestCol+1);
                    else if (bestCol == sizeTrees[detectorsI].size_ptr->size() - 1) finalProbSize[detectorsI] = sizeTrees[detectorsI].prob_ptr->at(bestCol-1);
                    else finalProbSize[detectorsI] = (sizeTrees[detectorsI].prob_ptr->at(bestCol-1) + sizeTrees[detectorsI].prob_ptr->at(bestCol+1)) / 2.0;
                }

            }

            //std::cout << "Final Hcal Probability: " << finalProbH << std::endl;

            if (energy[0] == 0.0 && energy[1] == 0.0 && energy[2] == 0.0 && energy[3] == 0.0) // Ignore Ecal probability if no energy is deposited in Ecal
            {
                finalProbEp[0] = 1;
                finalProbSize[0] = 1;
            }

            for (int k = 1; k < 8; k++) // Reset any probabilities that are 0 or negative to a small value to avoid issues with log-likelihood calculation
            {
                if (finalProbEp[k] <= 0) finalProbEp[k] = 1e-6;
                if (finalProbSize[k] <= 0) finalProbSize[k] = 1e-6;
            }

            double logL_muon = 0.0;
            int measurementsCount = 0;

            for (int k = 0; k < 8; k++) {
                if (energy[k] > 0) {
                    logL_muon += TMath::Log(finalProbEp[k]);
                    logL_muon += TMath::Log(finalProbSize[k]);
                    measurementsCount += 2;
                }
            }
            if (energy[4] == 0.0 && energy[5] == 0.0 && energy[6] == 0.0 && energy[7] == 0.0) // Penalise signals with no energy in the Hcal by adding a small probability
            {
                logL_muon += TMath::Log(1e-6);
                logL_muon += TMath::Log(1e-6);
                measurementsCount += 2;
            }
            if (measurementsCount > 0) logL_muon /= (measurementsCount/2.0);


            if (truePDG == 13) // If it's a muon
            {
                h_muonLikelihood->Fill(logL_muon);
                h_muonLikelihoodVsp->Fill(particleMomentum, logL_muon);
                h_muonLikelihoodVsEta->Fill(particleEta, logL_muon);
                //std::cout << "Log-Likelihood: " << logL_muon << std::endl;
            }

            if (truePDG != 13) // If it's not a muon
            {
                h_nonMuonLikelihood->Fill(logL_muon);
                h_nonMuonLikelihoodVsp->Fill(particleMomentum, logL_muon);
                h_nonMuonLikelihoodVsEta->Fill(particleEta, logL_muon);
                //std::cout << "Log-Likelihood: " << logL_muon << std::endl;
            }


        }
    }

    TCanvas *c_muonLikelihood = new TCanvas("c_muonLikelihood", "Muon Likelihood", 800, 600);
    h_muonLikelihood->SetMarkerStyle(20);
    h_muonLikelihood->SetMarkerColor(sixColourScheme[0]);
    h_nonMuonLikelihood->SetMarkerStyle(20);
    h_nonMuonLikelihood->SetMarkerColor(sixColourScheme[1]);
    if (h_muonLikelihood->Integral() >= h_nonMuonLikelihood->Integral())
    {
        h_muonLikelihood->Draw("HIST P");
        h_nonMuonLikelihood->Draw("HIST P SAME");
    }
    else
    {
        h_nonMuonLikelihood->Draw("HIST P");
        h_muonLikelihood->Draw("HIST P SAME");
    }
    c_muonLikelihood->Update();

    TCanvas *c_muonLikelihoodVsp = new TCanvas("c_muonLikelihoodVsp", "Muon Likelihood vs p", 800, 600);
    h_muonLikelihoodVsp->Draw("COLZ");
    c_muonLikelihoodVsp->Update();

    TCanvas *c_nonMuonLikelihoodVsp = new TCanvas("c_nonMuonLikelihoodVsp", "Non-Muon Likelihood vs p", 800, 600);
    h_nonMuonLikelihoodVsp->Draw("COLZ");
    c_nonMuonLikelihoodVsp->Update();

    TCanvas *c_muonLikelihoodVsEta = new TCanvas("c_muonLikelihoodVsEta", "Muon Likelihood vs Eta", 800, 600);
    h_muonLikelihoodVsEta->Draw("COLZ");
    c_muonLikelihoodVsEta->Update();

    TCanvas *c_nonMuonLikelihoodVsEta = new TCanvas("c_nonMuonLikelihoodVsEta", "Non-Muon Likelihood vs Eta", 800, 600);
    h_nonMuonLikelihoodVsEta->Draw("COLZ");
    c_nonMuonLikelihoodVsEta->Update();

    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");
    h_muonLikelihood->Write();
    h_nonMuonLikelihood->Write();
    h_muonLikelihoodVsp->Write();
    h_nonMuonLikelihoodVsp->Write();
    h_muonLikelihoodVsEta->Write();
    h_nonMuonLikelihoodVsEta->Write();
    outFile->Close();

}