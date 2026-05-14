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

void lookupTableMaker()
{
    gROOT->SetBatch(kTRUE);

    TString infile="reconOut/mu_20-100GeV_reconOut_Pruned.root";

    double EndcapNHcal_Factor = 1.0/6.0;

    std::string outfilename = "mu_20-100GeV_lookupTable.root";

    // Set up input file chain
    TChain *mychain = new TChain("events");
    mychain->Add(infile);

    // Initialize reader
    TTreeReader tree_reader(mychain);

    // Get Particle Information
    TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<double> partMomX(tree_reader, "MCParticles.momentum.x");
    TTreeReaderArray<double> partMomY(tree_reader, "MCParticles.momentum.y");
    TTreeReaderArray<double> partMomZ(tree_reader, "MCParticles.momentum.z");
    TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");
    TTreeReaderArray<double> partMass(tree_reader, "MCParticles.mass");
    TTreeReaderArray<float> partCharge(tree_reader, "MCParticles.charge");
    TTreeReaderArray<unsigned int> partParb(tree_reader, "MCParticles.parents_begin");
    TTreeReaderArray<unsigned int> partPare(tree_reader, "MCParticles.parents_end");
    TTreeReaderArray<int> partParI(tree_reader, "_MCParticles_parents.index");

    // Get Reconstructed Track Information
    TTreeReaderArray<float> recoMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float> recoMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float> recoMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");
    TTreeReaderArray<int> trackPDG(tree_reader, "ReconstructedChargedParticles.PDG");
    TTreeReaderArray<float> trackMass(tree_reader, "ReconstructedChargedParticles.mass");
    TTreeReaderArray<float> trackCharge(tree_reader, "ReconstructedChargedParticles.charge");
    TTreeReaderArray<float> trackEng(tree_reader, "ReconstructedChargedParticles.energy");

    // Get Truth Seeded Reconstructed Track Information
    TTreeReaderArray<float> trecoMomX(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x");
    TTreeReaderArray<float> trecoMomY(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y");
    TTreeReaderArray<float> trecoMomZ(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z");
    TTreeReaderArray<int> ttrackPDG(tree_reader, "ReconstructedTruthSeededChargedParticles.PDG");
    TTreeReaderArray<float> ttrackMass(tree_reader, "ReconstructedTruthSeededChargedParticles.mass");
    TTreeReaderArray<float> ttrackCharge(tree_reader, "ReconstructedTruthSeededChargedParticles.charge");
    TTreeReaderArray<float> ttrackEng(tree_reader, "ReconstructedTruthSeededChargedParticles.energy");

    // Get Associations Between MCParticles and ReconstructedChargedParticles
    TTreeReaderArray<int> recoAssoc(tree_reader, "_ReconstructedChargedParticleAssociations_rec.index");
    TTreeReaderArray<int> simuAssoc(tree_reader, "_ReconstructedChargedParticleAssociations_sim.index");

    // Get B0 Information
    TTreeReaderArray<int> recoAssocB0(tree_reader, "_B0ECalClusterAssociations_rec.index");
    TTreeReaderArray<int> simuAssocB0(tree_reader, "_B0ECalClusterAssociations_sim.index");
    TTreeReaderArray<float> B0Eng(tree_reader, "B0ECalClusters.energy");
    TTreeReaderArray<float> B0z(tree_reader, "B0ECalClusters.position.z");

    // Get Forward Detector Information
    TTreeReaderArray<float> RPEng(tree_reader, "ForwardRomanPotRecParticles.energy");
    TTreeReaderArray<int> RPpdg(tree_reader, "ForwardRomanPotRecParticles.PDG");
    TTreeReaderArray<float> RPMomX(tree_reader, "ForwardRomanPotRecParticles.momentum.x");
    TTreeReaderArray<float> RPMomY(tree_reader, "ForwardRomanPotRecParticles.momentum.y");
    TTreeReaderArray<float> RPMomZ(tree_reader, "ForwardRomanPotRecParticles.momentum.z");

    TTreeReaderArray<float> OffMEng(tree_reader, "ForwardOffMRecParticles.energy");

    // Ecal Information
    TTreeReaderArray<int> simuAssocEcalBarrel(tree_reader, "_EcalBarrelClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocEcalBarrel(tree_reader, "_EcalBarrelClusterAssociations_rec.index");
    TTreeReaderArray<float> EcalBarrelEng(tree_reader, "EcalBarrelClusters.energy");
    
    TTreeReaderArray<int> simuAssocEcalBarrelImg(tree_reader, "_EcalBarrelImagingClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocEcalBarrelImg(tree_reader, "_EcalBarrelImagingClusterAssociations_rec.index");
    TTreeReaderArray<float> EcalBarrelImgEng(tree_reader, "EcalBarrelImagingClusters.energy");
    TTreeReaderArray<int> simuAssocEcalBarrelScFi(tree_reader, "_EcalBarrelScFiClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocEcalBarrelScFi(tree_reader, "_EcalBarrelScFiClusterAssociations_rec.index");
    TTreeReaderArray<float> EcalBarrelScFing(tree_reader, "EcalBarrelScFiClusters.energy");
    

    TTreeReaderArray<int> simuAssocEcalEndcapP(tree_reader, "_EcalEndcapPClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocEcalEndcapP(tree_reader, "_EcalEndcapPClusterAssociations_rec.index");    
    TTreeReaderArray<float> EcalEndcapPEng(tree_reader, "EcalEndcapPClusters.energy");

    TTreeReaderArray<int> simuAssocEcalEndcapN(tree_reader, "_EcalEndcapNClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocEcalEndcapN(tree_reader, "_EcalEndcapNClusterAssociations_rec.index");
    TTreeReaderArray<float> EcalEndcapNEng(tree_reader, "EcalEndcapNClusters.energy");

    // Hcal Information
    TTreeReaderArray<int> simuAssocHcalBarrel(tree_reader, "_HcalBarrelClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocHcalBarrel(tree_reader, "_HcalBarrelClusterAssociations_rec.index");
    TTreeReaderArray<float> HcalBarrelEng(tree_reader, "HcalBarrelClusters.energy");

    TTreeReaderArray<int> simuAssocHcalEndcapP(tree_reader, "_HcalEndcapPInsertClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocHcalEndcapP(tree_reader, "_HcalEndcapPInsertClusterAssociations_rec.index");    
    TTreeReaderArray<float> HcalEndcapPEng(tree_reader, "HcalEndcapPInsertClusters.energy");

    TTreeReaderArray<int> simuAssocLFHcal(tree_reader, "_LFHCALClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocLFHcal(tree_reader, "_LFHCALClusterAssociations_rec.index");    
    TTreeReaderArray<float> LFHcalEng(tree_reader, "LFHCALClusters.energy");

    TTreeReaderArray<int> simuAssocHcalEndcapN(tree_reader, "_HcalEndcapNClusterAssociations_sim.index");
    TTreeReaderArray<int> recoAssocHcalEndcapN(tree_reader, "_HcalEndcapNClusterAssociations_rec.index");
    TTreeReaderArray<float> HcalEndcapNEng(tree_reader, "HcalEndcapNClusters.energy");

    // Histograms for E/p and cluster size

    std::string detectors[8] = {"EcalBarrelImaging", "EcalBarrelScFi", "EcalEndcapP", "EcalEndcapN", "HcalBarrel", "HcalEndcapP", "LFHcal", "HcalEndcapN"};

    TH2D *h_EpTrue[8];

    h_EpTrue[0] = new TH2D("h_EpTrueEcalBarrelImaging","Ecal Energy vs Track Momentum for Muons; p;  E/p",50,0.,100.,100,0.,2.0);
    h_EpTrue[1] = new TH2D("h_EpTrueEcalBarrelScFi","Ecal Energy vs Track Momentum for Muons; p;  E/p",50,0.,100.,100,0.,2.0);
    h_EpTrue[2] = new TH2D("h_EpTrueEcalEndcapP","Ecal Energy vs Track Momentum for Muons; p;  E/p",50,0.,100.,100,0.,2.0);
    h_EpTrue[3] = new TH2D("h_EpTrueEcalEndcapN","Ecal Energy vs Track Momentum for Muons; p;  E/p",50,0.,100.,100,0.,2.0);

    h_EpTrue[4] = new TH2D("h_EpTrueHcalBarrel","Hcal Energy vs Track Momentum for Muons; p;  E/p",50,0.,100.,100,0.,2.0);
    h_EpTrue[5] = new TH2D("h_EpTrueHcalEndcapP","Hcal Energy vs Track Momentum for Muons; p;  E/p",50,0.,100.,100,0.,2.0);
    h_EpTrue[6] = new TH2D("h_EpTrueLFHcal","LFHcal Energy vs Track Momentum for Muons; p;  E/p",50,0.,100.,100,0.,2.0);
    h_EpTrue[7] = new TH2D("h_EpTrueHcalEndcapN","Hcal Energy vs Track Momentum for Muons; p;  E/p",50,0.,100.,100,0.,2.0);

    TH2D *h_SizeVsp[8];

    h_SizeVsp[0] = new TH2D("h_EcalBarrelImagingHitsSizeVsE","Ecal Barrel Imaging Hits Size vs Energy; Energy (GeV); Size",50,0.,100.,100,0.,100.);
    h_SizeVsp[1] = new TH2D("h_EcalBarrelScFiHitsSizeVsE","Ecal Barrel ScFi Hits Size vs Energy; Energy (GeV); Size",50,0.,100.,100,0.,100.);
    h_SizeVsp[2] = new TH2D("h_EcalEndcapPHitsSizeVsE","Ecal Endcap P Hits Size vs Energy; Energy (GeV); Size",50,0.,100.,100,0.,100.);
    h_SizeVsp[3] = new TH2D("h_EcalEndcapNHitsSizeVsE","Ecal Endcap N Hits Size vs Energy; Energy (GeV); Size",50,0.,100.,100,0.,100.);

    h_SizeVsp[4] = new TH2D("h_HcalBarrelHitsSizeVsE","Hcal Barrel Hits Size vs Energy; Energy (GeV); Size",50,0.,100.,100,0.,100.);
    h_SizeVsp[5] = new TH2D("h_HcalEndcapPHitsSizeVsE","Hcal Endcap P Hits Size vs Energy; Energy (GeV); Size",50,0.,100.,100,0.,100.);
    h_SizeVsp[6] = new TH2D("h_LFHcalHitsSizeVsE","LFHcal Hits Size vs Energy; Energy (GeV); Size",50,0.,100.,100,0.,100.);
    h_SizeVsp[7] = new TH2D("h_HcalEndcapNHitsSizeVsE","Hcal Endcap N Hits Size vs Energy; Energy (GeV); Size",50,0.,100.,100,0.,100.);

    double energy[8], size[8];

    int eventID = 0;

    while(tree_reader.Next()) // Loop over events
    {
        eventID++;

        if (eventID % 10 == 0)
        {
            fprintf (stderr, "%4.2f Percent\r ", eventID*100.0/mychain->GetEntries());
            fflush (stderr);
        }


        for (int i = 0; i < trackPDG.GetSize(); i++)
        {
            for (int j = 0; j < 8; j++)
            {
                energy[j] = 0.;
                size[j] = 0.;
            }

            TVector3 recoMom(recoMomX[i],recoMomY[i],recoMomZ[i]);
            ROOT::Math::PxPyPzEVector reco4Mom(recoMomX[i],recoMomY[i],recoMomZ[i], trackEng[i]);

            if (i >= simuAssoc.GetSize() || simuAssoc[i] >= partPdg.GetSize()) continue;
            
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
            if (EcalBarrelScFing.GetSize() == simuAssocEcalBarrelScFi.GetSize())
            {
                for (int jS = 0; jS < simuAssocEcalBarrelScFi.GetSize(); jS++) // Look for associations in the Ecal Barrel ScFi
                {
                    if (simuAssocEcalBarrelScFi[jS] == simuAssoc[i])
                    {
                        energy[1] += EcalBarrelScFing[jS];
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

            TVector3 trueMom(partMomX[simuAssoc[i]],partMomY[simuAssoc[i]],partMomZ[simuAssoc[i]]);

            if (energy[0] > 0.) h_EpTrue[0]->Fill(trueMom.Mag(), energy[0]/trueMom.Mag());
            if (energy[1] > 0.) h_EpTrue[1]->Fill(trueMom.Mag(), energy[1]/trueMom.Mag());
            if (energy[2] > 0.) h_EpTrue[2]->Fill(trueMom.Mag(), energy[2]/trueMom.Mag());
            if (energy[3] > 0.) h_EpTrue[3]->Fill(trueMom.Mag(), energy[3]/trueMom.Mag());
            if (energy[4] > 0.) h_EpTrue[4]->Fill(trueMom.Mag(), energy[4]/trueMom.Mag());
            if (energy[5] > 0.) h_EpTrue[5]->Fill(trueMom.Mag(), energy[5]/trueMom.Mag());
            if (energy[6] > 0.) h_EpTrue[6]->Fill(trueMom.Mag(), energy[6]/trueMom.Mag());
            if (energy[7] > 0.) h_EpTrue[7]->Fill(trueMom.Mag(), energy[7]/trueMom.Mag());

            if (size[0] > 0.) h_SizeVsp[0]->Fill(trueMom.Mag(), size[0]);
            if (size[1] > 0.) h_SizeVsp[1]->Fill(trueMom.Mag(), size[1]);
            if (size[2] > 0.) h_SizeVsp[2]->Fill(trueMom.Mag(), size[2]);
            if (size[3] > 0.) h_SizeVsp[3]->Fill(trueMom.Mag(), size[3]);
            if (size[4] > 0.) h_SizeVsp[4]->Fill(trueMom.Mag(), size[4]);
            if (size[5] > 0.) h_SizeVsp[5]->Fill(trueMom.Mag(), size[5]);
            if (size[6] > 0.) h_SizeVsp[6]->Fill(trueMom.Mag(), size[6]);
            if (size[7] > 0.) h_SizeVsp[7]->Fill(trueMom.Mag(), size[7]);

        }
    }

    // New Histograms for the lookup tables

    TH1D *h_EpTrueProbabilities[8][50];

    TH1D *h_SizeProbabilities[8][50];

    std::vector<std::vector<std::vector<double>>> EpTrueLookupTable(8, std::vector<std::vector<double>>(h_EpTrue[0]->GetNbinsX() + 1, std::vector<double>(h_EpTrue[0]->GetNbinsY() + 1, 0.0)));

    std::vector<std::vector<std::vector<double>>> SizeVsELookupTable(8, std::vector<std::vector<double>>(h_SizeVsp[0]->GetNbinsX() + 1, std::vector<double>(h_SizeVsp[0]->GetNbinsY() + 1, 0.0)));

    std::vector<double> pValues(h_EpTrue[0]->GetNbinsX() + 1, 0.0);
    std::vector<double> EpValues(h_EpTrue[0]->GetNbinsY() + 1, 0.0);
    std::vector<double> sizeValues(h_SizeVsp[0]->GetNbinsY() + 1, 0.0);

    for (int detectorI = 0; detectorI < 8; detectorI++)
    {
        for (int i = 1; i < h_EpTrue[detectorI]->GetNbinsX(); i++) {
            pValues.at(i) = h_EpTrue[detectorI]->GetXaxis()->GetBinCenter(i);
            h_EpTrueProbabilities[detectorI][i] = h_EpTrue[detectorI]->ProjectionY(Form("probabilities_det%d_bin%d", detectorI, i), i, i);
            for (int j = 0; j < h_EpTrue[detectorI]->GetNbinsY(); j++) {
                if (h_EpTrueProbabilities[detectorI][i]->Integral() == 0) continue;
                double truep = h_EpTrue[detectorI]->GetXaxis()->GetBinCenter(i);
                double trueEp = h_EpTrue[detectorI]->GetYaxis()->GetBinCenter(j);
                double count = h_EpTrue[detectorI]->GetBinContent(i, j);
                EpValues.at(j) = trueEp;
                if (count > 0) {
                    double probability = count / h_EpTrueProbabilities[detectorI][i]->Integral();
                    int recBinX = h_EpTrue[detectorI]->GetXaxis()->FindBin(truep);
                    int recBinY = h_EpTrue[detectorI]->GetYaxis()->FindBin(trueEp);
                    if (recBinX > 0 && recBinY > 0) {
                        EpTrueLookupTable[detectorI][recBinX][recBinY] = probability;
                    }
                }
            }
            h_EpTrueProbabilities[detectorI][i]->Delete();
        }

        for (int i = 1; i < h_SizeVsp[detectorI]->GetNbinsX(); i++) {
            h_SizeProbabilities[detectorI][i] = h_SizeVsp[detectorI]->ProjectionY(Form("probabilities_det%d_bin%d", detectorI, i), i, i);
            for (int j = 0; j < h_SizeVsp[detectorI]->GetNbinsY(); j++) {
                if (h_SizeProbabilities[detectorI][i]->Integral() == 0) continue;
                double truep = h_SizeVsp[detectorI]->GetXaxis()->GetBinCenter(i);
                double trueEp = h_SizeVsp[detectorI]->GetYaxis()->GetBinCenter(j);
                double count = h_SizeVsp[detectorI]->GetBinContent(i, j);
                sizeValues.at(j) = trueEp;
                if (count > 0) {
                    double probability = count / h_SizeProbabilities[detectorI][i]->Integral();
                    int recBinX = h_SizeVsp[detectorI]->GetXaxis()->FindBin(truep);
                    int recBinY = h_SizeVsp[detectorI]->GetYaxis()->FindBin(trueEp);
                    if (recBinX > 0 && recBinY > 0) {
                        SizeVsELookupTable[detectorI][recBinX][recBinY] = probability;
                    }
                }
            }
            h_SizeProbabilities[detectorI][i]->Delete();
        }
    }


    TFile *ofile = TFile::Open(outfilename.c_str(),"RECREATE");

    for (int i = 0; i < 8; i++)
    {
        h_EpTrue[i]->Write();
        h_SizeVsp[i]->Write();
    }

    double p_val;
    std::vector<double> ep_array;
    std::vector<double> size_array;
    std::vector<double> prob_array;

    TTree *EcalBarrelImagingEpLookupTree = new TTree("EcalBarrelImagingEpLookupTree", "Ecal Barrel Imaging Ep Lookup Tree");
    EcalBarrelImagingEpLookupTree->Branch("p", &p_val, "p_val/D");
    EcalBarrelImagingEpLookupTree->Branch("trueEp", &ep_array);
    EcalBarrelImagingEpLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        ep_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < EpValues.size(); j++) {
            ep_array.push_back(EpValues.at(j));
            prob_array.push_back(EpTrueLookupTable[0][i][j]);
        }
        EcalBarrelImagingEpLookupTree->Fill();
    }

    EcalBarrelImagingEpLookupTree->Write();

    TTree *EcalBarrelScFiEpLookupTree = new TTree("EcalBarrelScFiEpLookupTree", "Ecal Barrel ScFi Ep Lookup Tree");
    EcalBarrelScFiEpLookupTree->Branch("p", &p_val, "p_val/D");
    EcalBarrelScFiEpLookupTree->Branch("trueEp", &ep_array);
    EcalBarrelScFiEpLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        ep_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < EpValues.size(); j++) {
            ep_array.push_back(EpValues.at(j));
            prob_array.push_back(EpTrueLookupTable[1][i][j]);
        }
        EcalBarrelScFiEpLookupTree->Fill();
    }

    EcalBarrelScFiEpLookupTree->Write();

    TTree *EcalEndcapPEpLookupTree = new TTree("EcalEndcapPEpLookupTree", "Ecal Endcap P Ep Lookup Tree");
    EcalEndcapPEpLookupTree->Branch("p", &p_val, "p_val/D");
    EcalEndcapPEpLookupTree->Branch("trueEp", &ep_array);
    EcalEndcapPEpLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        ep_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < EpValues.size(); j++) {
            ep_array.push_back(EpValues.at(j));
            prob_array.push_back(EpTrueLookupTable[2][i][j]);
        }
        EcalEndcapPEpLookupTree->Fill();
    }

    EcalEndcapPEpLookupTree->Write();

    TTree *EcalEndcapNEpLookupTree = new TTree("EcalEndcapNEpLookupTree", "Ecal Endcap N Ep Lookup Tree");
    EcalEndcapNEpLookupTree->Branch("p", &p_val, "p_val/D");
    EcalEndcapNEpLookupTree->Branch("trueEp", &ep_array);
    EcalEndcapNEpLookupTree->Branch("probability", &prob_array);    

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        ep_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < EpValues.size(); j++) {
            ep_array.push_back(EpValues.at(j));
            prob_array.push_back(EpTrueLookupTable[3][i][j]);
        }
        EcalEndcapNEpLookupTree->Fill();
    }

    EcalEndcapNEpLookupTree->Write();

    TTree *HcalBarrelEpLookupTree = new TTree("HcalBarrelEpLookupTree", "Hcal Barrel Ep Lookup Tree");
    HcalBarrelEpLookupTree->Branch("p", &p_val, "p_val/D");
    HcalBarrelEpLookupTree->Branch("trueEp", &ep_array);
    HcalBarrelEpLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        ep_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < EpValues.size(); j++) {
            ep_array.push_back(EpValues.at(j));
            prob_array.push_back(EpTrueLookupTable[4][i][j]);
        }
        HcalBarrelEpLookupTree->Fill();
    }

    HcalBarrelEpLookupTree->Write();

    TTree *HcalEndcapPEpLookupTree = new TTree("HcalEndcapPEpLookupTree", "Hcal Endcap P Ep Lookup Tree");
    HcalEndcapPEpLookupTree->Branch("p", &p_val, "p_val/D");
    HcalEndcapPEpLookupTree->Branch("trueEp", &ep_array);
    HcalEndcapPEpLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        ep_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < EpValues.size(); j++) {
            ep_array.push_back(EpValues.at(j));
            prob_array.push_back(EpTrueLookupTable[5][i][j]);
        }
        HcalEndcapPEpLookupTree->Fill();
    }

    HcalEndcapPEpLookupTree->Write();

    TTree *LFHcalEpLookupTree = new TTree("LFHcalEpLookupTree", "LFHcal Ep Lookup Tree");
    LFHcalEpLookupTree->Branch("p", &p_val, "p_val/D");
    LFHcalEpLookupTree->Branch("trueEp", &ep_array);
    LFHcalEpLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        ep_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < EpValues.size(); j++) {
            ep_array.push_back(EpValues.at(j));
            prob_array.push_back(EpTrueLookupTable[6][i][j]);
        }
        LFHcalEpLookupTree->Fill();
    }   

    LFHcalEpLookupTree->Write();

    TTree *HcalEndcapNEpLookupTree = new TTree("HcalEndcapNEpLookupTree", "Hcal Endcap N Ep Lookup Tree");
    HcalEndcapNEpLookupTree->Branch("p", &p_val, "p_val/D");
    HcalEndcapNEpLookupTree->Branch("trueEp", &ep_array);
    HcalEndcapNEpLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        ep_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < EpValues.size(); j++) {
            ep_array.push_back(EpValues.at(j));
            prob_array.push_back(EpTrueLookupTable[7][i][j]);
        }
        HcalEndcapNEpLookupTree->Fill();
    }

    HcalEndcapNEpLookupTree->Write();

    TTree *EcalBarrelImagingSizeLookupTree = new TTree("EcalBarrelImagingSizeLookupTree", "Ecal Barrel Imaging Size Lookup Tree");
    EcalBarrelImagingSizeLookupTree->Branch("p", &p_val, "p_val/D");
    EcalBarrelImagingSizeLookupTree->Branch("size", &size_array);
    EcalBarrelImagingSizeLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        size_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < sizeValues.size(); j++) {
            size_array.push_back(sizeValues.at(j));
            prob_array.push_back(SizeVsELookupTable[0][i][j]);
        }
        EcalBarrelImagingSizeLookupTree->Fill();
    }

    EcalBarrelImagingSizeLookupTree->Write();

    TTree *EcalBarrelScFiSizeLookupTree = new TTree("EcalBarrelScFiSizeLookupTree", "Ecal Barrel ScFi Size Lookup Tree");
    EcalBarrelScFiSizeLookupTree->Branch("p", &p_val, "p_val/D");
    EcalBarrelScFiSizeLookupTree->Branch("size", &size_array);
    EcalBarrelScFiSizeLookupTree->Branch("probability", &prob_array);   

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        size_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < sizeValues.size(); j++) {
            size_array.push_back(sizeValues.at(j));
            prob_array.push_back(SizeVsELookupTable[1][i][j]);
        }
        EcalBarrelScFiSizeLookupTree->Fill();
    }

    EcalBarrelScFiSizeLookupTree->Write();

    TTree *EcalEndcapPSizeLookupTree = new TTree("EcalEndcapPSizeLookupTree", "Ecal Endcap P Size Lookup Tree");
    EcalEndcapPSizeLookupTree->Branch("p", &p_val, "p_val/D");
    EcalEndcapPSizeLookupTree->Branch("size", &size_array);
    EcalEndcapPSizeLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        size_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < sizeValues.size(); j++) {
            size_array.push_back(sizeValues.at(j));
            prob_array.push_back(SizeVsELookupTable[2][i][j]);
        }
        EcalEndcapPSizeLookupTree->Fill();
    }

    EcalEndcapPSizeLookupTree->Write();

    TTree *EcalEndcapNSizeLookupTree = new TTree("EcalEndcapNSizeLookupTree", "Ecal Endcap N Size Lookup Tree");
    EcalEndcapNSizeLookupTree->Branch("p", &p_val, "p_val/D");
    EcalEndcapNSizeLookupTree->Branch("size", &size_array);
    EcalEndcapNSizeLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        size_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < sizeValues.size(); j++) {
            size_array.push_back(sizeValues.at(j));
            prob_array.push_back(SizeVsELookupTable[3][i][j]);
        }
        EcalEndcapNSizeLookupTree->Fill();
    }

    EcalEndcapNSizeLookupTree->Write();

    TTree *HcalBarrelSizeLookupTree = new TTree("HcalBarrelSizeLookupTree", "Hcal Barrel Size Lookup Tree");
    HcalBarrelSizeLookupTree->Branch("p", &p_val, "p_val/D");
    HcalBarrelSizeLookupTree->Branch("size", &size_array);
    HcalBarrelSizeLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        size_array.clear();
        prob_array.clear(); 
        for (size_t j = 1; j < sizeValues.size(); j++) {
            size_array.push_back(sizeValues.at(j));
            prob_array.push_back(SizeVsELookupTable[4][i][j]);
        }
        HcalBarrelSizeLookupTree->Fill();
    }

    HcalBarrelSizeLookupTree->Write();

    TTree *HcalEndcapPSizeLookupTree = new TTree("HcalEndcapPSizeLookupTree", "Hcal Endcap P Size Lookup Tree");
    HcalEndcapPSizeLookupTree->Branch("p", &p_val, "p_val/D");
    HcalEndcapPSizeLookupTree->Branch("size", &size_array);
    HcalEndcapPSizeLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        size_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < sizeValues.size(); j++) {
            size_array.push_back(sizeValues.at(j));
            prob_array.push_back(SizeVsELookupTable[5][i][j]);
        }
        HcalEndcapPSizeLookupTree->Fill();
    }

    HcalEndcapPSizeLookupTree->Write();

    TTree *LFHcalSizeLookupTree = new TTree("LFHcalSizeLookupTree", "LFHcal Size Lookup Tree");
    LFHcalSizeLookupTree->Branch("p", &p_val, "p_val/D");
    LFHcalSizeLookupTree->Branch("size", &size_array);
    LFHcalSizeLookupTree->Branch("probability", &prob_array);   

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        size_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < sizeValues.size(); j++) {
            size_array.push_back(sizeValues.at(j));
            prob_array.push_back(SizeVsELookupTable[6][i][j]);
        }
        LFHcalSizeLookupTree->Fill();
    }

    LFHcalSizeLookupTree->Write();

    TTree *HcalEndcapNSizeLookupTree = new TTree("HcalEndcapNSizeLookupTree", "Hcal Endcap N Size Lookup Tree");
    HcalEndcapNSizeLookupTree->Branch("p", &p_val, "p_val/D");
    HcalEndcapNSizeLookupTree->Branch("size", &size_array);
    HcalEndcapNSizeLookupTree->Branch("probability", &prob_array);

    for (size_t i = 1; i < pValues.size(); i++) {
        p_val = pValues.at(i);
        size_array.clear();
        prob_array.clear();
        for (size_t j = 1; j < sizeValues.size(); j++) {
            size_array.push_back(sizeValues.at(j));
            prob_array.push_back(SizeVsELookupTable[7][i][j]);
        }
        HcalEndcapNSizeLookupTree->Fill();
    }

    HcalEndcapNSizeLookupTree->Write();

    ofile->Close();


}