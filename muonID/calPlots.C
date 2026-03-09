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

void calPlots()
{

    TString infile="reconOut/mu_20GeV_reconOut.root";
    std::string outfilename = "outputs/muCalPlotsOutput.root";

    // Set output file for the histograms
    TFile *ofile = TFile::Open(outfilename.c_str(),"RECREATE");

    // Set up input file chain
    TChain *mychain = new TChain("events");
    mychain->Add(infile);
    //mychain->Add(infileMuPlus);

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
    TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");
    TTreeReaderArray<int> trackPDG(tree_reader, "ReconstructedChargedParticles.PDG");
    TTreeReaderArray<float> trackMass(tree_reader, "ReconstructedChargedParticles.mass");
    TTreeReaderArray<float> trackCharge(tree_reader, "ReconstructedChargedParticles.charge");
    TTreeReaderArray<float> trackEng(tree_reader, "ReconstructedChargedParticles.energy");

    // Get Truth Seeded Reconstructed Track Information
    TTreeReaderArray<float> ttrackMomX(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x");
    TTreeReaderArray<float> ttrackMomY(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y");
    TTreeReaderArray<float> ttrackMomZ(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z");
    TTreeReaderArray<int> ttrackPDG(tree_reader, "ReconstructedTruthSeededChargedParticles.PDG");
    TTreeReaderArray<float> ttrackMass(tree_reader, "ReconstructedTruthSeededChargedParticles.mass");
    TTreeReaderArray<float> ttrackCharge(tree_reader, "ReconstructedTruthSeededChargedParticles.charge");
    TTreeReaderArray<float> ttrackEng(tree_reader, "ReconstructedTruthSeededChargedParticles.energy");

    // Get Associations Between MCParticles and ReconstructedChargedParticles
    TTreeReaderArray<unsigned int> recoAssoc(tree_reader, "ReconstructedChargedParticleAssociations.recID");
    TTreeReaderArray<unsigned int> simuAssoc(tree_reader, "ReconstructedChargedParticleAssociations.simID");

    // Get B0 Information
    TTreeReaderArray<unsigned int> recoAssocB0(tree_reader, "B0ECalClusterAssociations.recID");
    TTreeReaderArray<unsigned int> simuAssocB0(tree_reader, "B0ECalClusterAssociations.simID");
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
    TTreeReaderArray<unsigned int> simuAssocEcalBarrel(tree_reader, "EcalBarrelClusterAssociations.simID");
    TTreeReaderArray<unsigned int> recoAssocEcalBarrel(tree_reader, "EcalBarrelClusterAssociations.recID");
    TTreeReaderArray<float> EcalBarrelEng(tree_reader, "EcalBarrelClusters.energy");

    TTreeReaderArray<float> EcalBarrelImagingHitsEng(tree_reader, "EcalBarrelImagingRecHits.energy");
    TTreeReaderArray<float> EcalBarrelImagingHitsX(tree_reader, "EcalBarrelImagingRecHits.position.x");
    TTreeReaderArray<float> EcalBarrelImagingHitsY(tree_reader, "EcalBarrelImagingRecHits.position.y");
    TTreeReaderArray<float> EcalBarrelImagingHitsZ(tree_reader, "EcalBarrelImagingRecHits.position.z");

    TTreeReaderArray<float> EcalBarrelScFiHitsEng(tree_reader, "EcalBarrelScFiRecHits.energy");
    TTreeReaderArray<float> EcalBarrelScFiHitsX(tree_reader, "EcalBarrelScFiRecHits.position.x");
    TTreeReaderArray<float> EcalBarrelScFiHitsY(tree_reader, "EcalBarrelScFiRecHits.position.y");
    TTreeReaderArray<float> EcalBarrelScFiHitsZ(tree_reader, "EcalBarrelScFiRecHits.position.z");

    TTreeReaderArray<unsigned int> simuAssocEcalEndcapP(tree_reader, "EcalEndcapPClusterAssociations.simID");
    TTreeReaderArray<unsigned int> recoAssocEcalEndcapP(tree_reader, "EcalEndcapPClusterAssociations.recID");    
    TTreeReaderArray<float> EcalEndcapPEng(tree_reader, "EcalEndcapPClusters.energy");

    TTreeReaderArray<float> EcalEndcapPHitsEng(tree_reader, "EcalEndcapPRecHits.energy");
    TTreeReaderArray<float> EcalEndcapPHitsX(tree_reader, "EcalEndcapPRecHits.position.x");
    TTreeReaderArray<float> EcalEndcapPHitsY(tree_reader, "EcalEndcapPRecHits.position.y");
    TTreeReaderArray<float> EcalEndcapPHitsZ(tree_reader, "EcalEndcapPRecHits.position.z");

    TTreeReaderArray<unsigned int> simuAssocEcalEndcapN(tree_reader, "EcalEndcapNClusterAssociations.simID");
    TTreeReaderArray<unsigned int> recoAssocEcalEndcapN(tree_reader, "EcalEndcapNClusterAssociations.recID");
    TTreeReaderArray<float> EcalEndcapNEng(tree_reader, "EcalEndcapNClusters.energy");

    TTreeReaderArray<float> EcalEndcapNHitsEng(tree_reader, "EcalEndcapNRecHits.energy");
    TTreeReaderArray<float> EcalEndcapNHitsX(tree_reader, "EcalEndcapNRecHits.position.x");
    TTreeReaderArray<float> EcalEndcapNHitsY(tree_reader, "EcalEndcapNRecHits.position.y");
    TTreeReaderArray<float> EcalEndcapNHitsZ(tree_reader, "EcalEndcapNRecHits.position.z");

    // Hcal Information
    TTreeReaderArray<unsigned int> simuAssocHcalBarrel(tree_reader, "HcalBarrelClusterAssociations.simID");
    TTreeReaderArray<unsigned int> recoAssocHcalBarrel(tree_reader, "HcalBarrelClusterAssociations.recID");
    TTreeReaderArray<float> HcalBarrelEng(tree_reader, "HcalBarrelClusters.energy");

    TTreeReaderArray<float> HcalBarrelHitsEng(tree_reader, "HcalBarrelRecHits.energy");
    TTreeReaderArray<float> HcalBarrelHitsX(tree_reader, "HcalBarrelRecHits.position.x");
    TTreeReaderArray<float> HcalBarrelHitsY(tree_reader, "HcalBarrelRecHits.position.y");
    TTreeReaderArray<float> HcalBarrelHitsZ(tree_reader, "HcalBarrelRecHits.position.z");

    TTreeReaderArray<unsigned int> simuAssocHcalEndcapP(tree_reader, "HcalEndcapPInsertClusterAssociations.simID");
    TTreeReaderArray<unsigned int> recoAssocHcalEndcapP(tree_reader, "HcalEndcapPInsertClusterAssociations.recID");    
    TTreeReaderArray<float> HcalEndcapPEng(tree_reader, "HcalEndcapPInsertClusters.energy");

    TTreeReaderArray<float> HcalEndcapPHitsEng(tree_reader, "HcalEndcapPInsertRecHits.energy");
    TTreeReaderArray<float> HcalEndcapPHitsX(tree_reader, "HcalEndcapPInsertRecHits.position.x");
    TTreeReaderArray<float> HcalEndcapPHitsY(tree_reader, "HcalEndcapPInsertRecHits.position.y");
    TTreeReaderArray<float> HcalEndcapPHitsZ(tree_reader, "HcalEndcapPInsertRecHits.position.z");

    TTreeReaderArray<unsigned int> simuAssocLFHcal(tree_reader, "LFHCALClusterAssociations.simID");
    TTreeReaderArray<unsigned int> recoAssocLFHcal(tree_reader, "LFHCALClusterAssociations.recID");    
    TTreeReaderArray<float> LFHcalEng(tree_reader, "LFHCALClusters.energy");

    TTreeReaderArray<float> HcalLFHitsEng(tree_reader, "LFHCALRecHits.energy");
    TTreeReaderArray<float> HcalLFHitsX(tree_reader, "LFHCALRecHits.position.x");
    TTreeReaderArray<float> HcalLFHitsY(tree_reader, "LFHCALRecHits.position.y");
    TTreeReaderArray<float> HcalLFHitsZ(tree_reader, "LFHCALRecHits.position.z");

    TTreeReaderArray<unsigned int> simuAssocHcalEndcapN(tree_reader, "HcalEndcapNClusterAssociations.simID");
    TTreeReaderArray<unsigned int> recoAssocHcalEndcapN(tree_reader, "HcalEndcapNClusterAssociations.recID");
    TTreeReaderArray<float> HcalEndcapNEng(tree_reader, "HcalEndcapNClusters.energy");

    TTreeReaderArray<float> HcalEndcapNHitsEng(tree_reader, "HcalEndcapNRecHits.energy");
    TTreeReaderArray<float> HcalEndcapNHitsX(tree_reader, "HcalEndcapNRecHits.position.x");
    TTreeReaderArray<float> HcalEndcapNHitsY(tree_reader, "HcalEndcapNRecHits.position.y");
    TTreeReaderArray<float> HcalEndcapNHitsZ(tree_reader, "HcalEndcapNRecHits.position.z");




    // Define histograms

    TH2D *h_EcalHitsXY = new TH2D("h_EcalHitsXY","Ecal Hits XY; x (cm); y (cm)",125,-2500.,2500.,125,-2500.,2500.);
    TH2D *h_EcalHitsZY = new TH2D("h_EcalHitsZY","Ecal Hits ZY; z (cm); y (cm)",150,-3000.,3000.,125,-2500.,2500.);
    TH2D *h_EcalHitsZX = new TH2D("h_EcalHitsZX","Ecal Hits ZX; z (cm); x (cm)",150,-3000.,3000.,125,-2500.,2500.);

    TH2D *h_EcalHitsXYevent = new TH2D("h_EcalHitsXYevent","Ecal Hits XY; x (cm); y (cm)",125,-2500.,2500.,125,-2500.,2500.);
    TH2D *h_EcalHitsZYevent = new TH2D("h_EcalHitsZYevent","Ecal Hits ZY; z (cm); y (cm)",150,-3000.,3000.,125,-2500.,2500.);
    TH2D *h_EcalHitsZXevent = new TH2D("h_EcalHitsZXevent","Ecal Hits ZX; z (cm); x (cm)",150,-3000.,3000.,125,-2500.,2500.);

    TH2D *h_HcalHitsXY = new TH2D("h_HcalHitsXY","Hcal Hits XY; x (cm); y (cm)",125,-2500.,2500.,125,-2500.,2500.);
    TH2D *h_HcalHitsZY = new TH2D("h_HcalHitsZY","Hcal Hits ZY; z (cm); y (cm)",150,-3000.,3000.,125,-2500.,2500.);
    TH2D *h_HcalHitsZX = new TH2D("h_HcalHitsZX","Hcal Hits ZX; z (cm); x (cm)",150,-3000.,3000.,125,-2500.,2500.);

    TH2D *h_HcalHitsXYevent = new TH2D("h_HcalHitsXYevent","Hcal Hits XY; x (cm); y (cm)",125,-2500.,2500.,125,-2500.,2500.);
    TH2D *h_HcalHitsZYevent = new TH2D("h_HcalHitsZYevent","Hcal Hits ZY; z (cm); y (cm)",150,-3000.,3000.,125,-2500.,2500.);
    TH2D *h_HcalHitsZXevent = new TH2D("h_HcalHitsZXevent","Hcal Hits ZX; z (cm); x (cm)",150,-3000.,3000.,125,-2500.,2500.);

    TH1D *h_EcalBarrelImagingHitsDeviation = new TH1D("h_EcalBarrelImagingHitsDeviation","Ecal Barrel Imaging Hits Deviation; Deviation (cm)",125,0.,5000.);
    TH1D *h_EcalBarrelScFiHitsDeviation = new TH1D("h_EcalBarrelScFiHitsDeviation","Ecal Barrel ScFi Hits Deviation; Deviation (cm)",125,0.,5000.);
    TH1D *h_EcalEndcapPHitsDeviation = new TH1D("h_EcalEndcapPHitsDeviation","Ecal Endcap P Hits Deviation; Deviation (cm)",125,0.,5000.);
    TH1D *h_EcalEndcapNHitsDeviation = new TH1D("h_EcalEndcapNHitsDeviation","Ecal Endcap N Hits Deviation; Deviation (cm)",125,0.,5000.);

    TH1D *h_HcalBarrelHitsDeviation = new TH1D("h_HcalBarrelHitsDeviation","Hcal Barrel Hits Deviation; Deviation (cm)",125,0.,5000.);
    TH1D *h_HcalEndcapPHitsDeviation = new TH1D("h_HcalEndcapPHitsDeviation","Hcal Endcap P Hits Deviation; Deviation (cm)",125,0.,5000.);
    TH1D *h_HcalEndcapNHitsDeviation = new TH1D("h_HcalEndcapNHitsDeviation","Hcal Endcap N Hits Deviation; Deviation (cm)",125,0.,5000.);
    TH1D *h_HcalLFHitsDeviation = new TH1D("h_HcalLFHitsDeviation","LFHcal Hits Deviation; Deviation (cm)",125,0.,5000.);

    int eventID = 0;
    double maxEcalImgDeviation = 0.;
    double maxEcalScFiDeviation = 0.;
    double maxEcalPEndcapDeviation = 0.;
    double maxEcalNEndcapDeviation = 0.;
    double maxHcalDeviation = 0.;
    double maxHcalPEndcapDeviation = 0.;
    double maxHcalNEndcapDeviation = 0.;
    double maxHcalLFDeviation = 0.;

    while(tree_reader.Next()) // Loop over events
    {
      eventID++;

      maxEcalImgDeviation = 0.;
      maxEcalScFiDeviation = 0.;
      maxEcalPEndcapDeviation = 0.;
      maxEcalNEndcapDeviation = 0.;
      maxHcalDeviation = 0.;
      maxHcalPEndcapDeviation = 0.;
      maxHcalNEndcapDeviation = 0.;
      maxHcalLFDeviation = 0.;

      h_EcalHitsXYevent->Reset();
      h_EcalHitsZYevent->Reset();
      h_EcalHitsZXevent->Reset(); 
      h_HcalHitsXYevent->Reset();
      h_HcalHitsZYevent->Reset();
      h_HcalHitsZXevent->Reset();

      if (eventID % 10 == 0)
      {
        fprintf (stderr, "%4.2f Percent\r ", eventID*100.0/mychain->GetEntries());
        fflush (stderr);
      }

      for (int iEH = 0; iEH < EcalBarrelImagingHitsEng.GetSize(); iEH++)
      {
        h_EcalHitsXY->Fill(EcalBarrelImagingHitsX[iEH],EcalBarrelImagingHitsY[iEH],EcalBarrelImagingHitsEng[iEH]);
        h_EcalHitsZY->Fill(EcalBarrelImagingHitsZ[iEH],EcalBarrelImagingHitsY[iEH],EcalBarrelImagingHitsEng[iEH]);
        h_EcalHitsZX->Fill(EcalBarrelImagingHitsZ[iEH],EcalBarrelImagingHitsX[iEH],EcalBarrelImagingHitsEng[iEH]);
        h_EcalHitsXYevent->Fill(EcalBarrelImagingHitsX[iEH],EcalBarrelImagingHitsY[iEH],EcalBarrelImagingHitsEng[iEH]);
        h_EcalHitsZYevent->Fill(EcalBarrelImagingHitsZ[iEH],EcalBarrelImagingHitsY[iEH],EcalBarrelImagingHitsEng[iEH]);
        h_EcalHitsZXevent->Fill(EcalBarrelImagingHitsZ[iEH],EcalBarrelImagingHitsX[iEH],EcalBarrelImagingHitsEng[iEH]);

        for (int iEH2 = 0; iEH2 < EcalBarrelImagingHitsEng.GetSize(); iEH2++)
        {
            if (iEH == iEH2) continue;
            double deviation = sqrt(pow(EcalBarrelImagingHitsX[iEH]-EcalBarrelImagingHitsX[iEH2],2) + pow(EcalBarrelImagingHitsY[iEH]-EcalBarrelImagingHitsY[iEH2],2) + pow(EcalBarrelImagingHitsZ[iEH]-EcalBarrelImagingHitsZ[iEH2],2));
            if (deviation > maxEcalImgDeviation) maxEcalImgDeviation = deviation;
        }

      }

      if (maxEcalImgDeviation > 0.) h_EcalBarrelImagingHitsDeviation->Fill(maxEcalImgDeviation);

      for (int iEH = 0; iEH < EcalBarrelScFiHitsEng.GetSize(); iEH++)
      {
        h_EcalHitsXY->Fill(EcalBarrelScFiHitsX[iEH],EcalBarrelScFiHitsY[iEH],EcalBarrelScFiHitsEng[iEH]);
        h_EcalHitsZY->Fill(EcalBarrelScFiHitsZ[iEH],EcalBarrelScFiHitsY[iEH],EcalBarrelScFiHitsEng[iEH]);
        h_EcalHitsZX->Fill(EcalBarrelScFiHitsZ[iEH],EcalBarrelScFiHitsX[iEH],EcalBarrelScFiHitsEng[iEH]);
        h_EcalHitsXYevent->Fill(EcalBarrelScFiHitsX[iEH],EcalBarrelScFiHitsY[iEH],EcalBarrelScFiHitsEng[iEH]);
        h_EcalHitsZYevent->Fill(EcalBarrelScFiHitsZ[iEH],EcalBarrelScFiHitsY[iEH],EcalBarrelScFiHitsEng[iEH]);
        h_EcalHitsZXevent->Fill(EcalBarrelScFiHitsZ[iEH],EcalBarrelScFiHitsX[iEH],EcalBarrelScFiHitsEng[iEH]);
        
        for (int iEH2 = 0; iEH2 < EcalBarrelScFiHitsEng.GetSize(); iEH2++)
        {
            if (iEH == iEH2) continue;
            double deviation = sqrt(pow(EcalBarrelScFiHitsX[iEH]-EcalBarrelScFiHitsX[iEH2],2) + pow(EcalBarrelScFiHitsY[iEH]-EcalBarrelScFiHitsY[iEH2],2) + pow(EcalBarrelScFiHitsZ[iEH]-EcalBarrelScFiHitsZ[iEH2],2));
            if (deviation > maxEcalScFiDeviation) maxEcalScFiDeviation = deviation;
        }

      }

      if (maxEcalScFiDeviation > 0.) h_EcalBarrelScFiHitsDeviation->Fill(maxEcalScFiDeviation);

      for (int iEH = 0; iEH < EcalEndcapPHitsEng.GetSize(); iEH++)
      {
        h_EcalHitsXY->Fill(EcalEndcapPHitsX[iEH],EcalEndcapPHitsY[iEH],EcalEndcapPHitsEng[iEH]);
        h_EcalHitsZY->Fill(EcalEndcapPHitsZ[iEH],EcalEndcapPHitsY[iEH],EcalEndcapPHitsEng[iEH]);
        h_EcalHitsZX->Fill(EcalEndcapPHitsZ[iEH],EcalEndcapPHitsX[iEH],EcalEndcapPHitsEng[iEH]);
        h_EcalHitsXYevent->Fill(EcalEndcapPHitsX[iEH],EcalEndcapPHitsY[iEH],EcalEndcapPHitsEng[iEH]);
        h_EcalHitsZYevent->Fill(EcalEndcapPHitsZ[iEH],EcalEndcapPHitsY[iEH],EcalEndcapPHitsEng[iEH]);
        h_EcalHitsZXevent->Fill(EcalEndcapPHitsZ[iEH],EcalEndcapPHitsX[iEH],EcalEndcapPHitsEng[iEH]);
        
        for (int iEH2 = 0; iEH2 < EcalEndcapPHitsEng.GetSize(); iEH2++)
        {
            if (iEH == iEH2) continue;
            double deviation = sqrt(pow(EcalEndcapPHitsX[iEH]-EcalEndcapPHitsX[iEH2],2) + pow(EcalEndcapPHitsY[iEH]-EcalEndcapPHitsY[iEH2],2) + pow(EcalEndcapPHitsZ[iEH]-EcalEndcapPHitsZ[iEH2],2));
            if (deviation > maxEcalPEndcapDeviation) maxEcalPEndcapDeviation = deviation;
        }

      }

      if (maxEcalPEndcapDeviation > 0.) h_EcalEndcapPHitsDeviation->Fill(maxEcalPEndcapDeviation);

      for (int iEH = 0; iEH < EcalEndcapNHitsEng.GetSize(); iEH++)
      {
        h_EcalHitsXY->Fill(EcalEndcapNHitsX[iEH],EcalEndcapNHitsY[iEH],EcalEndcapNHitsEng[iEH]);
        h_EcalHitsZY->Fill(EcalEndcapNHitsZ[iEH],EcalEndcapNHitsY[iEH],EcalEndcapNHitsEng[iEH]);
        h_EcalHitsZX->Fill(EcalEndcapNHitsZ[iEH],EcalEndcapNHitsX[iEH],EcalEndcapNHitsEng[iEH]);
        h_EcalHitsXYevent->Fill(EcalEndcapNHitsX[iEH],EcalEndcapNHitsY[iEH],EcalEndcapNHitsEng[iEH]);
        h_EcalHitsZYevent->Fill(EcalEndcapNHitsZ[iEH],EcalEndcapNHitsY[iEH],EcalEndcapNHitsEng[iEH]);
        h_EcalHitsZXevent->Fill(EcalEndcapNHitsZ[iEH],EcalEndcapNHitsX[iEH],EcalEndcapNHitsEng[iEH]);
        
        for (int iEH2 = 0; iEH2 < EcalEndcapNHitsEng.GetSize(); iEH2++)
        {
            if (iEH == iEH2) continue;
            double deviation = sqrt(pow(EcalEndcapNHitsX[iEH]-EcalEndcapNHitsX[iEH2],2) + pow(EcalEndcapNHitsY[iEH]-EcalEndcapNHitsY[iEH2],2) + pow(EcalEndcapNHitsZ[iEH]-EcalEndcapNHitsZ[iEH2],2));
            if (deviation > maxEcalNEndcapDeviation) maxEcalNEndcapDeviation = deviation;
        }

      }

      if (maxEcalNEndcapDeviation > 0.) h_EcalEndcapNHitsDeviation->Fill(maxEcalNEndcapDeviation);

      for (int iHH = 0; iHH < HcalBarrelHitsEng.GetSize(); iHH++)
      {
        h_HcalHitsXY->Fill(HcalBarrelHitsX[iHH],HcalBarrelHitsY[iHH],HcalBarrelHitsEng[iHH]);
        h_HcalHitsZY->Fill(HcalBarrelHitsZ[iHH],HcalBarrelHitsY[iHH],HcalBarrelHitsEng[iHH]);
        h_HcalHitsZX->Fill(HcalBarrelHitsZ[iHH],HcalBarrelHitsX[iHH],HcalBarrelHitsEng[iHH]);
        h_HcalHitsXYevent->Fill(HcalBarrelHitsX[iHH],HcalBarrelHitsY[iHH],HcalBarrelHitsEng[iHH]);
        h_HcalHitsZYevent->Fill(HcalBarrelHitsZ[iHH],HcalBarrelHitsY[iHH],HcalBarrelHitsEng[iHH]);
        h_HcalHitsZXevent->Fill(HcalBarrelHitsZ[iHH],HcalBarrelHitsX[iHH],HcalBarrelHitsEng[iHH]);

        for (int iHH2 = 0; iHH2 < HcalBarrelHitsEng.GetSize(); iHH2++)
        {
            if (iHH == iHH2) continue;
            double deviation = sqrt(pow(HcalBarrelHitsX[iHH]-HcalBarrelHitsX[iHH2],2) + pow(HcalBarrelHitsY[iHH]-HcalBarrelHitsY[iHH2],2) + pow(HcalBarrelHitsZ[iHH]-HcalBarrelHitsZ[iHH2],2));
            if (deviation > maxHcalDeviation) maxHcalDeviation = deviation;
        }

      }


      if (maxHcalDeviation > 0.) h_HcalBarrelHitsDeviation->Fill(maxHcalDeviation);

      for (int iHH = 0; iHH < HcalEndcapPHitsEng.GetSize(); iHH++)
      {
        h_HcalHitsXY->Fill(HcalEndcapPHitsX[iHH],HcalEndcapPHitsY[iHH],HcalEndcapPHitsEng[iHH]);
        h_HcalHitsZY->Fill(HcalEndcapPHitsZ[iHH],HcalEndcapPHitsY[iHH],HcalEndcapPHitsEng[iHH]);
        h_HcalHitsZX->Fill(HcalEndcapPHitsZ[iHH],HcalEndcapPHitsX[iHH],HcalEndcapPHitsEng[iHH]);
        h_HcalHitsXYevent->Fill(HcalEndcapPHitsX[iHH],HcalEndcapPHitsY[iHH],HcalEndcapPHitsEng[iHH]);
        h_HcalHitsZYevent->Fill(HcalEndcapPHitsZ[iHH],HcalEndcapPHitsY[iHH],HcalEndcapPHitsEng[iHH]);
        h_HcalHitsZXevent->Fill(HcalEndcapPHitsZ[iHH],HcalEndcapPHitsX[iHH],HcalEndcapPHitsEng[iHH]);

        for (int iHH2 = 0; iHH2 < HcalEndcapPHitsEng.GetSize(); iHH2++)
        {
            if (iHH == iHH2) continue;
            double deviation = sqrt(pow(HcalEndcapPHitsX[iHH]-HcalEndcapPHitsX[iHH2],2) + pow(HcalEndcapPHitsY[iHH]-HcalEndcapPHitsY[iHH2],2) + pow(HcalEndcapPHitsZ[iHH]-HcalEndcapPHitsZ[iHH2],2));
            if (deviation > maxHcalPEndcapDeviation) maxHcalPEndcapDeviation = deviation;
        }

      }

      if (maxHcalPEndcapDeviation > 0.) h_HcalEndcapPHitsDeviation->Fill(maxHcalPEndcapDeviation);

      for (int iHH = 0; iHH < HcalEndcapNHitsEng.GetSize(); iHH++)
      {
        h_HcalHitsXY->Fill(HcalEndcapNHitsX[iHH],HcalEndcapNHitsY[iHH],HcalEndcapNHitsEng[iHH]);
        h_HcalHitsZY->Fill(HcalEndcapNHitsZ[iHH],HcalEndcapNHitsY[iHH],HcalEndcapNHitsEng[iHH]);
        h_HcalHitsZX->Fill(HcalEndcapNHitsZ[iHH],HcalEndcapNHitsX[iHH],HcalEndcapNHitsEng[iHH]);
        h_HcalHitsXYevent->Fill(HcalEndcapNHitsX[iHH],HcalEndcapNHitsY[iHH],HcalEndcapNHitsEng[iHH]);
        h_HcalHitsZYevent->Fill(HcalEndcapNHitsZ[iHH],HcalEndcapNHitsY[iHH],HcalEndcapNHitsEng[iHH]);
        h_HcalHitsZXevent->Fill(HcalEndcapNHitsZ[iHH],HcalEndcapNHitsX[iHH],HcalEndcapNHitsEng[iHH]);

        for (int iHH2 = 0; iHH2 < HcalEndcapNHitsEng.GetSize(); iHH2++)
        {
            if (iHH == iHH2) continue;
            double deviation = sqrt(pow(HcalEndcapNHitsX[iHH]-HcalEndcapNHitsX[iHH2],2) + pow(HcalEndcapNHitsY[iHH]-HcalEndcapNHitsY[iHH2],2) + pow(HcalEndcapNHitsZ[iHH]-HcalEndcapNHitsZ[iHH2],2));
            if (deviation > maxHcalNEndcapDeviation) maxHcalNEndcapDeviation = deviation;
        }

      }

      if (maxHcalNEndcapDeviation > 0.) h_HcalEndcapNHitsDeviation->Fill(maxHcalNEndcapDeviation);

      for (int iHH = 0; iHH < HcalLFHitsEng.GetSize(); iHH++)
      {
        h_HcalHitsXY->Fill(HcalLFHitsX[iHH],HcalLFHitsY[iHH],HcalLFHitsEng[iHH]);
        h_HcalHitsZY->Fill(HcalLFHitsZ[iHH],HcalLFHitsY[iHH],HcalLFHitsEng[iHH]);
        h_HcalHitsZX->Fill(HcalLFHitsZ[iHH],HcalLFHitsX[iHH],HcalLFHitsEng[iHH]);
        h_HcalHitsXYevent->Fill(HcalLFHitsX[iHH],HcalLFHitsY[iHH],HcalLFHitsEng[iHH]);
        h_HcalHitsZYevent->Fill(HcalLFHitsZ[iHH],HcalLFHitsY[iHH],HcalLFHitsEng[iHH]);
        h_HcalHitsZXevent->Fill(HcalLFHitsZ[iHH],HcalLFHitsX[iHH],HcalLFHitsEng[iHH]);

        for (int iHH2 = 0; iHH2 < HcalLFHitsEng.GetSize(); iHH2++)
        {
            if (iHH == iHH2) continue;
            double deviation = sqrt(pow(HcalLFHitsX[iHH]-HcalLFHitsX[iHH2],2) + pow(HcalLFHitsY[iHH]-HcalLFHitsY[iHH2],2) + pow(HcalLFHitsZ[iHH]-HcalLFHitsZ[iHH2],2));
            if (deviation > maxHcalLFDeviation) maxHcalLFDeviation = deviation;
        }

      }

      if (maxHcalLFDeviation > 0.) h_HcalLFHitsDeviation->Fill(maxHcalLFDeviation);

      /*
      if(EcalBarrelHitsEng.GetSize() > 10 && HcalBarrelHitsEng.GetSize() > 10)
      {
        TCanvas *c1 = new TCanvas("c1","Ecal and Hcal Barrel Hits",1200,800);
        c1->Divide(2,3);
        c1->cd(1);
        h_EcalBarrelHitsXYevent->Draw("COLZ");
        c1->cd(2);
        h_HcalBarrelHitsXYevent->Draw("COLZ");
        c1->cd(3);
        h_EcalBarrelHitsZXevent->Draw("COLZ");
        c1->cd(4);
        h_HcalBarrelHitsZXevent->Draw("COLZ");
        c1->cd(5);
        h_EcalBarrelHitsZYevent->Draw("COLZ");
        c1->cd(6);
        h_HcalBarrelHitsZYevent->Draw("COLZ");
        c1->Update();
        std::cout << "Press Enter to continue..." << std::endl;
        std::cin.get();
        delete c1;
      }
      */

    }

    // Write output histograms to file
    ofile->cd();
    ofile->mkdir("calClusterShapePlots");
    ofile->cd("calClusterShapePlots");
    h_EcalHitsXY->Write();
    h_EcalHitsZY->Write();
    h_EcalHitsZX->Write();
    h_HcalHitsXY->Write();
    h_HcalHitsZY->Write();
    h_HcalHitsZX->Write();
    ofile->cd("..");
    ofile->mkdir("deviationPlots");
    ofile->cd("deviationPlots");
    h_EcalBarrelImagingHitsDeviation->Write();
    h_EcalBarrelScFiHitsDeviation->Write();
    h_EcalEndcapPHitsDeviation->Write();
    h_EcalEndcapNHitsDeviation->Write();
    h_HcalBarrelHitsDeviation->Write();
    h_HcalEndcapPHitsDeviation->Write();
    h_HcalEndcapNHitsDeviation->Write();
    h_HcalLFHitsDeviation->Write();
    ofile->cd("..");
    ofile->Close();
}