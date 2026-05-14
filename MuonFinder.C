#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <numeric>
#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TVector3.h>
#include <TMath.h>
#include <TF1.h>

#include "DVMP_JPsi_Analysis.h"

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

double IsMuon(TVector3 trackMom, int simuID, EpLookupTree epTrees[8], SizeLookupTree sizeTrees[8],
            TTreeReaderArray<float>& EcalBarrelImgEng, TTreeReaderArray<float>& EcalBarrelScFiEng, TTreeReaderArray<float>& EcalEndcapPEng, TTreeReaderArray<float>& EcalEndcapNEng, 
            TTreeReaderArray<float>& HcalBarrelEng, TTreeReaderArray<float>& HcalEndcapPEng, TTreeReaderArray<float>& LFHcalEng, TTreeReaderArray<float>& HcalEndcapNEng, 
            TTreeReaderArray<int>& simuAssocEcalBarrelImg, TTreeReaderArray<int>& simuAssocEcalBarrelScFi, TTreeReaderArray<int>& simuAssocEcalEndcapP, TTreeReaderArray<int>& simuAssocEcalEndcapN,
            TTreeReaderArray<int>& simuAssocHcalBarrel, TTreeReaderArray<int>& simuAssocHcalEndcapP, TTreeReaderArray<int>& simuAssocLFHcal, TTreeReaderArray<int>& simuAssocHcalEndcapN)
{

    double EndcapNHcal_Factor = 1.0/6.0;

    std::string detectorNames[8] = {"EcalBarrelImaging", "EcalBarrelScFi", "EcalEndcapP", "EcalEndcapN", "HcalBarrel", "HcalEndcapP", "LFHcal", "HcalEndcapN"};

    double pMin = epTrees[0].tree->GetMinimum("p");

    double pMax = epTrees[0].tree->GetMaximum("p");

    int nRows = epTrees[0].tree->GetEntries();

    double energy[8], size[8];

    double particleMomentum = 0.;
    double particleEta = 0.;


    if (TMath::Abs(particleEta) > 4) return 0.0; // Outside of calorimeter acceptance
    if (TMath::Abs(particleEta) >= 1.0 && TMath::Abs(particleEta) <= 1.3) return 0.0; // Remove overlap region between barrel and endcap calorimeters
    if (particleMomentum < 0.05 || particleMomentum > 100.0) return 0.0;


    if (EcalBarrelImgEng.GetSize() == simuAssocEcalBarrelImg.GetSize())
    {
        for (int jI = 0; jI < simuAssocEcalBarrelImg.GetSize(); jI++) // Look for associations in the Ecal Barrel Imaging
        {
            if (simuAssocEcalBarrelImg[jI] == simuID)
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
            if (simuAssocEcalBarrelScFi[jS] == simuID)
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
            if (simuAssocEcalEndcapP[jP] == simuID)
            {
                energy[2] += EcalBarrelScFiEng[jP];
                size[2] += 1.;
            }
        }
    }
    if (EcalEndcapNEng.GetSize() == simuAssocEcalEndcapN.GetSize())
    {
        for (int jN = 0; jN < simuAssocEcalEndcapN.GetSize(); jN++) // Look for associations in the Ecal Endcap N
        {
            if (simuAssocEcalEndcapN[jN] == simuID)
            {
                energy[3] += EcalEndcapNEng[jN];
                size[3] += 1.;
            }
        }
    }

    if (HcalBarrelEng.GetSize() == simuAssocHcalBarrel.GetSize())
    {
        for (int jB = 0; jB < simuAssocHcalBarrel.GetSize(); jB++) // Look for associations in the Ecal Barrel
        {
            if (simuAssocHcalBarrel[jB] == simuID)
            {
                energy[4] += HcalBarrelEng[jB];
                size[4] += 1.;
            }
        }
    }
            
    if (HcalEndcapPEng.GetSize() == simuAssocHcalEndcapP.GetSize())
    {
        for (int jP = 0; jP < simuAssocHcalEndcapP.GetSize(); jP++) // Look for associations in the Ecal Endcap P
        {
            if (simuAssocHcalEndcapP[jP] == simuID)
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
            if (simuAssocLFHcal[jL] == simuID)
            {
                energy[6] += LFHcalEng[jL];
                size[6] += 1.;
            }
        }
    }

            
    if (HcalEndcapNEng.GetSize() == simuAssocHcalEndcapN.GetSize())
    {
        for (int jN = 0; jN < simuAssocHcalEndcapN.GetSize(); jN++) // Look for associations in the Ecal Endcap N
        {
            if (simuAssocHcalEndcapN[jN] == simuID)
            {
                energy[7] += HcalEndcapNEng[jN];
                size[7] += 1.;
            }
        }
    }

    if (energy[0] == 0.0 && energy[1] == 0.0 && energy[2] == 0.0 && energy[3] == 0.0 && energy[4] == 0.0 && energy[5] == 0.0 && energy[6] == 0.0 && energy[7] == 0.0) return 0.0; // Skip if no calorimeter energy is associated

    double finalProbEp[8] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
    double finalProbSize[8] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};

    double minDeltaP = 999.0;

    int bestRow = 0;

    for (int detectorsI = 0; detectorsI < 8; detectorsI++)
    {
        if (energy[detectorsI] == 0.0) continue; // Skip if no energy is deposited in this detector
        for (int iEP = 0; iEP < nRows; iEP++) 
        {
            epTrees[detectorsI].tree->GetEntry(iEP);
            double deltaP = TMath::Abs(epTrees[detectorsI].p - particleMomentum);
                    
            if (deltaP < minDeltaP) 
            {
                minDeltaP = deltaP;
                bestRow = iEP;
            } 
            else 
            {
                break; 
            }
        }
        epTrees[detectorsI].tree->GetEntry(bestRow);
        double minDeltaEp = 999.0;
        int bestCol = 0;

        for (size_t jEP = 0; jEP < epTrees[detectorsI].Ep_ptr->size(); jEP++) 
        {
            double deltaEp = TMath::Abs(epTrees[detectorsI].Ep_ptr->at(jEP) - energy[detectorsI]/particleMomentum);
                    
            if (deltaEp < minDeltaEp) 
            {
                minDeltaEp = deltaEp;
                bestCol = jEP;
            } 
            else 
            {
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

        for (int iSize = 0; iSize < nRows; iSize++) 
        {
            sizeTrees[detectorsI].tree->GetEntry(iSize);
            double deltaP = TMath::Abs(sizeTrees[detectorsI].p - particleMomentum);
                    
            if (deltaP < minDeltaP) 
            {
                minDeltaP = deltaP;
                bestRow = iSize;
            } 
            else 
            {
                break; 
            }
        }

        sizeTrees[detectorsI].tree->GetEntry(bestRow);
        minDeltaEp = 999.0;
        bestCol = 0;

        for (size_t jSize = 0; jSize < sizeTrees[detectorsI].size_ptr->size(); jSize++) 
        {
            double deltaEp = TMath::Abs(sizeTrees[detectorsI].size_ptr->at(jSize) - size[detectorsI]);
                    
            if (deltaEp < minDeltaEp) 
            {
                minDeltaEp = deltaEp;
                bestCol = jSize;
            } 
            else 
            {
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

        for (int k = 0; k < 8; k++) 
        {
            if (energy[k] > 0) 
            {
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

        return logL_muon;

}