#!/usr/bin/python 
## @file readAncFiles.py Build histograms and covariance matrix from supplementary files
## Manoa Andriamirado, IIT 2023

import sys
import os
import ROOT
import numpy as np
import matplotlib.pyplot as plt

from math import sqrt, pow

baselineEdges = np.array([6.65826,7.35108,7.49064,7.67506,7.78471,7.96415,8.06882,8.18346,8.30308,8.52737,9.2501])

def buildCovMatrix(path):
    inputCovMatrix = os.path.join(path, "1.2_Osc_CovarianceMatrix.txt")
    print("    - Building covariance matrix from {}".format(inputCovMatrix))

    try:
        covMatrixData = np.loadtxt(inputCovMatrix, delimiter=',')
    except OSError:
        print("    ERROR: Could not find {}".format(inputCovMatrix))
        sys.exit()

    return covMatrixData

def extractHist(h1D):
    xZ = []

    nXs = h1D.GetNbinsX()
    
    for nX in range(1, nXs + 1):
        x = h1D.GetXaxis().GetBinCenter(nX)
        Z = h1D.GetBinContent(nX)

        xZ.append([x, Z])
    
    return xZ


def writeToTxtFile(hName, txtName):
    hist = rootFile.Get(hName)

    if hist == None:
        print("Could not find " + str(hName))
        sys.exit(2)

    typeHistogram = 1 if type(hist) is ROOT.TH1D else 2

    if typeHistogram == 2:
        nYs = hist.GetNbinsY()

        for nY in range(1, nYs):
            with open(txtName + str(nY) + ".txt", "w") as file_object:
                for x, Z in extractHist(hist.ProjectionX("temp", nY, nY)):
                    file_object.write("{0:.1f},{1:.3f}\n".format(x, Z))

    elif typeHistogram == 1:
        with open(txtName + ".txt", "w") as file_object:
            for x, Z in extractHist(hist):
                file_object.write("{0:.1f},{1:.3f}\n".format(x, Z))


def buildPrompt(path):
    print("    - Building prompt signals")  
    
    segBinning = np.linspace(-0.5, 153.5, 155, np.double)

    hLvsERef   = ROOT.TH2D("LvsERef", "Energy vs Baseline", 16, 0.8, 7.2, 10, baselineEdges.astype(np.double))
    hLvsERefEr = ROOT.TH2D("LvsERefEr", "Energy vs Baseline", 16, 0.8, 7.2, 10, baselineEdges.astype(np.double))

    segmentIndexMap = SegmentToBaselineMapping(path)

    for segNumber in range(0, 154):
        # Check if fiducialized segment
        segX = segNumber % 14
        segZ = segNumber // 14

        ## 11 x 14 segmented detector -> Remove outer segments
        if (segX == 0 or segX == 13):
            continue
        if (segZ == 0 or segZ == 10):
            continue

        inputPromptSeg = os.path.join(path, "1.4_Osc_Prompt{}.txt".format(segNumber))
        try:
            promptSegData = np.loadtxt(inputPromptSeg, dtype=np.double, delimiter=',')
        except OSError:
            print("    WARNING: Could not find {}".format(inputPromptSeg))
            continue
        
        baselineIndex = segmentIndexMap[segNumber]
        
        for eBinCenter, SigBinContent, ErrBinContent in zip(promptSegData[:,0].tolist(), promptSegData[:,1].tolist(), promptSegData[:,2].tolist()):
            lBinCenter = hLvsERef.GetYaxis().GetBinCenter(baselineIndex)

            hLvsERef.Fill(eBinCenter, lBinCenter, SigBinContent)
            hLvsERefEr.Fill(eBinCenter, lBinCenter, ErrBinContent * ErrBinContent)
        
    # Sum squared error
    for l in range(1, hLvsERef.GetNbinsY()+1):
        for e in range(1, hLvsERef.GetNbinsX()+1):
            hLvsERefEr.SetBinContent(e, l, sqrt(hLvsERefEr.GetBinContent(e, l)))

    return (hLvsERef, hLvsERefEr)
        

def buildNullOscillation(path):
    print("    - Building null oscillation")

    hLvsENull = ROOT.TH2D("LvsENull", "Energy vs Baseline", 16, 0.8, 7.2, 10, baselineEdges.astype(np.double, order='C'))

    for baselineIndex in range(1, hLvsENull.GetNbinsY()+1):
        inputNullSeg = os.path.join(path, "1.6_Osc_NullOscPred{}.txt".format(baselineIndex))
        try:
            nullSegData = np.loadtxt(inputNullSeg, delimiter=',')
        except OSError:
            print("  WARNING: Could not find {}".format(inputNullSeg))
            continue
        
        for energyIndex, binContent in enumerate(nullSegData[:,1].tolist(), start=1):
            hLvsENull.SetBinContent(energyIndex, baselineIndex, binContent)

    return hLvsENull

def SegmentToBaselineMapping(path):
    inputFile = os.path.join(path, "1.1_Osc_SegmentMap.txt")
    try:
        mappingData = np.loadtxt(inputFile, delimiter=',')
    except OSError:
        print("    WARNING: Could not find {}".format(inputFile))
        sys.exit()

    # Construct a dictionary where the keys are the segment index
    # and value as the baseline index 
    segNumber = mappingData[:,2].astype(np.int).tolist()
    baselineIndex = mappingData[:,0].astype(np.int).tolist()
    segmentMapping = dict(zip(segNumber, baselineIndex))

    return segmentMapping

def RatioDataNullPlot(path):
    print("++ Creating Data/Prediction ratio plot ++")

    # Reconstruct data from supp files
    hLvsERef, hLvsERefErr = buildPrompt(path)
    hLvsENull = buildNullOscillation(path) 
    covMatrix = buildCovMatrix(path)

    # Relativize Prediction P_{e,l}*(M_e/P_e)
    for i in range(1, hLvsERef.GetNbinsX()+1):
        ratio = hLvsERef.Integral(i, i, 1, -1)/hLvsENull.Integral(i, i, 1, -1)

        for j in range(1, hLvsERef.GetNbinsY()+1):
            hLvsENull.SetBinContent(i, j, hLvsENull.GetBinContent(i, j)*ratio)
    
    listHist = []

    # Create histogram
    for h in range(1, hLvsERef.GetNbinsX()+1):
        hTmp = hLvsERef.ProjectionY("h"+str(h), 1, -1)
        minBin = hLvsERef.GetXaxis().GetBinLowEdge(h)
        maxBin = minBin + hLvsERef.GetXaxis().GetBinWidth(h)

        hTmp.SetTitle("[{0:.1f}, {1:.1f}[ MeV; L [m]; Data/Prediction".format(minBin, maxBin))
        listHist.append(hTmp)
    
    # Fill histogram
    for e in range(1, hLvsERef.GetNbinsX()+1):
        for l in range(1, hLvsERef.GetNbinsY()+1):
            binContentData = hLvsERef.GetBinContent(e, l)
            binContentNull = hLvsENull.GetBinContent(e, l)
            
            # covariance matrix index
            iCovMatrix = (l-1)*hLvsERef.GetNbinsX() + (e-1)
            errorData = pow(hLvsERefErr.GetBinContent(e, l)/binContentData, 2)
            errorNull = pow(sqrt(covMatrix[iCovMatrix][iCovMatrix])/binContentNull, 2)
            
            ratio = binContentData/binContentNull
            error = sqrt(errorData + errorNull)
            
            listHist[e-1].SetBinContent(l, ratio)
            listHist[e-1].SetBinError(l, error)
    
    # Create canvas
    can = ROOT.TCanvas("can", "can", 800, 800)
    ROOT.gStyle.SetOptStat(0)
    can.Divide(4, 4)
    
    for h in range(len(listHist)):
        can.cd(h+1)
        
        listHist[h].SetMarkerStyle(20)
        listHist[h].SetMarkerSize(0.7)
        #listHist[h].GetYaxis().SetRangeUser(0, 2)
        listHist[h].Draw("PE")

    can.SaveAs("RatioDataNull.pdf")

def Usage():
    print("ERROR!")
    print("Usage: ./ReadAncillaryFiles.py path/to/ancfiles")
    sys.exit()

if __name__ == "__main__":
    # Check argument
    if len(sys.argv) != 2:
        Usage()

    RatioDataNullPlot(sys.argv[1])

