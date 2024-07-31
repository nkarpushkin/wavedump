/* Copyright (C) 2021 Institute for Nuclear Research, Moscow
 SPDX-License-Identifier: GPL-3.0-only
 Authors: Nikolay Karpushkin [committer] */

 /** @file BmnDigiContainerTemplate.cxx
  ** @author Nikolay Karpushkin <karpushkin@inr.ru>
  ** @date 18.01.2022
  **
  ** Code for Data class for Bmn digi container template
  **/

#include "BmnDigiContainerTemplate.h"

#include <TBuffer.h> // for TBuffer
#include <TClass.h> // for BmnDigiContainerTemplate::IsA()
#include <TString.h> // for Form, TString

#include <string> // for basic_string

  // --- Default constructor
BmnDigiContainerTemplate::BmnDigiContainerTemplate()
  : TObject()
  , fAmpl()
  , fZL()
  , fIntegral()
  , fTimeMax()
  , fTimeCross()
  , fTimeCFD()
  , fToT()

  , fFitAmpl()
  , fFitZL()
  , fFitIntegral()
  , fFitR2()
  , fFitTimeMax()

  , fWfm()
  , fFitWfm() {}


// --- Copy constructor
BmnDigiContainerTemplate::BmnDigiContainerTemplate(const BmnDigiContainerTemplate& other)
  : TObject()
  , fAmpl(other.fAmpl)
  , fZL(other.fZL)
  , fIntegral(other.fIntegral)
  , fTimeMax(other.fTimeMax)
  , fTimeCross(other.fTimeCross)
  , fTimeCFD(other.fTimeCFD)
  , fToT(other.fToT)

  , fFitAmpl(other.fFitAmpl)
  , fFitZL(other.fFitZL)
  , fFitIntegral(other.fFitIntegral)
  , fFitR2(other.fFitR2)
  , fFitTimeMax(other.fFitTimeMax)

  , fWfm(other.fWfm)
  , fFitWfm(other.fFitWfm) {}


// --- Move constructor
BmnDigiContainerTemplate::BmnDigiContainerTemplate(BmnDigiContainerTemplate&& other)
  : TObject()
  , fAmpl(other.fAmpl)
  , fZL(other.fZL)
  , fIntegral(other.fIntegral)
  , fTimeMax(other.fTimeMax)
  , fTimeCross(other.fTimeCross)
  , fTimeCFD(other.fTimeCFD)
  , fToT(other.fToT)

  , fFitAmpl(other.fFitAmpl)
  , fFitZL(other.fFitZL)
  , fFitIntegral(other.fFitIntegral)
  , fFitR2(other.fFitR2)
  , fFitTimeMax(other.fFitTimeMax)

  , fWfm(other.fWfm)
  , fFitWfm(other.fFitWfm) {}


// --- Assignment operator
BmnDigiContainerTemplate& BmnDigiContainerTemplate::operator=(const BmnDigiContainerTemplate& other) {
  if (this != &other) {
    fAmpl = other.fAmpl;
    fZL = other.fZL;
    fIntegral = other.fIntegral;
    fTimeMax = other.fTimeMax;
    fTimeCross = other.fTimeCross;
    fTimeCFD = other.fTimeCFD;
    fToT = other.fToT;

    fFitAmpl = other.fFitAmpl;
    fFitZL = other.fFitZL;
    fFitIntegral = other.fFitIntegral;
    fFitR2 = other.fFitR2;
    fFitTimeMax = other.fFitTimeMax;

    fWfm = other.fWfm;
    fFitWfm = other.fFitWfm;
  }
  return *this;
}


// --- Move assignment operator
BmnDigiContainerTemplate& BmnDigiContainerTemplate::operator=(BmnDigiContainerTemplate&& other) {
  if (this != &other) {
    fAmpl = other.fAmpl;
    fZL = other.fZL;
    fIntegral = other.fIntegral;
    fTimeMax = other.fTimeMax;
    fTimeCross = other.fTimeCross;
    fTimeCFD = other.fTimeCFD;
    fToT = other.fToT;

    fFitAmpl = other.fFitAmpl;
    fFitZL = other.fFitZL;
    fFitIntegral = other.fFitIntegral;
    fFitR2 = other.fFitR2;
    fFitTimeMax = other.fFitTimeMax;

    fWfm = other.fWfm;
    fFitWfm = other.fFitWfm;
  }
  return *this;
}

void BmnDigiContainerTemplate::reset() {
  fAmpl = 0; /// Amplitude from waveform [adc counts]
  fZL = 0; /// ZeroLevel from waveform [adc counts]
  fIntegral = 0; /// Energy deposition from waveform [adc counts]
  fTimeMax = 0; /// Time of maximum in waveform [adc samples]
  fTimeCross = 0;
  fTimeCFD = 0;
  fToT = 0; /// Time over threshold [adc samples]

  fFitAmpl = 0.; /// Amplitude from fit of waveform [adc counts]
  fFitZL = 0.; /// ZeroLevel from fit of waveform [adc counts]
  fFitIntegral = 0.; /// Energy deposition from fit of waveform [adc counts]
  fFitR2 = 2.; /// Quality of waveform fit [] -- good near 0
  fFitTimeMax = -1.; /// Time of maximum in fit of waveform [adc samples]

  fWfm.clear();
  fFitWfm.clear();
}

void BmnDigiContainerTemplate::DrawWfm() {
    if (fWfm.empty()) return;

    // Create a TCanvas
    TCanvas* canv_ptr = new TCanvas();

    // Create a TMultiGraph to hold multiple TGraphs
    TMultiGraph* mg_ptr = new TMultiGraph();
    mg_ptr->SetTitle("Waveform and Fit");

    // Create a TGraph for the waveform
    std::vector<float> points(fWfm.size());
    std::iota(std::begin(points), std::end(points), 0); // Fill with 0, 1, ..., wfm.size() - 1.
    TGraph* tgr_ptr = new TGraph(fWfm.size(), &points[0], &fWfm[0]);
    mg_ptr->Add(tgr_ptr);
    tgr_ptr->SetTitle("Waveform");
    tgr_ptr->SetLineColor(kBlack);

    // Create a TGraph for the fit if available
    if (!fFitWfm.empty()) {
        TGraph* tgr_ptr_fit = new TGraph(fFitWfm.size(), &points[0], &fFitWfm[0]);
        mg_ptr->Add(tgr_ptr_fit);
        tgr_ptr_fit->SetLineColor(kRed);
        tgr_ptr_fit->SetLineWidth(2);
        tgr_ptr_fit->SetTitle("Fit");
    }

    // Draw the TMultiGraph
    mg_ptr->Draw("AL");
    mg_ptr->GetXaxis()->SetTitle("Time");
    mg_ptr->GetYaxis()->SetTitle("Amplitude");
    mg_ptr->GetYaxis()->SetTitleOffset(1.2);

    // Update canvas
    canv_ptr->Update();
}

ClassImp(BmnDigiContainerTemplate)
