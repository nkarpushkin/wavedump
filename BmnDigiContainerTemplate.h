/* Copyright (C) 2021 Institute for Nuclear Research, Moscow
   SPDX-License-Identifier: GPL-3.0-only
   Authors: Nikolay Karpushkin [committer] */

/** \file BmnDigiContainerTemplate.h
 ** \author Nikolay Karpushkin <karpushkin@inr.ru>
 ** \date 18.01.2022
 **/

/** \class BmnDigiContainerTemplate
 ** \brief Data class for Bmn digi container template
 ** \version 1.0
 **/

#ifndef BmnDigiContainerTemplate_H
#define BmnDigiContainerTemplate_H 1

#include "TROOT.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <Rtypes.h>      // for THashConsistencyHolder, ClassDefNV

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include <stdint.h>
#include <numeric>
#include <string>  // for string


class BmnDigiContainerTemplate : public TObject {

public:
  /**@brief Default constructor.
       **/
  BmnDigiContainerTemplate();


  /**  Copy constructor **/
  BmnDigiContainerTemplate(const BmnDigiContainerTemplate&);


  /** Move constructor  **/
  BmnDigiContainerTemplate(BmnDigiContainerTemplate&&);


  /** Assignment operator  **/
  BmnDigiContainerTemplate& operator=(const BmnDigiContainerTemplate&);


  /** Move Assignment operator  **/
  BmnDigiContainerTemplate& operator=(BmnDigiContainerTemplate&&);


  /** Destructor **/
  virtual ~BmnDigiContainerTemplate()
  {
    std::vector<float>().swap(fWfm);
    std::vector<float>().swap(fFitWfm);
  }


  /** @brief Class name (static)
       ** @return BmnDigiContainerTemplate
       **/
  virtual const char* GetClassName() { return "BmnDigiContainerTemplate"; }


  /** @brief Fit R2 quality
       ** @return Fit R2 quality
       **/
  float GetFitR2() const { return fFitR2; };


  /** @brief Waveform
       ** @return Signal Waveform
       **/
  std::vector<float> GetWfm() const { return fWfm; }


  void reset();

  void DrawWfm();
  void DeleteCanvases() { gROOT->GetListOfCanvases()->Delete(); }

  float fAmpl           = 0;  /// Amplitude from waveform [adc counts]
  float fZL             = 0;  /// ZeroLevel from waveform [adc counts]
  float fIntegral       = 0;  /// Energy deposition from waveform [adc counts]
  float fTimeMax        = 0;  /// Time of maximum in waveform [adc samples]
  float fTimeCross      = 0;  /// Time cross threshold [adc samples]
  float fTimeCFD      = 0;  /// Time cross threshold [adc samples]
  float fToT            = 0;  /// Time over threshold [adc samples]

  float fFitAmpl      = 0.;  /// Amplitude from fit of waveform [adc counts]
  float fFitZL        = 0.;  /// ZeroLevel from fit of waveform [adc counts]
  float fFitIntegral  = 0.;  /// Energy deposition from fit of waveform [adc counts]
  float fFitR2        = 2.;  /// Quality of waveform fit [] -- good near 0
  float fFitTimeMax   = -1.; /// Time of maximum in fit of waveform [adc samples]

  std::vector<float> fWfm;
  std::vector<float> fFitWfm;

  template<class Archive>
  void serialize(Archive& ar, const unsigned int /*version*/)
  {
    ar& fAmpl;
    ar& fZL;
    ar& fIntegral;
    ar& fTimeMax;
    ar& fTimeCross;
    ar& fTimeCFD;
    ar& fToT;

    ar& fFitAmpl;
    ar& fFitZL;
    ar& fFitIntegral;
    ar& fFitR2;
    ar& fFitTimeMax;

    ar& fWfm;
    ar& fFitWfm;
  }

private:
  /// BOOST serialization interface
  friend class boost::serialization::access;


  ClassDefNV(BmnDigiContainerTemplate, 3);
};

#endif  // BmnDigiContainerTemplate_H
