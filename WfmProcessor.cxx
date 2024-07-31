#include "WfmProcessor.h"
#include <functional> 

void WfmProcessor::ProcessWfm(std::vector<float> &wfm, DigiData *digi)
{
  auto type_copy = digi->type;
  switch (type_copy) {
  case SignalType::Blank:
  case SignalType::Trigger:
    // Common processing for Blank and Trigger cases
    break;

  case SignalType::LogicNegative:
  case SignalType::AnalogNegative:
  {
    // Invert
    if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". Inverting";
    std::transform(wfm.begin(), wfm.end(), wfm.begin(),
      std::bind(std::multiplies<float>(), -1.0, std::placeholders::_1));
    if (type_copy == SignalType::LogicNegative) 
      type_copy = SignalType::LogicPositive;
    if (type_copy == SignalType::AnalogNegative) 
      type_copy = SignalType::AnalogPositive;
    // Fall through intentionally
  }

  case SignalType::LogicPositive:
  {
    if (fdigiPars.isWriteWfm) digi->container.fWfm = wfm;
    assert(fdigiPars.gateBegin > 0 && fdigiPars.gateEnd > 0);
    if (fdigiPars.gateBegin >= wfm.size()) {
      if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". waveform too short: accessing " <<
        fdigiPars.gateBegin << "/" << wfm.size() << ". Check calibration file ";
      fdigiPars.gateBegin = wfm.size() - 1;
    }
    if (fdigiPars.gateEnd >= wfm.size()) {
      if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". waveform too short: accessing " <<
        fdigiPars.gateEnd << "/" << wfm.size() << ". Check calibration file ";
      fdigiPars.gateEnd = wfm.size() - 1;
    }
    //Zero level calculation
    if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". ZL calc";
    digi->container.fZL = std::accumulate(wfm.begin(), wfm.begin() + fdigiPars.gateBegin, 0.0) / fdigiPars.gateBegin;

    //MAX and Integral calculation including borders
    if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". MAX & INT search";
    digi->container.fIntegral = std::accumulate(wfm.begin() + fdigiPars.gateBegin, wfm.begin() + fdigiPars.gateEnd + 1,
      -digi->container.fZL * (fdigiPars.gateEnd - fdigiPars.gateBegin + 1));
    auto const max_iter = std::max_element(wfm.begin() + fdigiPars.gateBegin, wfm.begin() + fdigiPars.gateEnd + 1);
    digi->container.fAmpl = *max_iter - digi->container.fZL;
    digi->container.fTimeMax = (int)std::distance(wfm.begin(), max_iter);

    //TOT calculation
    float searchAmpl = digi->container.fZL + digi->container.fAmpl * 0.5;
    auto cross1 = GoToLevel(wfm, searchAmpl, digi->container.fTimeMax, -1);
    auto cross2 = GoToLevel(wfm, searchAmpl, digi->container.fTimeMax, 1);
    digi->container.fToT = cross2 - cross1;
    digi->container.fTimeCross = cross1;
    digi->container.fTimeCFD = cross1;
    break;
  }

  default:
  case SignalType::AnalogPositive:
  {
    if (fdigiPars.isWriteWfm) digi->container.fWfm = wfm;
    assert(fdigiPars.gateBegin > 0 && fdigiPars.gateEnd > 0);
    if (fdigiPars.gateBegin >= wfm.size()) {
      if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". waveform too short: accessing " <<
        fdigiPars.gateBegin << "/" << wfm.size() << ". Check calibration file ";
      fdigiPars.gateBegin = wfm.size() - 1;
    }
    if (fdigiPars.gateEnd >= wfm.size()) {
      if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". waveform too short: accessing " <<
        fdigiPars.gateEnd << "/" << wfm.size() << ". Check calibration file ";
      fdigiPars.gateEnd = wfm.size() - 1;
    }
    //Zero level calculation
    if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". ZL calc";
    digi->container.fZL = std::accumulate(wfm.begin(), wfm.begin() + fdigiPars.gateBegin, 0.0) / fdigiPars.gateBegin;

    //MAX and Integral calculation including borders
    if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". MAX & INT search";
    digi->container.fIntegral = std::accumulate(wfm.begin() + fdigiPars.gateBegin, wfm.begin() + fdigiPars.gateEnd + 1,
      -digi->container.fZL * (fdigiPars.gateEnd - fdigiPars.gateBegin + 1));
    auto const max_iter = std::max_element(wfm.begin() + fdigiPars.gateBegin, wfm.begin() + fdigiPars.gateEnd + 1);
    digi->container.fAmpl = *max_iter - digi->container.fZL;
    digi->container.fTimeMax = (int)std::distance(wfm.begin(), max_iter);

    //TOT calculation
    float searchAmpl = digi->container.fZL + 20./0.488; // 20 mV 
    auto cross1 = GoToLevel(wfm, searchAmpl, digi->container.fTimeMax, -1);
    auto cross2 = GoToLevel(wfm, searchAmpl, digi->container.fTimeMax, 1);
    digi->container.fToT = cross2 - cross1;
    digi->container.fTimeCross = cross1;

    searchAmpl = digi->container.fZL + digi->container.fAmpl * 0.2; // 20% 
    auto cross = GoToLevel(wfm, searchAmpl, digi->container.fTimeMax, -1);
    digi->container.fTimeCFD = cross;

    //Prony fitting procedure
    if (digi->fitflag) {
      if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". Fitting";

      //##################################
      const int model_order = 5;
      const int exponents = 3;
      //##################################

      PsdSignalFitting::PronyFitter Pfitter(model_order, exponents, fdigiPars.gateBegin, fdigiPars.gateEnd);
      //Pfitter.SetDebugMode(1);
      Pfitter.SetWaveform(wfm, digi->container.fZL);
      int best_signal_begin = Pfitter.ChooseBestSignalBeginHarmonics(fdigiPars.gateBegin, digi->container.fTimeMax);

      Pfitter.SetSignalBegin(best_signal_begin);
      Pfitter.CalculateFitHarmonics();
      Pfitter.CalculateFitAmplitudes();
      complex<float>* harmonics = Pfitter.GetHarmonics();
      complex<float>* amplituds = Pfitter.GetAmplitudes();
      std::vector<std::pair<float,float>> result;
      for(int i = 0; i<=model_order; ++i){
        result.push_back(std::make_pair(-1.0/log(real(harmonics[i])), real(amplituds[i])));
      }
      digi->tau = result;

      digi->container.fFitIntegral = Pfitter.GetIntegral(fdigiPars.gateBegin, fdigiPars.gateEnd);
      digi->container.fFitAmpl = Pfitter.GetMaxAmplitude() - Pfitter.GetZeroLevel();
      float fit_R2 = Pfitter.GetRSquare(fdigiPars.gateBegin, fdigiPars.gateEnd);
      digi->container.fFitR2 = (fit_R2 > 2.0) ? 2.0 : fit_R2;
      digi->container.fFitZL = Pfitter.GetZeroLevel();
      digi->container.fFitTimeMax = Pfitter.GetSignalMaxTime();

      if (false) printf("fit integral %.0f integral %.0f R2 %.3f\n", digi->container.fFitIntegral, digi->container.fIntegral, digi->container.fFitR2);
      if (fdigiPars.isWriteWfm)
        digi->container.fFitWfm = Pfitter.GetFitWfm();

    }
    break;
  }
  }


}

void WfmProcessor::MeanRMScalc(std::vector<float> wfm, float *Mean, float *RMS, int begin, int end, int step)
{
  begin = (begin < 0) ? 0 : begin;
  if (begin > end)
  {
    float swap = end;
    end = begin;
    begin = swap;
  };
  step = TMath::Abs(step);
  *Mean = *RMS = 0.;
  int Delta = 0;
  for (int n = begin; n <= end; n += step)
  {
    *Mean += wfm[n];
    Delta++;
  }
  *Mean /= (float)Delta;
  for (int n = begin; n <= end; n += step)
    *RMS += (wfm[n] - *Mean) * (wfm[n] - *Mean);
  *RMS = sqrt(*RMS / ((float)Delta));
  //printf("AMPL %.2f, RMS %.2f\n",*Mean,*RMS);
}

float WfmProcessor::GoToLevel(std::vector<float>& wfm, float Level, int point, int iterator) {
  float ResultTime;
  while (point >= 0 && point < (int)wfm.size() && point+iterator >= 0 && point+iterator < (int)wfm.size()) {
    if ((Level - wfm.at(point)) * (Level - wfm.at(point+iterator)) <= 0) {
      ResultTime = LevelBy2Points(point, wfm.at(point), point + iterator, wfm.at(point+iterator), Level);
      return ResultTime;
    }
    point += iterator;
  } //while

  point = -1;
  return 0;
}

float WfmProcessor::LevelBy2Points(float X1, float Y1, float X2, float Y2, float Y0)
{
    //            [X2, Y2] 0
    //                   *
    //                 *
    //   Y0--------- 0
    //             *  X0 (returned)
    //  [X1, Y1] 0
    return (X1*Y0 - X1*Y2 - X2*Y0 + X2*Y1) / (Y1-Y2);
}


ClassImp(WfmProcessor)
