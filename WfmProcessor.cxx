#include "WfmProcessor.h"
#include <functional> 

void WfmProcessor::ProcessWfm(std::vector<float> &wfm, DigiData *digi)
{
  if(isDebug) printf("Signal Type %d\n", digi->type);
  if(digi->type == SignalType::Blank) return;

  if (fdigiPars.isWriteWfm) digi->container.fWfm = wfm;
  assert(digi->gates.first > 0 && digi->gates.second > 0);
  unsigned int gateBegin = digi->gates.first;
  unsigned int gateEnd = digi->gates.second;
  
  if (gateBegin >= wfm.size()) {
    if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". waveform too short: accessing " <<
      gateBegin << "/" << wfm.size() << ". Check calibration file ";
    gateBegin = wfm.size() - 1;
  }
  if (gateEnd >= wfm.size()) {
    if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". waveform too short: accessing " <<
      gateEnd << "/" << wfm.size() << ". Check calibration file ";
    gateEnd = wfm.size() - 1;
  }

  //Zero level calculation
  if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". ZL calc\n";
  float ZL = std::accumulate(wfm.begin(), wfm.begin() + gateBegin, 0.0) / gateBegin;
  std::transform(wfm.begin(), wfm.end(), wfm.begin(), [ZL](float value) { return value - ZL; });
  digi->container.fZL = ZL;

  auto type_copy = digi->type;
  switch (type_copy) {
  case SignalType::LogicNegative:
  case SignalType::AnalogNegative:
  {
    // Invert
    if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". Inverting\n";
    std::transform(wfm.begin(), wfm.end(), wfm.begin(),
      std::bind(std::multiplies<float>(), -1.0, std::placeholders::_1));
    if (type_copy == SignalType::LogicNegative) {
      type_copy = SignalType::LogicPositive;
      if(isDebug) printf("Signal Type is switched to %d\n", type_copy);
    }
    if (type_copy == SignalType::AnalogNegative) {
      type_copy = SignalType::AnalogPositive;
      if(isDebug) printf("Signal Type is switched to %d\n", type_copy);
    }
    // Fall through intentionally
  }

  case SignalType::LogicPositive:
  case SignalType::AnalogPositive:
  {
    //MAX and Integral calculation including borders
    if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". MAX & INT search\n";
    digi->container.fIntegral = std::accumulate(wfm.begin() + gateBegin, wfm.begin() + gateEnd + 1, 0.0);
    auto const max_iter = std::max_element(wfm.begin() + gateBegin, wfm.begin() + gateEnd + 1);
    digi->container.fAmpl = *max_iter;
    digi->container.fTimeMax = (int)std::distance(wfm.begin(), max_iter);

    // TOT calculation for LogicPositive
    if (type_copy == SignalType::LogicPositive) {
      float searchAmpl = digi->container.fAmpl * 0.5;
      auto cross1 = GoToLevel(wfm, searchAmpl, digi->container.fTimeMax, -1);
      auto cross2 = GoToLevel(wfm, searchAmpl, digi->container.fTimeMax, 1);
      digi->container.fToT = cross2 - cross1;
      digi->container.fTimeCross = cross1;
      digi->container.fTimeCFD = cross1;
    }

    // TOT calculation for AnalogPositive
    if (type_copy == SignalType::AnalogPositive) {
      float searchAmpl = 20./0.488; // 20 mV 
      auto cross1 = GoToLevel(wfm, searchAmpl, digi->container.fTimeMax, -1);
      auto cross2 = GoToLevel(wfm, searchAmpl, digi->container.fTimeMax, 1);
      digi->container.fToT = cross2 - cross1;
      digi->container.fTimeCross = cross1;

      searchAmpl = digi->container.fAmpl * 0.2; // 20% 
      auto cross = GoToLevel(wfm, searchAmpl, digi->container.fTimeMax, -1);
      digi->container.fTimeCFD = cross;
    }

    //Prony fitting procedure
    if (type_copy == SignalType::AnalogPositive && digi->fitflag) {
      if (isDebug) std::cout << "WfmProcessor : Filling " << digi->container.GetClassName() << ". Fitting\n";

      //##################################
      const int model_order = 5;
      const int exponents = 3;
      //##################################

      PsdSignalFitting::PronyFitter Pfitter(model_order, exponents, gateBegin, gateEnd);
      //Pfitter.SetDebugMode(1);
      Pfitter.SetWaveform(wfm, 0.0);
      int best_signal_begin = Pfitter.ChooseBestSignalBeginHarmonics(digi->container.fTimeCFD-5, digi->container.fTimeCFD);

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

      digi->container.fFitIntegral = Pfitter.GetIntegral(gateBegin, gateEnd);
      digi->container.fFitAmpl = Pfitter.GetMaxAmplitude() - Pfitter.GetZeroLevel();
      float fit_R2 = Pfitter.GetRSquare(gateBegin, gateEnd);
      digi->container.fFitR2 = (fit_R2 > 2.0) ? 2.0 : fit_R2;
      digi->container.fFitZL = Pfitter.GetZeroLevel();
      digi->container.fFitTimeMax = Pfitter.GetSignalMaxTime();
      if (false) printf("fit integral %.0f integral %.0f R2 %.3f\n", digi->container.fFitIntegral, digi->container.fIntegral, digi->container.fFitR2);

      auto fit = Pfitter.GetFitWfm();
      if (digi->type == SignalType::AnalogPositive && fdigiPars.isWriteWfm) {
        std::transform(fit.begin(), fit.end(), fit.begin(), [ZL](float value) { return value + ZL; });
        digi->container.fFitWfm = fit;
      }
      if (digi->type == SignalType::AnalogNegative && fdigiPars.isWriteWfm) {
        std::transform(fit.begin(), fit.end(), fit.begin(),
          std::bind(std::multiplies<float>(), -1.0, std::placeholders::_1));
        std::transform(fit.begin(), fit.end(), fit.begin(), [ZL](float value) { return value + ZL; });
        digi->container.fFitWfm = fit;
      }
    }
    break;
  }

  default:
  {
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
