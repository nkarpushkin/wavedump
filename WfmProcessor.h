#ifndef WfmProcessor_H
#define WfmProcessor_H

#include <iostream>
#include <numeric>
#include <vector>
#include <complex>

#include "DigiData.hpp"
#include "PronyFitter.cxx"

using SignalType = DigiData::SignalType;

class WfmProcessor {
public:

    WfmProcessor() { isDebug = false; }
    ~WfmProcessor() {}

    struct digiPars {
        bool isWriteWfm;
        unsigned int gateBegin;
        unsigned int gateEnd;
        float threshold;
    } fdigiPars;

    void ProcessWfm(std::vector<float> &wfm, DigiData *digi);
    bool isDebug;

private:
    void MeanRMScalc(std::vector<float> wfm, float *Mean, float *RMS, int begin, int end, int step = 1);
    float GoToLevel(std::vector<float> &wfm, float Level, int point, int iterator);
    float LevelBy2Points(float X1, float Y1, float X2, float Y2, float Y0);

    ClassDef(WfmProcessor, 1);
};

#endif /* WfmProcessor_H */
