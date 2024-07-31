#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include "RawDataHandler.hpp"
//#include "BmnDigiContainerTemplate.cxx"
#include "WfmProcessor.cxx"

class DataProcessor : public WfmProcessor {


public:

    DataProcessor(std::vector<RawDataHandler::FileData>& dataFiles, std::pair<int,int> gate);
    ~DataProcessor();
    bool process(std::vector<RawDataHandler::FileData>& dataFiles);
    void setOutPath(std::string outPath);
    void setSignalTypes(std::map<string, DigiData::SignalType> sigTypes);
    void setFitFlags(std::map<string, bool> fitFlags);
    const std::map<string, DigiData> getOutDigi() {return outDigi;}
    // Add other processing methods...

private:
  std::map<string, DigiData> outDigi;
  TFile* outfile;
  TTree* outtree;

};

bool DataProcessor::process(std::vector<RawDataHandler::FileData>& dataFiles) {
  for(auto &it : outDigi) {
    it.second.container.reset();
  }

  bool flag = false;
  for (auto& file : dataFiles) {
    auto name = file.name;
    auto& wfm = file.data;

    auto ThisDigi = &outDigi.at(name);
    ProcessWfm(wfm, ThisDigi);
    if(ThisDigi->fitflag && ThisDigi->container.GetFitR2() < 3)
      flag = true;
  }

  outtree->Fill();
  return flag;
}

void DataProcessor::setOutPath(std::string outPath) {
  outfile = new TFile(outPath.c_str(), "RECREATE");
  outtree = new TTree("data", "data");
  for(auto &it : outDigi)
    outtree->Branch(it.first.c_str(), &it.second);
}

void DataProcessor::setSignalTypes(std::map<string, DigiData::SignalType> sigTypes) {
  for(auto &it : outDigi) {
    if (sigTypes.find(it.first) != sigTypes.end())
      it.second.type = sigTypes[it.first];    
  }
}

void DataProcessor::setFitFlags(std::map<string, bool> fitFlags) {
  for(auto &it : outDigi) {
    if (fitFlags.find(it.first) != fitFlags.end())
      it.second.fitflag = fitFlags[it.first];    
  }
}

DataProcessor::DataProcessor(std::vector<RawDataHandler::FileData>& dataFiles, std::pair<int,int> gate) : WfmProcessor() {
  DigiData empty;
  for(auto &file : dataFiles)
    outDigi[file.name] = empty;

  fdigiPars.gateBegin = gate.first;
  fdigiPars.gateEnd = gate.second;
}

DataProcessor::~DataProcessor() {
  if(outtree) {
    outtree->Write();
  }
  if(outfile) {
    outfile->Close();
    delete outfile;
  }
}
