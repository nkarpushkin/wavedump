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

    DataProcessor(std::vector<RawDataHandler::FileData>& dataFiles);
    ~DataProcessor();
    bool process(std::vector<RawDataHandler::FileData>& dataFiles);
    void setOutPath(std::string outPath);
    void setSignalTypes(std::map<string, DigiData::SignalType> sigTypes);
    void setCommonGates(std::pair<unsigned int, unsigned int> gate);
    void setGates(std::map<string, std::pair<unsigned int, unsigned int>> gates);
    void setCommonZLGates(std::pair<unsigned int, unsigned int> gate);
    void setZLGates(std::map<string, std::pair<unsigned int, unsigned int>> gates);
    void setFitFlags(std::map<string, bool> fitFlags);
    void setFitTau(std::map<std::string, std::vector<std::complex<float>>> fitTau);
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
    if (isDebug) printf("Working on %s\n", name.c_str());
    ProcessWfm(wfm, ThisDigi);
    if(ThisDigi->fitflag && ThisDigi->container.GetFitR2() < 3)
      flag = true;
  }

  outtree->Fill();
  return flag;
}

void DataProcessor::setOutPath(std::string outPath) {
  outfile = new TFile(outPath.c_str(), "RECREATE");
  outtree = new TTree("tree", "tree");
  for(auto &it : outDigi)
    outtree->Branch(it.first.c_str(), &it.second);
}

void DataProcessor::setSignalTypes(std::map<string, DigiData::SignalType> sigTypes) {
  for(auto &it : outDigi) {
    if (sigTypes.find(it.first) != sigTypes.end())
      it.second.type = sigTypes[it.first];    
  }
}

void DataProcessor::setCommonGates(std::pair<unsigned int, unsigned int> gate) {
  for(auto &it : outDigi)
    it.second.gates = gate;
}

void DataProcessor::setGates(std::map<string, std::pair<unsigned int, unsigned int>> gates) {
  for(auto &it : outDigi) {
    if (gates.find(it.first) != gates.end())
      it.second.gates = gates[it.first];    
  }
}

void DataProcessor::setCommonZLGates(std::pair<unsigned int, unsigned int> gate) {
  for(auto &it : outDigi)
    it.second.zlgates = gate;
}

void DataProcessor::setZLGates(std::map<string, std::pair<unsigned int, unsigned int>> gates) {
  for(auto &it : outDigi) {
    if (gates.find(it.first) != gates.end())
      it.second.zlgates = gates[it.first];    
  }
}

void DataProcessor::setFitFlags(std::map<string, bool> fitFlags) {
  for(auto &it : outDigi) {
    if (fitFlags.find(it.first) != fitFlags.end()) {
      it.second.fitflag = fitFlags[it.first];    
      if (isDebug) printf("Set %s to fit flag %d\n", it.first.c_str(), it.second.fitflag);
    }
  }
}

void DataProcessor::setFitTau(std::map<std::string, std::vector<std::complex<float>>> fitTau) {
  for (auto& it : outDigi) {
    auto found = fitTau.find(it.first);
    if (found != fitTau.end()) {
      it.second.fitflag = true;
      const auto& tauvec = found->second;
      it.second.tau.clear();
      it.second.tau.reserve(tauvec.size());
      for (const auto& c : tauvec) {
        it.second.tau.emplace_back(c.real(), c.imag());
      }
      if (isDebug) {
        printf("Set %s fit taus:", it.first.c_str());
        for (const auto& tau : it.second.tau)
          printf(" {%.0f,%.0f}", tau.first, tau.second);
        printf("\n");
      }
    }
  }
}


DataProcessor::DataProcessor(std::vector<RawDataHandler::FileData>& dataFiles) : WfmProcessor() {
  DigiData empty;
  for(auto &file : dataFiles)
    outDigi[file.name] = empty;
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
