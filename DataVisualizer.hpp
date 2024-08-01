#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <TCanvas.h>
#include <TH1F.h>
#include "RawDataHandler.hpp"
//#include "DataProcessor.hpp"

class DataVisualizer {
private:
  TCanvas* canvas;
  std::vector<TH1F*> graphs;
  std::vector<int> colors;
  std::vector<int> ignore;

public:
    DataVisualizer(int samples, std::pair<int,int> range);
    ~DataVisualizer();

    bool onStatus;
    void fillGraphs(std::vector<RawDataHandler::FileData>& dataFiles);
    void fillFit(const std::map<string, DigiData> outDigi, std::map<string, bool> fitFlags);
    int visualize();
    bool getStatus();

};

void DataVisualizer::fillGraphs(std::vector<RawDataHandler::FileData>& dataFiles) {
  if (!canvas) {
    canvas = new TCanvas("c1", "c1", 800, 600);
    if (!canvas) {
      std::cerr << "Failed to create canvas." << std::endl;
      return;
    }
  }

  // Clear previous graphs and free memory
  for (auto& graph : graphs) {
    delete graph;
  }
  graphs.clear();

  for (auto &file : dataFiles) {
    auto channel = file.channel;
    if(std::find(ignore.begin(), ignore.end(), channel) != ignore.end()) 
      continue;
    auto& wfm = file.data;
    TH1F* gr = new TH1F(Form("channel_%d", channel), Form("channel_%d; Sample ; Amplitude", channel), wfm.size(), 0, wfm.size());
    int color = kBlack;
    if (channel < colors.size()) color = colors.at(channel);
    gr->SetLineColor(color);
    gr->SetMarkerColor(color);
    for (size_t j = 0; j < wfm.size(); ++j)
      gr->SetBinContent(j, wfm.at(j));
   
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.5);
    gr->SetStats(0);
    graphs.push_back(gr); // Store the graph pointer
  }

  return;
}

void DataVisualizer::fillFit(const std::map<string, DigiData> outDigi, std::map<string, bool> fitFlags) {
  for(auto &it : fitFlags) {
    if(it.second) {
      auto& wfm = outDigi.at(it.first).container.fFitWfm;
      TH1F* gr = new TH1F(Form("fit_%s", it.first.c_str()), Form("fit_%s; Sample ; Amplitude", it.first.c_str()), wfm.size(), 0, wfm.size());
      gr->SetLineColor(kRed);
      gr->SetLineWidth(2);
      gr->SetStats(0); 
      for (size_t j = 0; j < wfm.size(); ++j)
        gr->SetBinContent(j, wfm.at(j));
      graphs.push_back(gr); // Store the graph pointer
    }
  }
}

int DataVisualizer::visualize() {
  for(auto& gr : graphs){
    TString name = gr->GetName();
    if(name.Contains("fit")) 
      gr->Draw("L same");
    else
      gr->Draw("ALP same");
  }

  //canvas->Modified();
  canvas->Update();
  string user;
  std::getline(std::cin, user);
  if (user.empty()) 
    return 0;
  if (user == "c") {
    onStatus = false;
    return 0;
  }
  if (user == "q") {
    onStatus = false;
    return 1;
  }
  try {
    int number = std::stoi(user);
    auto it = std::find(ignore.begin(), ignore.end(), number);
    if (it != ignore.end())
      ignore.erase(it);
    else
      ignore.push_back(number);
  }
  catch (const std::exception& e) {
    std::cerr << "Failed to convert string to integer: " << e.what() << std::endl;
    return 0;
  }

  return 0;
}

bool DataVisualizer::getStatus() {
  return onStatus;
}

DataVisualizer::DataVisualizer(int samples, std::pair<int,int> range) {
  onStatus = true;
  canvas = new TCanvas("c1", "c1", 800, 600);
  canvas->SetGrid();
  TH1F* tmp = new TH1F("tmp", "", samples, 0, samples);
  tmp->GetYaxis()->SetRangeUser(range.first, range.second); 
  tmp->SetStats(0);
  tmp->Draw();
  colors = {kRed, kBlue, kMagenta, kOrange, kGreen+2, kYellow+2, kCyan+2, kViolet+2, kBlack, kRed-5, kBlue-5, kMagenta-5, kOrange-5, kGreen-5, kYellow-5, kCyan-5, kViolet-5, kBlack-5};
}
DataVisualizer::~DataVisualizer() {
}
