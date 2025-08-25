#include "RawDataHandler.hpp"
#include "DataProcessor.hpp"
#include "DataVisualizer.hpp"
#include <dirent.h>

#include <fstream>
#include <cstring>
#include <cstddef>
#include <filesystem>

using SignalType = DigiData::SignalType;


std::vector<std::string> getFilesInFolder(const std::string& folderPath) {
  std::vector<std::string> files;
  DIR* dir;
  struct dirent* entry;

  // Open the directory
  dir = opendir(folderPath.c_str());
  if (dir == nullptr) {
    std::cerr << "Error opening directory: " << folderPath << std::endl;
    return files;
  }

  // Read directory entries
  while ((entry = readdir(dir)) != nullptr) {
    std::string filename = entry->d_name;
    if (filename.find(".dat") != std::string::npos) {
      std::cout << "Found dat file: " << filename << std::endl;
      files.push_back(filename);
    }
  }

  // Close the directory
  closedir(dir);

  return files;
}

int wavedump(const std::string& folderPath = "", const std::string& outPath = "result.root", int nEventsToProcess = -1, bool visualize = 0) {
  //constants here
  const int header_length = 6;
  const int wfm_length = 1024;
  const int visual_min_adc = 2000;
  const int visual_max_adc = 3000;
  const int gate_beg = 450;
  const int gate_end = 700;
  const int zl_beg = 400;
  const int zl_end = 440;

  RawDataHandler rawDataHandler;
  DataVisualizer dataVisualizer(wfm_length, std::make_pair(visual_min_adc, visual_max_adc));

  auto file_list = getFilesInFolder(folderPath);
  for (auto& filename : file_list)
    rawDataHandler.openFile(folderPath + "/" + filename);

  //auto size = std::filesystem::file_size((file_list + file_list.at(0).c_str()));
  //long nEvents = floor(size / wfm_length);
  //nEventsToProcess = (nEventsToProcess <= 0) ? nEvents : std::min(nEvents, nEventsToProcess);

  DataProcessor dataProcessor(rawDataHandler.dataFiles);
  dataProcessor.setCommonZLGates(std::make_pair(zl_beg, zl_end));
  dataProcessor.setCommonGates(std::make_pair(gate_beg, gate_end));
  dataProcessor.fdigiPars.isWriteWfm = false;
  dataProcessor.setOutPath(outPath);

  std::map<std::string, SignalType> sigTypes;
  // sigTypes["TR_0_0"] = SignalType::LogicNegative;
  // sigTypes["wave_0"] = SignalType::AnalogPositive;
  // sigTypes["wave_1"] = SignalType::AnalogPositive;
  // sigTypes["wave_2"] = SignalType::Blank;
  // sigTypes["wave_3"] = SignalType::LogicNegative;
  // sigTypes["wave_4"] = SignalType::AnalogPositive;
  // sigTypes["wave_5"] = SignalType::LogicNegative;
  // sigTypes["wave_6"] = SignalType::AnalogPositive;
  // sigTypes["wave_7"] = SignalType::Blank;
  sigTypes["TR_0_0"] = SignalType::AnalogNegative;
  sigTypes["TR_0_1"] = SignalType::AnalogNegative;
  sigTypes["wave_0"] = SignalType::AnalogNegative;
  sigTypes["wave_1"] = SignalType::AnalogNegative;
  sigTypes["wave_2"] = SignalType::AnalogNegative;
  sigTypes["wave_3"] = SignalType::AnalogNegative;
  sigTypes["wave_4"] = SignalType::AnalogNegative;
  sigTypes["wave_5"] = SignalType::AnalogNegative;
  sigTypes["wave_6"] = SignalType::AnalogNegative;
  sigTypes["wave_7"] = SignalType::AnalogNegative;
  sigTypes["wave_8"] = SignalType::AnalogNegative;
  sigTypes["wave_9"] = SignalType::AnalogNegative;
  sigTypes["wave_10"] = SignalType::AnalogNegative;
  sigTypes["wave_11"] = SignalType::AnalogNegative;
  sigTypes["wave_12"] = SignalType::AnalogNegative;
  sigTypes["wave_13"] = SignalType::AnalogNegative;
  sigTypes["wave_14"] = SignalType::AnalogNegative;
  sigTypes["wave_15"] = SignalType::AnalogNegative;
  dataProcessor.setSignalTypes(sigTypes);

  // std::map<std::string, bool> fitFlags;
  // fitFlags["TR_0_0"] = true;
  // fitFlags["TR_0_1"] = true;
  // fitFlags["wave_0"] = true;
  // fitFlags["wave_1"] = true;
  // fitFlags["wave_2"] = true;
  // fitFlags["wave_3"] = true;
  // fitFlags["wave_4"] = true;
  // fitFlags["wave_5"] = true;
  // fitFlags["wave_6"] = true;
  // fitFlags["wave_7"] = true;
  // fitFlags["wave_8"] = true;
  // fitFlags["wave_9"] = true;
  // fitFlags["wave_10"] = true;
  // fitFlags["wave_11"] = true;
  // fitFlags["wave_12"] = true;
  // fitFlags["wave_13"] = true;
  // fitFlags["wave_14"] = true;
  // fitFlags["wave_15"] = true;
  // dataProcessor.setFitFlags(fitFlags);

  std::map<std::string, std::vector<std::complex<float>>> fitTau;
  fitTau["TR_0_0"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["TR_0_1"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_0"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_1"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_2"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_3"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_4"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_5"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_6"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_7"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_8"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_9"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_10"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_11"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_12"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_13"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_14"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  fitTau["wave_15"] = { {2.64,0}, {8.84,0}, {20.12,0} };
  dataProcessor.setFitTau(fitTau);

  int status = 0;
  for (int ev = 0; ev < nEventsToProcess; ev++) {
    if (status != 0) break;
    if (ev % 1000 == 0)
      cout << "Event: " << ev << endl;
    status = rawDataHandler.skip(header_length);
    status = rawDataHandler.readEvent(wfm_length);
    if (visualize && dataVisualizer.getStatus()) dataVisualizer.fillGraphs(rawDataHandler.dataFiles);
    bool good_fit = dataProcessor.process(rawDataHandler.dataFiles);
    if (visualize && dataVisualizer.getStatus()) {
      if (good_fit) dataVisualizer.fillFit(dataProcessor.getOutDigi());
      status = dataVisualizer.visualize();
    }
  }

  return 0;
}
