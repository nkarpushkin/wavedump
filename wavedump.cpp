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
  const int visual_min_adc = -5000;
  const int visual_max_adc = 5000;
  const int gate_beg = 150;
  const int gate_end = 600;

  RawDataHandler rawDataHandler;
  DataVisualizer dataVisualizer(wfm_length, std::make_pair(visual_min_adc, visual_max_adc));

  auto file_list = getFilesInFolder(folderPath);
  for (auto& filename : file_list)
    rawDataHandler.openFile(folderPath + "/" + filename);

  //auto size = std::filesystem::file_size((file_list + file_list.at(0).c_str()));
  //long nEvents = floor(size / wfm_length);
  //nEventsToProcess = (nEventsToProcess <= 0) ? nEvents : std::min(nEvents, nEventsToProcess);

  DataProcessor dataProcessor(rawDataHandler.dataFiles);
  dataProcessor.setCommonGates(std::make_pair(gate_beg, gate_end));
  dataProcessor.setOutPath(outPath);

  std::map<std::string, SignalType> sigTypes;
  sigTypes["TR_0_0"] = SignalType::LogicNegative;
  sigTypes["wave_0"] = SignalType::AnalogPositive;
  sigTypes["wave_1"] = SignalType::AnalogPositive;
  sigTypes["wave_2"] = SignalType::Blank;
  sigTypes["wave_3"] = SignalType::LogicNegative;
  sigTypes["wave_4"] = SignalType::AnalogPositive;
  sigTypes["wave_5"] = SignalType::LogicNegative;
  sigTypes["wave_6"] = SignalType::AnalogPositive;
  sigTypes["wave_7"] = SignalType::Blank;
  dataProcessor.setSignalTypes(sigTypes);

  std::map<std::string, bool> fitFlags;
  fitFlags["wave_0"] = true;
  fitFlags["wave_1"] = true;
  fitFlags["wave_4"] = true;
  fitFlags["wave_6"] = true;
  dataProcessor.setFitFlags(fitFlags);

  int status = 0;
  for (int ev = 0; ev < nEventsToProcess; ev++) {
    if (status != 0) break;
    if (ev % 1000 == 0)
      cout << "Event: " << ev << endl;
    //status = rawDataHandler.skip(header_length);
    status = rawDataHandler.readEvent(wfm_length);
    if (visualize && dataVisualizer.getStatus()) dataVisualizer.fillGraphs(rawDataHandler.dataFiles);
    bool good_fit = dataProcessor.process(rawDataHandler.dataFiles);
    if (visualize && dataVisualizer.getStatus()) {
      if (good_fit) dataVisualizer.fillFit(dataProcessor.getOutDigi(), fitFlags);
      status = dataVisualizer.visualize();
    }
  }

  return 0;
}
