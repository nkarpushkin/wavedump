#ifndef RAWDATAHANDLER_H
#define RAWDATAHANDLER_H

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <regex>

class RawDataHandler {
public:
  RawDataHandler();
  ~RawDataHandler();

  struct FileData {
    int board;
    unsigned int channel;
    std::string name;
    std::ifstream fileStream;
    std::vector<float> data;
  };
  std::vector<FileData> dataFiles;

  string getFileNameWithoutExtension(const std::string& filePath);
  int openFile(const std::string& path);
  int readEvent(int length);
  int skip(int length);
  void closeFiles();
};

string RawDataHandler::getFileNameWithoutExtension(const std::string& filePath) {
  // Find the position of the last directory separator
  size_t lastSlashPos = filePath.find_last_of("/");

  // If no directory separator is found, use beginning of the string
  if (lastSlashPos == std::string::npos) {
    lastSlashPos = 0;
  }
  else {
    // Include the character after the last directory separator
    lastSlashPos++;
  }

  // Find the position of the last dot (.)
  size_t dotPos = filePath.rfind('.');

  // If a dot is found and it's after the last directory separator
  if (dotPos != std::string::npos && dotPos > lastSlashPos) {
    // Extract the substring between the last directory separator and the last dot
    return filePath.substr(lastSlashPos, dotPos - lastSlashPos);
  }
  else {
    // If no dot is found or it's before the last directory separator, return the whole substring after the last directory separator
    return filePath.substr(lastSlashPos);
  }
}

// Definition of the static dataFiles member variable
 //std::vector<RawDataHandler::FileData> RawDataHandler::dataFiles;
int RawDataHandler::openFile(const std::string& path) {
  FileData newFile;
  newFile.fileStream.open(path, std::ios::binary);
  newFile.name = getFileNameWithoutExtension(path);

  if (!newFile.fileStream.is_open()) {
    std::cerr << "Failed to open file: " << path << std::endl;
    // Handle error, possibly return or throw an exception
    return 1;
  }

  // Move the file pointer to byte with channel. TODO -- not working with channels > 7
  /*
  newFile.fileStream.seekg(12, std::ios::beg);
  if (newFile.fileStream.read(reinterpret_cast<char*>(&newFile.channel), sizeof(newFile.channel))) {
      std::cout << "File " << path << " Channel : " << newFile.channel << std::endl;
  } else {
      std::cerr << "Failed to read byte." << std::endl;
  }
  newFile.fileStream.seekg(0, std::ios::beg);
  */
  std::regex pattern("_([0-9]+)\\.dat");
  std::smatch match;
  if (std::regex_search(path, match, pattern)) {
    if (match.size() > 1) {
      // Extract the matched integer from the regex match
      std::string numberString = match[1].str();

      // Convert the substring to an integer
      try {
        int number = std::stoi(numberString);
        newFile.channel = number;
      }
      catch (const std::exception& e) {
        std::cerr << "Failed to convert string to integer: " << e.what() << std::endl;
        return 1;
      }
    }
    else {
      std::cerr << "Integer not found in the file path." << std::endl;
      return 1;
    }
  }
  else {
    std::cerr << "Regex pattern not matched in the file path." << std::endl;
    return 1;
  }
  if (path.find("TR_0_0") != std::string::npos) newFile.channel = 16;
  if (path.find("TR_0_1") != std::string::npos) newFile.channel = 17;
  dataFiles.push_back(std::move(newFile)); // Move newFile into dataFiles

  return 0;
}

int RawDataHandler::skip(int length) {
  for (auto& fileData : dataFiles) {
    fileData.data.clear();
    fileData.data.resize(length);
    for (int j = 0; j < length; j++) {
      int buf;
      fileData.fileStream.read(reinterpret_cast<char*>(&buf), sizeof(buf));
      if (fileData.fileStream.eof()) return 1;
    }
  }
  return 0;
}

int RawDataHandler::readEvent(int length) {
  for (auto& fileData : dataFiles) {
    fileData.data.clear();
    fileData.data.resize(length);
    for (int j = 0; j < length; j++) {
      float buf;
      fileData.fileStream.read(reinterpret_cast<char*>(&buf), sizeof(buf));
      if (fileData.fileStream.eof()) return 1;
      fileData.data.at(j) = buf;
    }
  }
  return 0;
}

void RawDataHandler::closeFiles() {
  // Close files after reading
  for (auto& fileData : dataFiles)
    fileData.fileStream.close();
}

RawDataHandler::RawDataHandler() {}
RawDataHandler::~RawDataHandler() {
  closeFiles();
}

#endif // RAWDATAHANDLER_H
