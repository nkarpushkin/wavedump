#ifndef DIGIDATA_H
#define DIGIDATA_H

#include <vector>
#include <complex>
#include "BmnDigiContainerTemplate.cxx"

class DigiData {
public:
    enum class SignalType { AnalogPositive, AnalogNegative, LogicPositive, LogicNegative, Blank, Trigger };

    DigiData();
    DigiData(const DigiData& other); // Copy constructor
    DigiData(DigiData&& other) noexcept; // Move constructor
    DigiData& operator=(const DigiData& other); // Copy assignment operator
    DigiData& operator=(DigiData&& other) noexcept; // Move assignment operator
    virtual ~DigiData();

    SignalType type;
    BmnDigiContainerTemplate container;
    bool fitflag;
    std::vector<std::pair<float,float>> tau;
};

DigiData::DigiData() {
    // Constructor implementation
}

DigiData::DigiData(const DigiData& other) {
    // Copy constructor implementation
    type = other.type;
    container = other.container;
    fitflag = other.fitflag;
    tau = other.tau;
}

DigiData::DigiData(DigiData&& other) noexcept {
    // Move constructor implementation
    type = std::move(other.type);
    container = std::move(other.container);
    tau = std::move(other.tau);
}

DigiData& DigiData::operator=(const DigiData& other) {
    // Copy assignment operator implementation
    if (this != &other) {
        type = other.type;
        container = other.container;
        fitflag = other.fitflag;
        tau = other.tau;
    }
    return *this;
}

DigiData& DigiData::operator=(DigiData&& other) noexcept {
    // Move assignment operator implementation
    if (this != &other) {
        type = std::move(other.type);
        container = std::move(other.container);
        fitflag = std::move(other.fitflag);
        tau = std::move(other.tau);
    }
    return *this;
}

DigiData::~DigiData() {
    // Destructor implementation
}

#endif /* DIGIDATA_H */
