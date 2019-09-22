#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm> 
#include <string>
#include "gasket.h"

void badInput(std::string received, std::string desired){
  std::cout << "Input file gave me " << received << " but i wanted " << desired << std::endl;
}

template <typename InputFile, typename String, typename Range>
bool readIntoVector(InputFile& inputFile, String dataName, Range& vec){
  std::string inStr; inputFile >> inStr;
  if (inStr != dataName){badInput(inStr,dataName); return false;}
  while(inputFile>>inStr){
    if (inStr == "----"){ break; }
    vec.push_back(std::stod(inStr)); 
  }
  return true;
}

template <typename InputFile, typename String, typename Float>
bool readIntoValue(InputFile& inputFile, String dataName, Float& value){
  std::string inStr; inputFile >> inStr;
  if (inStr != dataName){badInput(inStr,dataName); return false;}
  inputFile >> value;
  return true;
}



int inputOutput() {

  std::ifstream inputFile("/Users/ameliajo/cats/src/gasket/numbers.txt");
  if(!inputFile) { std::cout << "Cannot open input inputFile.\n"; return 1; }

  std::string inStr;
  std::vector<double> rhoGrid, rho, alphas, betas, oscEnergies, oscWeights;
  double T, freeGasWgt, continWgt, oscWgt;

  if (not readIntoVector(inputFile, "rhoGrid", rhoGrid)){ return 1; }
  if (not readIntoVector(inputFile, "rhoValues", rho)){ return 1; }
  if (not readIntoVector(inputFile, "alpha", alphas)){ return 1; }
  if (not readIntoVector(inputFile, "beta", betas)){ return 1; }
  if (not readIntoVector(inputFile, "oscEnergies", oscEnergies)){ return 1; }
  if (not readIntoVector(inputFile, "oscWeights", oscWeights)){ return 1; }

  if (not readIntoValue(inputFile, "temperature", T)){ return 1; }
  if (not readIntoValue(inputFile, "freeWgt", freeGasWgt)){ return 1; }
  if (not readIntoValue(inputFile, "continWgt", continWgt)){ return 1; }
  if (not readIntoValue(inputFile, "oscWgt", oscWgt)){ return 1; }

  inputFile.close();

  std::vector<double> t = ranges::view::iota(0,80) 
                        | ranges::view::transform([](int i){return i*0.1;});


  std::cout << std::endl;
  std::cout << "all good!" << std::endl;
  std::cout << std::endl;

  //gasket(rhoGrid, rho, 







  return 0;

}
