#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm> // for std::copy
#include <string>

void badInput(std::string received, std::string desired){
  std::cout << "Input file gave me " << received << " but i wanted " << desired << std::endl;
}

template <typename InputFile, typename String, typename Range>
bool readIntoVector(InputFile& inputFile, String dataName, Range vec){
  std::string inStr; inputFile >> inStr;
  if (inStr != dataName){badInput(inStr,dataName); return false;}
  while(inputFile>>inStr){
    if (inStr == "----"){ break; }
    vec.push_back(std::stod(inStr)); 
  }
  return true;
}


int main() {

  std::ifstream inputFile("numbers.txt");
  if(!inputFile) { std::cout << "Cannot open input inputFile.\n"; return 1; }

  std::string inStr;
  std::vector<double> rhoGrid, rho, alphas, betas;
  double T;

  if (not readIntoVector(inputFile, "rhoGrid", rhoGrid)){ return 1; }
  if (not readIntoVector(inputFile, "rhoValues", rho)){ return 1; }
  if (not readIntoVector(inputFile, "alpha", alphas)){ return 1; }
  if (not readIntoVector(inputFile, "beta", betas)){ return 1; }

  // Read in other parameters -------------------------------------------------
  inputFile >> inStr; if (inStr != "temperature"){badInput(inStr,"temperature"); return 1;}
  inputFile >> T;

  inputFile.close();
  return 0;

}
