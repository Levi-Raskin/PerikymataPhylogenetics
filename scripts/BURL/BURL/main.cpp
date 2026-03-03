#include "Mcmc.hpp"
#include "Msg.hpp"
#include "MetropolisCoupledMcmc.hpp"
#include "PerikymataHSPv4.hpp"
#include "ReadTSV.hpp"
#include "ReadCSV.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

#include <iostream>

void printHeader(void);

int main(int argc, const char* argv[]) {
    printHeader();
    UserSettings& settings = UserSettings::userSettings();
    settings.initializeSettings(argc, argv);
    settings.print();
    std::string readDatType = settings.getReadDataType();

    //read in raw data
    std::vector<std::string> rawReadDatNames;
    Eigen::MatrixXd readDat;
    if(readDatType == "csv" || readDatType == "CSV"){
        ReadCSV r = ReadCSV(settings.getInputFile(), true, true);
        rawReadDatNames = r.getRownames();
        readDat = r.getEigenMat();
    }else if (readDatType == "tsv"|| readDatType == "TSV"){
        ReadTSV r = ReadTSV(settings.getInputFile(), true, true);
        rawReadDatNames = r.getRownames();
        readDat = r.getEigenMat();
    }else{
        Msg::error("Data input file type is " + readDatType + " expecting tsv or csv");
    }
    
    //remove quotation marks if we have them
    for(int i = 0 ; i < rawReadDatNames.size(); i++){
        std::string trimmed_name = rawReadDatNames[i];
        trimmed_name.erase(0, trimmed_name.find_first_not_of('"'));
        trimmed_name.erase(trimmed_name.find_last_not_of('"') + 1);
        rawReadDatNames[i] =trimmed_name;
    }
    
    //read in tree file
    std::ifstream file(settings.getInputTree());
    if (!file.is_open())
        Msg::error("Error opening tree file: " + settings.getInputFile());
    std::string treeNewick;
    std::getline(file, treeNewick);
    file.close();
    Tree tree = Tree(treeNewick);
    
    int numChains = settings.getNumChains();
    unsigned long numCycles = settings.getChainLength();
    int pf = settings.getPrintFrequency();
    int sf = settings.getSampleFrequency();
    if(numChains > 1){
        std::cout << "Running Metropolis-coupled MCMC with " << numChains << " chains parallelized across " << settings.getNumThreads() << " threads \n";
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::vector<PhylogeneticModel*> perikymataModels;
        perikymataModels.resize(numChains);
        for(int i = 0; i < numChains; i++)
            perikymataModels[i] = new PerikymataHSPv4(&tree, rawReadDatNames, &readDat);

        MetropolisCoupledMcmc mcmc(numCycles, pf, sf, perikymataModels);
        mcmc.run();
    }else if (numChains == 1){
        std::cout << "Running standard MCMC \n";
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        PhylogeneticModel* perikymataModel = new PerikymataHSPv4(&tree, rawReadDatNames, &readDat);
        Mcmc mcmc(numCycles, pf, sf, perikymataModel);
        mcmc.run();
    }
    
    return 0;
}

void printHeader(void) {

    std::cout << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "   BURL - Between and within-group Uncertainty in Rates across Lineages" << std::endl;
    std::cout << "   * Levi Yoder Raskin (University of California, Berkeley)" << std::endl;
    std::cout << "   * John P. Huelsenbeck (University of California, Berkeley)" << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
}
