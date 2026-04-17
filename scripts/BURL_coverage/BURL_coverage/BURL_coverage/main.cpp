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
    settings.writeLog();
    
    
    settings.startTiming();
    
    int numChains = settings.getNumChains();
    unsigned long numCycles = settings.getChainLength();
    int pf = settings.getPrintFrequency();
    int sf = settings.getSampleFrequency();
    
    for(int i = 0; i < settings.getNumReps(); i++){
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::cout << "On rep " << i + 1 << " of " << settings.getNumReps() << "\n";
        std::cout << "Current coverage: " << 10 << "\n";
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        
        
        //simulate data
        
        
        
        //infer parameters from simulated data
        if(numChains > 1){
            std::vector<PhylogeneticModel*> perikymataModels;
            perikymataModels.resize(numChains);
            for(int i = 0; i < numChains; i++)
                perikymataModels[i] = new PerikymataHSPv4(&tree, rawReadDatNames, &readDat);

            MetropolisCoupledMcmc mcmc(numCycles, pf, sf, perikymataModels);
            mcmc.run();
        }else if (numChains == 1){
            PhylogeneticModel* perikymataModel = new PerikymataHSPv4(&tree, rawReadDatNames, &readDat);
            Mcmc mcmc(numCycles, pf, sf, perikymataModel);
            mcmc.run();
        }
    //check coverage
    
    }
    
    //overwrite trace file with coverage
    
    settings.endTiming();
    return 0;
}

void printHeader(void) {

    std::cout << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "   BURLc - validation of BURL through coverage checks" << std::endl;
    std::cout << "   * Levi Yoder Raskin (University of California, Berkeley)" << std::endl;
    std::cout << "   * John P. Huelsenbeck (University of California, Berkeley)" << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
}
