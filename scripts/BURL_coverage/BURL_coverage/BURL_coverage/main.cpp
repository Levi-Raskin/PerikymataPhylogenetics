#include "Mcmc.hpp"
#include "Msg.hpp"
#include "MetropolisCoupledMcmc.hpp"
#include "PerikymataHSPv4.hpp"
#include "ReadTSV.hpp"
#include "ReadCSV.hpp"
#include "SimulateData.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

#include <iostream>
#include <chrono>
#include <iomanip>
#include <sstream>

void printHeader(void);
std::string formatDuration(double seconds);

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
    
    SimulateData s = SimulateData();
    
    int numReps = settings.getNumReps();
    auto wallStart = std::chrono::steady_clock::now();

    for(int i = 0; i < numReps; i++){
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::cout << "On rep " << i + 1 << " of " << numReps << "\n";

        if (i > 0) {
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - wallStart).count();
            double avgPerRep = elapsed / i;
            double eta = avgPerRep * (numReps - i);
            std::cout << "Elapsed: " << formatDuration(elapsed)
                      << "  |  Avg/rep: " << formatDuration(avgPerRep)
                      << "  |  ETA: " << formatDuration(eta) << "\n";
        }

        std::cout << "-----------------------------------------------------------------------" << std::endl;
        
        // Simulate data
        s.simulateData();
        
        if(numChains > 1){
            std::vector<PhylogeneticModel*> perikymataModels;
            perikymataModels.resize(numChains);
            for(int i = 0; i < numChains; i++)
                perikymataModels[i] = new PerikymataHSPv4(s.getSimulatedTree(), s.getSimulatedRownames(), s.getSimulatedData());

            MetropolisCoupledMcmc mcmc(numCycles, pf, sf, perikymataModels);
            mcmc.run();
        } else if (numChains == 1){
            PhylogeneticModel* perikymataModel = new PerikymataHSPv4(s.getSimulatedTree(), s.getSimulatedRownames(), s.getSimulatedData());
            Mcmc mcmc(numCycles, pf, sf, perikymataModel);
            mcmc.run();
        }
        
        // Check coverage
        s.checkCredInt();
        s.print();
    }

    {
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - wallStart).count();
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::cout << "All " << numReps << " reps completed in " << formatDuration(elapsed) << std::endl;
        std::cout << "-----------------------------------------------------------------------" << std::endl;
    }
    
    s.writeCoverage();
    settings.endTiming();
    return 0;
}

std::string formatDuration(double seconds) {
    std::ostringstream oss;
    if (seconds < 60.0) {
        oss << std::fixed << std::setprecision(1) << seconds << "s";
    } else if (seconds < 3600.0) {
        int m = static_cast<int>(seconds) / 60;
        int s = static_cast<int>(seconds) % 60;
        oss << m << "m " << s << "s";
    } else {
        int h = static_cast<int>(seconds) / 3600;
        int m = (static_cast<int>(seconds) % 3600) / 60;
        int s = static_cast<int>(seconds) % 60;
        oss << h << "h " << m << "m " << s << "s";
    }
    return oss.str();
}

void printHeader(void) {
    std::cout << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "   BURLc - validation of BURL through coverage checks" << std::endl;
    std::cout << "   * Levi Yoder Raskin (University of California, Berkeley)" << std::endl;
    std::cout << "   * John P. Huelsenbeck (University of California, Berkeley)" << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
}
