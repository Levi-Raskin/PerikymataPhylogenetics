#include <iostream>
#include <string>
#include <vector>
#include <thread>

#include "Msg.hpp"
#include "UserSettings.hpp"


void UserSettings::checkSettings(void) {

    if (settingsInitialized == false)
        Msg::error("Settings are not initialized");
}

void UserSettings::initializeSettings(int argc, const char* argv[]) {

    if (settingsInitialized == true)
        {
        Msg::warning("Settings have already been initialized");
        return;
        }

    inputFile = "";
    inputTree = "";
    outputFile = "";
    chainLength = 100;
    numChains = 4;
    numThreads = 4;
    printFrequency = 1;
    sampleFrequency = 1;

    // put the arguments into a vector of strings to make things easier
    std::vector<std::string> arguments;
    for (int i=0; i<argc; i++)
        arguments.push_back(std::string(argv[i]));
        
    executablePath = arguments[0];
    std::string arg = "";
    for (int i=1; i<arguments.size(); i++)
        {
        if (arg == "")
            {
            arg = arguments[i];
            }
        else
            {
            if (arg == "-i"){
                inputFile = arguments[i];
                readDatDatatype = inputFile.substr(inputFile.length() - 3);
            }else if (arg == "-it")
                inputTree = arguments[i];
            else if (arg == "-o")
                outputFile = arguments[i];
            else if (arg == "-n")
                chainLength = atoi(arguments[i].c_str());
            else if (arg == "-p")
                printFrequency = atoi(arguments[i].c_str());
            else if (arg == "-s")
                sampleFrequency = atoi(arguments[i].c_str());
            else if (arg == "-c")
                numChains = atoi(arguments[i].c_str());
            else if (arg == "-nt")
                numThreads = atoi(arguments[i].c_str());
            else if (arg == "-help"){
                std::cout << "Usage: BURL [options]\n";
                std::cout << "Options:\n";
                std::cout << "  -i  <file>    Input file        (expecting TSV or CSV)\n";
                std::cout << "                  *Input file requires rownames labeling species and colnames labeling traits\n";
                std::cout << "  -it <file>    Input tree file   (expecting .txt file)\n";
                std::cout << "  -o  <file>    Output file       (Outputs TSV only)\n";
                std::cout << "  -n  <int>     Chain length\n";
                std::cout << "  -p  <int>     Print frequency\n";
                std::cout << "  -s  <int>     Sample frequency\n";
                std::cout << "  -c  <int>     Number of chains  (1 for MCMC, 2+ for Metropolis-coupled MCMC)\n";
                std::cout << "  -nt <int>     Number of threads\n";
            }
            else
                Msg::error("Unknown argument \"" + arg + "\"");
            arg = "";
            }
        }
    settingsInitialized = true;
    
    if(numChains < 1){
        std::cout << "Requested " << numChains << " chains, minimum required is one; replacing with one. \n";
        numChains = 1;
    }
    
    if(numThreads < 1){
        std::cout << "Requested " << numThreads << " threads, minimum required is one; replacing with one. \n";
        numThreads = 1;
    }
    
    if(chainLength < 10){
        std::cout << "Requested " << chainLength << " MCMC cycles; replacing with 10. Expect nonconvergence!. \n";
        chainLength = 10;
    }
    if(printFrequency < 1){
        std::cout << "Requested " << printFrequency << " console updates, minimum required is one; replacing with one. \n";
        printFrequency = 1;
    }
    if(sampleFrequency < 1){
        std::cout << "Requested " << sampleFrequency << " samples from chain, minimum required is one; replacing with one. \n";
        sampleFrequency = 1;
    }
    
    if(numThreads > numChains){
        std::cout << "Requested " << numThreads << " threads, but this program only supports at most an equal number of threads as chains.\n";
        numThreads = numChains;
        std::cout << "Adjusted number of threads to " << numThreads << ".\n";
    }
    
    int maxNumThreads = std::thread::hardware_concurrency();
    if(numThreads >= maxNumThreads){
        std::cout << "Requested " << numThreads << " threads but you only have available " << maxNumThreads << "\n";
        numThreads = maxNumThreads - 1;
        std::cout << "Automatically replacing with " << numThreads << "\n";
    }
    
}

void UserSettings::print(void) {

    checkSettings();
    std::cout << "Input file name:                       " << inputFile << std::endl;
    std::cout << "Input tree file name:                  " << inputTree << std::endl;
    std::cout << "Output file name:                      " << outputFile << std::endl;
    std::cout << "Chain length:                          " << chainLength << std::endl;
    std::cout << "Number of chains:                      " << numChains << std::endl;
    std::cout << "Number of threads:                     " << numThreads << std::endl;
    std::cout << "Print-to-screen frequency:             " << printFrequency << std::endl;
    std::cout << "Chain sampling frequency:              " << sampleFrequency << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;

}
