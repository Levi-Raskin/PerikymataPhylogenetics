#include <iostream>
#include <fstream>
#include <set>
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
    if (settingsInitialized == true) {
        Msg::warning("Settings have already been initialized");
        return;
    }

    // Defaults
    outputFile      = "";
    chainLength     = 20000;
    numChains       = 10;
    numThreads      = 10;
    printFrequency  = 1000;
    sampleFrequency = 1000;
    logTransformData = false;

    std::vector<std::string> arguments;
    for (int i = 0; i < argc; i++)
        arguments.push_back(std::string(argv[i]));

    executablePath = arguments[0];

    // Known flags and whether they take a value
    std::set<std::string> knownFlags = {
        "-o", "-c", "-nreps", "-ntips", "-ntraits", "-nimp", "-nobs", "-p", "-i", "-help", "-h"
    };
    std::set<std::string> valueFlags = {
        "-o", "-c", "-nreps", "-ntips", "-ntraits", "-nimp", "-nobs", "-p", "-i"
    };

    for (int i = 1; i < (int)arguments.size(); i++) {
        std::string arg = arguments[i];

        // Check it looks like a flag
        if (arg.empty())
            Msg::error("Empty argument at position " + std::to_string(i));

        if (knownFlags.find(arg) == knownFlags.end())
            Msg::error("Unknown flag \"" + arg + "\". Use -help to see valid options.");

        // Help flag (no value)
        if (arg == "-help" || arg == "-h") {
            printHelp();
            return;
        }

        // All remaining flags require a value — check it exists
        if (valueFlags.count(arg)) {
            if (i + 1 >= (int)arguments.size())
                Msg::error("Flag \"" + arg + "\" requires a value but none was provided.");

            std::string val = arguments[++i];

            // Catch accidentally passing another flag as a value
            if (knownFlags.count(val))
                Msg::error("Flag \"" + arg + "\" expects a value, but got another flag \"" + val + "\".");

            if (arg == "-o") {
                outputFile = val;
                
            }else if (arg == "-p") {
                if (val == "T" || val == "t" || val == "true")
                    withPhylogeny = true;
                else if (val == "F" || val == "f" || val == "false")
                    withPhylogeny = false;
                else
                    Msg::error("Flag \"-p\" expects T/true or F/false, but got \"" + val + "\".");
            }else if (arg == "-i") {
                if (val == "T" || val == "t" || val == "true")
                    withIntraspecific = true;
                else if (val == "F" || val == "f" || val == "false")
                    withIntraspecific = false;
                else
                    Msg::error("Flag \"-i\" expects T/true or F/false, but got \"" + val + "\".");
            }else {
                // Integer-valued flags
                // Check all characters are digits (allowing leading minus for negative detection)
                bool isNegative = (val[0] == '-');
                std::string digits = isNegative ? val.substr(1) : val;
                bool isInt = !digits.empty() && std::all_of(digits.begin(), digits.end(), ::isdigit);

                if (!isInt)
                    Msg::error("Flag \"" + arg + "\" expects an integer, but got \"" + val + "\".");

                int intVal = std::stoi(val);
                if (arg == "-nreps")        numReps     = intVal;
                else if (arg == "-c")       chainLength = intVal;
                else if (arg == "-ntips")   numTips     = intVal;
                else if (arg == "-ntraits") numTraits   = intVal;
                else if (arg == "-nimp")    numImp      = intVal;
                else if (arg == "-nobs")    numObs      = intVal;
            }
        }
    }

    settingsInitialized = true;

    // ── Post-parse validation  ──────────────────────────────

    int maxNumThreads = (int)std::thread::hardware_concurrency();
    if (maxNumThreads <= 0) {
        Msg::warning("Could not determine hardware thread count; defaulting to 1 thread.");
        maxNumThreads = 1;
    }

    if (numChains < 1) {
        Msg::warning("Chains must be >= 1; resetting to 1.");
        numChains = 1;
    }

    if (chainLength < 10) {
        Msg::warning("Chain length " + std::to_string(chainLength) + " is very short; resetting to 10. Expect non-convergence!");
        chainLength = 10;
    }

    if (printFrequency < 1) {
        Msg::warning("Print frequency must be >= 1; resetting to 1.");
        printFrequency = 1;
    }

    if (sampleFrequency < 1) {
        Msg::warning("Sample frequency must be >= 1; resetting to 1.");
        sampleFrequency = 1;
    }

    if (numThreads >= maxNumThreads) {
        Msg::warning("Requested " + std::to_string(numThreads) +
                     " threads, but only " + std::to_string(maxNumThreads) +
                     " available; capping at " + std::to_string(maxNumThreads - 1) + ".");
        numThreads = maxNumThreads - 1;
    }

    if (numThreads > numChains) {
        Msg::warning("Threads (" + std::to_string(numThreads) +
                     ") cannot exceed chains (" + std::to_string(numChains) +
                     "); reducing threads to match.");
        numThreads = numChains;
    }

    if (numThreads < 1) {
        Msg::warning("Threads must be >= 1; resetting to 1.");
        numThreads = 1;
    }
    
    if(numReps <= 5){
        Msg::warning("Only checking coverage in five or fewer simulated data; results are meaningless.");
    }

    if (outputFile.empty())
        Msg::warning("No output file specified (-o). Use -help for usage.");

}

void UserSettings::print(void) {

    checkSettings();
    std::cout << "Output file name:                      " << outputFile << std::endl;
    std::cout << "Chain length:                          " << chainLength << std::endl;
    std::cout << "Number of chains:                      " << numChains << std::endl;
    std::cout << "Number of threads:                     " << numThreads << std::endl;
    std::cout << "Print-to-screen frequency:             " << printFrequency << std::endl;
    std::cout << "Chain sampling frequency:              " << sampleFrequency << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "Number of repititions:                 " << numReps << std::endl;
    std::cout << "Number of tips:                        " << numTips << std::endl;
    std::cout << "Number of traits:                      " << numTraits << std::endl;
    std::cout << "Number of missing observations:        " << numImp << std::endl;
    std::cout << "Number of observations per tip:        " << numObs << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
    std::cout << "With phylogeny:                        " << withPhylogeny << std::endl;
    std::cout << "With intraspecific:                    " << withIntraspecific << std::endl;
    std::cout << "-----------------------------------------------------------------------" << std::endl;
}

void UserSettings::printHelp(void){
    std::cout << "Usage: BURLc [options]\n";
    std::cout << "Options:\n";
    std::cout << "  -o          <file>      Output file without filepath ending (Outputs coverage in TSV only)\n";
    std::cout << "  -nreps      <int>       Number of repetitions for coverage check\n";
    std::cout << "  -ntips      <int>       Number of tips to simulate tree with\n";
    std::cout << "  -ntraits    <int>       Number of traits to simulate\n";
    std::cout << "  -nimp       <int>       Number of missing observations (global)\n";
    std::cout << "  -nobs       <int>       Number of observations per species\n";
    std::cout << "  -p          <T/F>       Inference with phylogeny (true/T, false/F)\n";
    std::cout << "  -i          <T/F>       Inference with intraspecific variation (true/T, false/F) \n";
}

void UserSettings::writeLog(void){
//    const std::string tsv = ".tsv";
//    const std::string txt = "Log.txt";
    logFile = outputFile + "Log.txt";
//    logFile.replace(outputFile.size() - tsv.size(), tsv.size(), txt);
    
    std::ofstream log(logFile);
    if (!log.is_open())
        Msg::error("Could not open log file: " + logFile);

    log << "Output file name:                      " << outputFile << "\n";
    log << "Chain length:                          " << chainLength << "\n";
    log << "Number of chains:                      " << numChains << "\n";
    log << "Number of threads:                     " << numThreads << "\n";
    log << "Print-to-screen frequency:             " << printFrequency << "\n";
    log << "Chain sampling frequency:              " << sampleFrequency << "\n";
    log << "---------------------------------------" << "\n";
    log << "Number of repetitions:                 " << numReps << "\n";
    log << "Number of tips:                        " << numTips << "\n";
    log << "Number of traits:                      " << numTraits << "\n";
    log << "Number of missing observations:        " << numImp << "\n";
    log << "Number of observations per tip:        " << numObs << "\n";
    log << "---------------------------------------" << "\n";


    log.close();

}

void UserSettings::startTiming(void){
    startTime = std::chrono::steady_clock::now();
}

void UserSettings::endTiming(void){
    std::ofstream log(logFile, std::ios::app);
    auto endTime = std::chrono::steady_clock::now();
    double elapsedMinutes = std::chrono::duration<double, std::ratio<60>>(endTime - startTime).count();
    log << "Time elapsed (minutes):                " << std::fixed << std::setprecision(4) << elapsedMinutes << "\n";
    log.close();
}
