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
    inputFile       = "";
    inputTree       = "";
    outputFile      = "";
    chainLength     = 100;
    numChains       = 4;
    numThreads      = 4;
    printFrequency  = 1;
    sampleFrequency = 1;
    logTransformData = false;

    std::vector<std::string> arguments;
    for (int i = 0; i < argc; i++)
        arguments.push_back(std::string(argv[i]));

    executablePath = arguments[0];

    // Known flags and whether they take a value
    std::set<std::string> knownFlags = {
        "-i", "-it", "-o", "-n", "-p", "-s", "-c", "-nt", "-log", "-help", "-h"
    };
    std::set<std::string> valueFlags = {
        "-i", "-it", "-o", "-n", "-p", "-s", "-c", "-nt", "-log"
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

            if (arg == "-i") {
                inputFile = val;
                std::string ext = (val.length() >= 3) ? val.substr(val.length() - 3) : "";
                if (ext != "tsv" && ext != "csv")
                    Msg::warning("Input file \"" + val + "\" does not have a .tsv or .csv extension.");
                readDatDatatype = ext;

            } else if (arg == "-it") {
                inputTree = val;
                std::string ext = (val.length() >= 3) ? val.substr(val.length() - 3) : "";
                if (ext != "txt")
                    Msg::warning("Input tree file \"" + val + "\" does not have a .txt extension.");

            } else if (arg == "-o") {
                outputFile = val;

            }else if (arg == "-log") {
                if (val == "T" || val == "true")
                    logTransformData = true;
                else if (val == "F" || val == "false")
                    logTransformData = false;
                else
                    Msg::error("Flag \"-log\" expects T/true or F/false, but got \"" + val + "\".");

            } else {
                // Integer-valued flags
                // Check all characters are digits (allowing leading minus for negative detection)
                bool isNegative = (val[0] == '-');
                std::string digits = isNegative ? val.substr(1) : val;
                bool isInt = !digits.empty() && std::all_of(digits.begin(), digits.end(), ::isdigit);

                if (!isInt)
                    Msg::error("Flag \"" + arg + "\" expects an integer, but got \"" + val + "\".");

                int intVal = std::stoi(val);

                if (arg == "-n")        chainLength     = intVal;
                else if (arg == "-p")   printFrequency  = intVal;
                else if (arg == "-s")   sampleFrequency = intVal;
                else if (arg == "-c")   numChains       = intVal;
                else if (arg == "-nt")  numThreads      = intVal;
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

    if (inputFile.empty())
        Msg::warning("No input file specified (-i). Use -help for usage.");
    if (outputFile.empty())
        Msg::warning("No output file specified (-o). Use -help for usage.");

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

void UserSettings::printHelp(void){
    std::cout << "Usage: BURL [options]\n";
    std::cout << "Options:\n";
    std::cout << "  -i  <file>    Input file        (expecting TSV or CSV)\n";
    std::cout << "                  *Input file requires rownames labeling species and colnames labeling traits\n";
    std::cout << "  -it <file>    Input tree file   (expecting one line text file .txt, .nwk tested)\n";
    std::cout << "  -o  <file>    Output file       (Outputs TSV only)\n";
    std::cout << "  -n  <int>     Chain length\n";
    std::cout << "  -p  <int>     Print frequency\n";
    std::cout << "  -s  <int>     Sample frequency\n";
    std::cout << "  -c  <int>     Number of chains  (1 for MCMC, 2+ for Metropolis-coupled MCMC)\n";
    std::cout << "  -nt <int>     Number of threads\n";
    std::cout << "  -log <T/F>    Log-transform data (T/true or F/false)\n";
}

void UserSettings::writeLog(void){
    const std::string tsv = ".tsv";
    const std::string txt = "Log.txt";
    logFile = outputFile;
    logFile.replace(outputFile.size() - tsv.size(), tsv.size(), txt);
    
    std::ofstream log(logFile);
    if (!log.is_open())
        Msg::error("Could not open log file: " + logFile);

    log << "Input file name:                       " << inputFile << "\n";
    log << "Input tree file name:                  " << inputTree << "\n";
    log << "Output file name:                      " << outputFile << "\n";
    log << "Chain length:                          " << chainLength << "\n";
    log << "Number of chains:                      " << numChains << "\n";
    log << "Number of threads:                     " << numThreads << "\n";
    log << "Print-to-screen frequency:             " << printFrequency << "\n";
    log << "Chain sampling frequency:              " << sampleFrequency << "\n";
    log << "Log-transform data:                    " << (logTransformData ? "true" : "false") << "\n";

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
