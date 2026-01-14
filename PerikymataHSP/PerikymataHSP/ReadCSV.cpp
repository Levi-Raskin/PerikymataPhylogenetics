#include "ReadCSV.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

ReadCSV::ReadCSV(std::string filepath, bool rn, bool cn) : readingString(false){
    
    //file path is full string leading to where the CSV file is
    //data class is string indicating expected data type-- double, int supported right now
    //missingDatSym is what is the expected symbol for missing data in this csv
    //bool rn indicates whether we are expecting rownames or not
    //bool cn indicates whether we are expecting colnames or not
    
    std::ifstream file(filepath);
    
    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
    }
    
    doubleData.clear();
    
    std::string line;
    int rowCounter = 0;
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string cell;
        int colCounter = 0;
        while (std::getline(ss, cell, ',')) {
            if((rn == true && colCounter == 0) && ((cn == false) || (cn == true && rowCounter > 0)))
                rownames.push_back(cell);
            else if(cn == true && rowCounter == 0)
                colnames.push_back(cell);
            else {
                try {
                    double cellDat = stod(cell);  // Convert string to double
                    row.push_back(cellDat);       // Store the result
                } catch (const std::invalid_argument&) {
                    // Handle non-numeric string
                    row.push_back(std::numeric_limits<double>::quiet_NaN());
                } catch (const std::out_of_range&) {
                    // Handle value out of range for double
                    row.push_back(std::numeric_limits<double>::quiet_NaN());
                }
            }
            
                
            colCounter++;
        }
        rowCounter++;
        if(row.size() > 0)
            doubleData.push_back(row);
    }
    file.close();
}

ReadCSV::ReadCSV(std::string filepath, bool rn, bool cn, bool string) : readingString(true) {
    
    //file path is full string leading to where the CSV file is
    //data class is string indicating expected data type-- double, int supported right now
    //missingDatSym is what is the expected symbol for missing data in this csv
    //bool rn indicates whether we are expecting rownames or not
    //bool cn indicates whether we are expecting colnames or not
    
    std::ifstream file(filepath);
    
    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
    }
    
    stringData.clear();
    
    std::string line;
    int rowCounter = 0;
    while (std::getline(file, line)) {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;
        int colCounter = 0;
        while (std::getline(ss, cell, ',')) {
            if((rn == true && colCounter == 0) && ((cn == false) || (cn == true && rowCounter > 0)))
                rownames.push_back(cell);
            else if(cn == true && rowCounter == 0)
                colnames.push_back(cell);
            else {
                row.push_back(cell);
            }
            
                
            colCounter++;
        }
        rowCounter++;
        if(row.size() > 0)
            stringData.push_back(row);
    }
    file.close();
}

Eigen::MatrixXd ReadCSV::getEigenMat(void){
    scratch = Eigen::MatrixXd::Zero(doubleData.size(), doubleData[0].size());
    for(int i = 0; i < doubleData.size(); i++){
        for(int j = 0; j < doubleData[0].size(); j++){
            scratch(i, j) = doubleData[i][j];
        }
    }
    return scratch;
}

void ReadCSV::print(void){
    if(readingString == false){
        for(int i = 0; i < doubleData.size(); i++){
            for(int j = 0; j < doubleData[0].size(); j++){
                std::cout << doubleData[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }else{
        for(int i = 0; i < stringData.size(); i++){
            for(int j = 0; j < stringData[0].size(); j++){
                std::cout << stringData[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
    
        
}
