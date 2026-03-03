#include "Eigen/Dense"
#include "Msg.hpp"
#include "WriteTSV.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>

WriteTSV::WriteTSV(void){
}

WriteTSV::WriteTSV(std::string filepath, bool overwrite) : filepath(filepath), numRows(0), numCols(0), rnFlag(false){
    if(overwrite == true)
        std::remove(filepath.c_str());
    fout.open(filepath, std::ios::out | std::ios::app);
}

WriteTSV::~WriteTSV(void){
    fout.close();
}

void WriteTSV::addColumnNamesTSV(std::vector<std::string> cn){
    if(numRows == 0){
        for(int i = 0; i < cn.size(); i++){
            if(i < (cn.size() - 1))
                fout << cn[i] << "\t";
            else
                fout << cn[i];
        }
        fout << "\n";
        numRows++;
    }else{
        Msg::error("adding column names after data has alread been entered");
    }
    fout.flush();
}

void WriteTSV::addFilepath(std::string fp, bool overwrite){
    if (fout.is_open()) {
        fout.close();
    }
    filepath = fp;
    if(overwrite == true)
        std::remove(filepath.c_str());
    fout.open(filepath, std::ios::out | std::ios::app);
    numRows = 0;
    numCols = 0;
    rnFlag = false;
    if (!fout.is_open())
        Msg::error("File stream is not open: " + filepath);
}

void WriteTSV::addRownamesTSV(std::vector<std::string> rn){
    rownames = rn;
    rnFlag = true;
}

void WriteTSV::appendDataTSV(Eigen::MatrixXd x){
    if (!fout.is_open())
        Msg::error("Attempting to write to TSV file before opening it.");

    std::vector<double> row;
    row.clear();
    for(int i = 0; i < x.rows(); i++){
        for(int j = 0; j < x.cols(); j++)
            row.push_back(x(i, j));
        appendDataTSV(row);
        row.clear();
    }
}

void WriteTSV::appendDataTSV(std::vector<double> data){
    if (!fout.is_open())
        Msg::error("Attempting to write to TSV file before opening it.");

    if(rnFlag == false){
        for(int i = 0; i < data.size(); i++){
        if(i < (data.size()-1))
            fout << data[i] << "\t";
        else
            fout << data[i];
        }
        fout << "\n";
        numRows++;
    }else{
        if (numRows >= rownames.size())
            Msg::error("Not enough row names provided");
        fout << rownames[numRows] << "\t";
        for(int i = 0; i < data.size(); i++){
        if(i < (data.size()-1))
            fout << data[i] << "\t";
        else
            fout << data[i];
        }
        fout << "\n";
        numRows++;
    }
    fout.flush();
}

void WriteTSV::appendDataTSV(std::vector<std::string> data){
    if (!fout.is_open())
        Msg::error("Attempting to write to TSV file before opening it.");
    if(rnFlag == false){
        for(int i = 0; i < data.size(); i++){
        if(i < (data.size()-1))
            fout << data[i] << "\t";
        else
            fout << data[i];
        }
        fout << "\n";
        numRows++;
    }else{
        if (numRows >= rownames.size())
            Msg::error("Not enough row names provided");
        fout << rownames[numRows] << "\t";
        for(int i = 0; i < data.size(); i++){
        if(i < (data.size()-1))
            fout << data[i] << "\t";
        else
            fout << data[i];
        }
        fout << "\n";
        numRows++;
    }
    fout.flush();
}

void WriteTSV::appendDataTSV(std::string data){
    fout << data;
    fout << "\n";
    numRows++;
    fout.flush();
}

void WriteTSV::closeTSV(void){
    fout.close();
}
