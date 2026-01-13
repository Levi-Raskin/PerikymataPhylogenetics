#include "Msg.hpp"
#include "WriteCSV.hpp"

#include <cstdio>
#include <fstream>

WriteCSV::WriteCSV(std::string filepath, bool overwrite) : filepath(filepath), numRows(0), numCols(0){
    if(overwrite == true)
        std::remove(filepath.c_str());
    fout.open(filepath, std::ios::out | std::ios::app);
}

WriteCSV::~WriteCSV(void){
    fout.close();
}

void WriteCSV::addColumnNamesCSV(std::vector<std::string> cn){
    if(numRows == 0){
        for(int i = 0; i < cn.size(); i++){
            if(i < (cn.size() - 1))
                fout << cn[i] << ", ";
            else
                fout << cn[i];
        }
        fout << "\n";
        numRows++;
    }else{
        Msg::error("adding column names after data has alread been entered");
    }
}

void WriteCSV::appendDataCSV(std::vector<double> data){
    for(int i = 0; i < data.size(); i++){
        if(i < (data.size()-1))
            fout << data[i] << ", ";
        else
            fout << data[i];
    }
    fout << "\n";
    numRows++;
}

void WriteCSV::appendDataCSV(double data){
    fout << data;
    fout << "\n";
    numRows++;
}

void WriteCSV::appendDataCSV(std::string data){
    fout << data;
    fout << "\n";
    numRows++;
}


void WriteCSV::closeCSV(void){
    fout.close();
}
