#ifndef WriteCSV_hpp
#define WriteCSV_hpp

#include <fstream>
#include <string>
#include <vector>

class WriteCSV{
    public:
                        WriteCSV(std::string filepath, bool overwrite);
                       ~WriteCSV(void);
        void            addColumnNamesCSV(std::vector<std::string> cn);
        void            appendDataCSV(std::vector<double> data);
        void            appendDataCSV(double data);
        void            appendDataCSV(std::string data);
        void            closeCSV(void);
        
    private:
        std::fstream    fout;
        int             numCols;
        int             numRows;
        std::string     filepath;
};

#endif
