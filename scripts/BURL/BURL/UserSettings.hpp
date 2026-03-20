#ifndef UserSettings_hpp
#define UserSettings_hpp

#include <chrono>
#include <string>


class UserSettings {

    public:
        static UserSettings&        userSettings(void)
                                        {
                                        static UserSettings us;
                                        return us;
                                        }
        void                        appendOutputFile(std::string s) { outputFile = outputFile + s;}
        void                        initializeSettings(int argc, const char* argv[]);
        unsigned long               getChainLength(void) { checkSettings(); return chainLength; }
        std::string                 getInputFile(void) { checkSettings(); return inputFile; }
        std::string                 getInputTree(void) { checkSettings(); return inputTree; }
        bool                        getLogTransformData(void) { checkSettings(); return logTransformData; }
        int                         getNumChains(void) { return numChains; }
        int                         getNumThreads(void) { return numThreads; }
        std::string                 getOutputFile(void) { checkSettings(); return outputFile; }
        int                         getPrintFrequency(void) { checkSettings(); return printFrequency; }
        std::string                 getReadDataType(void) { checkSettings(); return readDatDatatype; }
        int                         getSampleFrequency(void) { checkSettings(); return sampleFrequency; }
        void                        print(void);
        void                        printHelp(void);
        void                        startTiming(void);
        void                        endTiming(void);
        void                        writeLog(void);

    private:
                                    UserSettings(void) { settingsInitialized = false; }
                                    UserSettings(const UserSettings& u);
        UserSettings&               operator=(const UserSettings&);
        void                        checkSettings(void);
        std::chrono::steady_clock::time_point startTime;
        std::string                 executablePath;
        std::string                 inputFile;
        std::string                 inputTree;
        std::string                 logFile;
        std::string                 outputFile;
        std::string                 readDatDatatype;
        unsigned long               chainLength;
        int                         numChains;
        int                         numThreads;
        int                         printFrequency;
        int                         sampleFrequency;
        bool                        settingsInitialized;
        bool                        logTransformData;
};

#endif
