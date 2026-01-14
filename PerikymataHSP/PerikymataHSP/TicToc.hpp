#ifndef TicToc_hpp
#define TicToc_hpp

#include <chrono>
#include <string>

class TicToc{
  
    public:
                TicToc(std::string h);
        void    end(void);
    private:
        std::chrono::high_resolution_clock::time_point  startTime;
        std::string                                     header;
};

#endif
