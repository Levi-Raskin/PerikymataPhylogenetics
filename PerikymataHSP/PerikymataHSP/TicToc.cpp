#include "TicToc.hpp"
#include <chrono>
#include <iostream>

TicToc::TicToc(std::string h) : header(h) {

    startTime = std::chrono::high_resolution_clock::now();
}

void TicToc::end(void){

    std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
    
    std::chrono::nanoseconds duration = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime);
    if(duration.count()  < 9999 ){
        std::cout << header << ": " << duration.count() << " nanoseconds" <<std::endl;
    }else if (duration.count()  < 9999999){
        std::cout << header << ": " << 0.001*duration.count() << " microseconds" <<std::endl;
    }else if (duration.count()  < 9999999999){
        std::cout << header << ": " << 0.000001*duration.count() << " milliseconds" <<std::endl;
    }else{
        std::cout << header << ": " << 0.000000001*duration.count() << " seconds" <<std::endl;
    }
}
