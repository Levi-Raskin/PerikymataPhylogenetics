#ifndef ThreadPool_hpp
#define ThreadPool_hpp

#include <condition_variable>
#include <functional>
#include <iostream>
#include <mutex>
#include <queue>
#include <thread>

class ThreadPool {
    public:
                                        ThreadPool(int numThreads);
                                       ~ThreadPool();
    void                                enqueue(std::function<void()> task);
    private:
    std::vector<std::thread>            threads;
    std::queue<std::function<void()>>   tasks;
    std::mutex                          queue_mutex;
    std::condition_variable             cv;
    bool                                stop;
};

#endif
