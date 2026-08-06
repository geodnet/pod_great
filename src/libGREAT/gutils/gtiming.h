#ifndef GTIMING_H
#define GTIMING_H

#include <string>

#include <chrono>

namespace great
{

template <typename F>
inline void timing_part(const std::string& info,double& time,F&& func)
{

    std::chrono::high_resolution_clock::time_point end_t, beg_t;
    beg_t = std::chrono::high_resolution_clock::now();

    func();

    end_t = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - beg_t).count() / 1000.0;
    //printf("[INFO] %s Time: %lf sec\n", info.c_str(), elapsed);

    time+= elapsed;
}

template <typename F>
inline void timing(const std::string& info,F&& func)
{
    double result = 0.0;
    timing_part(info,result,std::forward<F>(func));
}


}
#endif // !GTIMING_H