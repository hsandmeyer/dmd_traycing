#ifndef MEASURE_TIME_H
#define MEASURE_TIME_H
#include <chrono>
#include <functional>
#include <iostream>
#include <ostream>
#include <typeinfo>

using Seconds = std::chrono::duration<double>;

template <typename Function, typename OStream, typename... Args>
auto measure_time(OStream &out, Function &&toTime, Args &&... a)
{
    auto start{std::chrono::steady_clock::now()};
    std::invoke(std::forward<Function>(toTime), std::forward<Args>(a)...);
    auto stop{std::chrono::steady_clock::now()};
    out << "Time used for " << typeid(toTime).name() << ": "
        << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << std::endl;

    return stop - start;
}

#endif
