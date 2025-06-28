// profiler.h
#pragma once // Ensures this header is included only once per compilation unit

#include <windows.h> // For QueryPerformanceCounter, QueryPerformanceFrequency, Sleep
#include <iostream>
#include <string>
#include <map>
#include <chrono>
#include <iomanip> // For std::fixed and std::setprecision

// Define a struct to hold statistics for each function
struct FunctionStats {
    long long totalTimeUs; // Total time in microseconds
    long long callCount;

    FunctionStats(); // Declare constructor
};

// Declarations of global variables and functions (using extern for globals)
extern std::map<std::string, FunctionStats> g_functionProfilingData;
extern long long g_performanceFrequency;

// Declarations of functions
void InitializePerformanceFrequency();
void PrintProfilingReport();

// RAII class for timing functions
class FunctionTimer {
public:
    FunctionTimer(const char* functionName);
    ~FunctionTimer();

private:
    const char* m_functionName;
    long long m_startTime;
};

// Macro to easily instrument functions
// Use __FUNCSIG__ for MSVC, __PRETTY_FUNCTION__ for GCC/Clang
#ifdef _MSC_VER
#define PROFILE_FUNCTION FunctionTimer timer(__FUNCSIG__);
#else
#define PROFILE_FUNCTION FunctionTimer timer(__PRETTY_FUNCTION__);
#endif