// profiler.cpp
#include "MyTimer.h" // Include your profiling header

// Definitions of the FunctionStats constructor
FunctionStats::FunctionStats() : totalTimeUs(0), callCount(0) {}

// Definitions of global variables
// These are defined WITHOUT 'extern' here, meaning memory is allocated here.
std::map<std::string, FunctionStats> g_functionProfilingData;
long long g_performanceFrequency = 0;

// Definition of the InitializePerformanceFrequency function
void InitializePerformanceFrequency() {
    if (g_performanceFrequency == 0) {
        LARGE_INTEGER li;
        if (!QueryPerformanceFrequency(&li)) {
            // Error handling, though unlikely to fail on modern Windows
            g_performanceFrequency = 1; // Prevent division by zero
        }
        else {
            g_performanceFrequency = li.QuadPart;
        }
    }
}

// Definitions of FunctionTimer methods
FunctionTimer::FunctionTimer(const char* functionName) : m_functionName(functionName) {
    InitializePerformanceFrequency(); // Ensure frequency is set
    LARGE_INTEGER li;
    QueryPerformanceCounter(&li);
    m_startTime = li.QuadPart;
}

FunctionTimer::~FunctionTimer() {
    LARGE_INTEGER li;
    QueryPerformanceCounter(&li);
    long long endTime = li.QuadPart;

    long long durationCounts = endTime - m_startTime;
    long long durationUs = (durationCounts * 1000000) / g_performanceFrequency;

    // Thread-safety note: For multi-threaded applications, you'd need a mutex here
    // std::lock_guard<std::mutex> lock(g_profilingMutex); // Assuming a mutex exists
    g_functionProfilingData[m_functionName].totalTimeUs += durationUs;
    g_functionProfilingData[m_functionName].callCount++;
}

// Definition of the PrintProfilingReport function
void PrintProfilingReport() {
    std::cout << "\n--- Function Profiling Report ---\n";
    std::cout << std::fixed << std::setprecision(3); // For better output formatting

    std::cout << std::left << std::setw(60) << "Function Name"
        << std::right << std::setw(15) << "Calls"
        << std::right << std::setw(20) << "Total Time (ms)"
        << std::right << std::setw(20) << "Avg Time (ms/call)"
        << std::endl;
    std::cout << std::string(115, '-') << std::endl;

    for (const auto& pair : g_functionProfilingData) {
        std::string funcName = pair.first;
		funcName = funcName.substr(funcName.find("__cdecl") + 8);
		funcName = funcName.substr(0, funcName.find('(')); // Trim to function name only


        const FunctionStats& stats = pair.second;

        double totalTimeMs = static_cast<double>(stats.totalTimeUs) / 1000.0;
        double avgTimeMs = (stats.callCount > 0) ? (totalTimeMs / stats.callCount) : 0.0;

        std::cout << std::left << std::setw(60) << funcName
            << std::right << std::setw(15) << stats.callCount
            << std::right << std::setw(20) << totalTimeMs
            << std::right << std::setw(20) << avgTimeMs
            << std::endl;
    }
    std::cout << "---------------------------------\n";
}