// profiler.cpp
#include "MyTimer.h" // Include your profiling header
#include <algorithm>

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
	g_functionProfilingData[m_functionName].timeSamples.push_back(durationUs);
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
		<< std::right << std::setw(20) << "Min Time (ms)"
		<< std::right << std::setw(20) << "Max Time (ms)"
		<< std::right << std::setw(20) << "Std"
		<< std::endl;
	std::cout << std::string(175, '-') << std::endl;

	for (const auto& pair : g_functionProfilingData) {
		std::string funcName = pair.first;

		if (funcName.find("::") != std::string::npos) {
			funcName = funcName.substr(funcName.find("__cdecl") + 8);
		}
		else {
			if (funcName.find("class") != std::string::npos) {
				funcName.erase(funcName.find("class"), 6);
				funcName.insert(funcName.find("__cdecl"), "::");
			}
			funcName.erase(funcName.find("__cdecl"), 8);
		}

		funcName = funcName.substr(0, funcName.find('(')); // Trim to function name only
		funcName.erase(remove_if(funcName.begin(), funcName.end(), isspace), funcName.end());


		const FunctionStats& stats = pair.second;

		double totalTimeMs = static_cast<double>(stats.totalTimeUs) / 1000.0;
		double avgTimeMs = (stats.callCount > 0) ? (totalTimeMs / stats.callCount) : 0.0;
		std::vector<long long> timeSamples = stats.timeSamples;


		// calculate min, max, and median, std
		std::sort(timeSamples.begin(), timeSamples.end());
		double minTime = static_cast<double>(timeSamples.front()) / 1000.0;
		double maxTime = static_cast<double>(timeSamples.back()) / 1000.0;
		//double medianTime = (timeSamples.size() % 2 == 0)
		//	                    ? (static_cast<double>(timeSamples[timeSamples.size() / 2 - 1] + timeSamples[
		//		                    timeSamples.size() / 2]) / 2.0)
		//	                    : static_cast<double>(timeSamples[timeSamples.size() / 2]);
		////double sum = std::accumulate(timeSamples.begin(), timeSamples.end(), 0LL);
		//double meanTime = totalTimeMs / stats.callCount;
		double variance = 0.0;
		for (const auto& sample : timeSamples) {
			double sampleMs = static_cast<double>(sample) / 1000.0;
			variance += (sampleMs - avgTimeMs) * (sampleMs - avgTimeMs);
		}
		variance /= static_cast<double>(timeSamples.size());
		double stdDevTime = std::sqrt(variance);
		// print
		/*std::cout << "Time samples for " << funcName << ":\n";
		std::cout << "  Min: " << minTime << " us\n";
		std::cout << "  Max: " << maxTime << " us\n";
		std::cout << "  Median: " << medianTime << " us\n";
		std::cout << "  Mean: " << meanTime << " us\n";
		std::cout << "  Std Dev: " << stdDevTime << " us\n";*/

		std::cout << std::left << std::setw(60) << funcName
			<< std::right << std::setw(15) << stats.callCount
			<< std::right << std::setw(20) << totalTimeMs
			<< std::right << std::setw(20) << avgTimeMs
			<< std::right << std::setw(20) << minTime
			<< std::right << std::setw(20) << maxTime
			<< std::right << std::setw(20) << stdDevTime
			<< std::endl;

		/*if (stats.callCount <= 100) {
			for (size_t i = 0; i < stats.timeSamples.size(); ++i) {
				std::cout << stats.timeSamples[i] << " ";
			}
			std::cout << "\n";
		}*/
	}

	std::cout << std::string(175, '-') << std::endl;
}
