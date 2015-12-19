#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>  // NULL
#include <vector>
#include <iostream>
 
class profiler {
 public:
  // members
  struct timeval start;
  struct timeval end;
  double sum;
 
  // methods
  void Start();
  void Stop();
  void PauseSeconds();
  void PauseMicroseconds();
  void Continue();
  void StartPause();
  double GetTime();
  double GetTimeInSec();
  profiler() { sum = 0; }
};

// ProfilerCostas implementations.
void profiler::StartPause() {
  gettimeofday(&start, NULL);
  this->PauseMicroseconds();
}
 
void profiler::Start() {
  gettimeofday(&start, NULL);
}
 
void profiler::Stop() {
  gettimeofday(&end, NULL);
}
 
void profiler::PauseSeconds() {
  gettimeofday(&end, NULL);
 
  // seconds
  double tmpStart = start.tv_sec + (start.tv_usec/1000000.0);
  double tmpEnd = end.tv_sec + (end.tv_usec/1000000.0);
  double timeInSec = tmpEnd - tmpStart;
  sum += timeInSec;
}
 
void profiler::PauseMicroseconds() {
  gettimeofday(&end, NULL);
  // microseconds
  double tmpStart = start.tv_sec*1000000 + (start.tv_usec);
  double tmpEnd = end.tv_sec*1000000  + (end.tv_usec);
  double timeInSec = tmpEnd - tmpStart;
  sum += timeInSec;
}

void profiler::Continue() {
  gettimeofday(&start, NULL);
}
 
double profiler::GetTime() {
  return sum;
}

double profiler::GetTimeInSec() {
  double tmpStart = start.tv_sec + (start.tv_usec/1000000.0);
  double tmpEnd = end.tv_sec + (end.tv_usec/1000000.0);
  double timeInSec = tmpEnd - tmpStart;
  return timeInSec;
}
