/**
 * Copyright (c) 2022 John Robinson.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

// dbscan.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#include "dbscan.h"

// set these to the type to be used for the test.  dkey_t for the point coordinate type
typedef int64_t dkey_t;
typedef int64_t dval_t;

/*
 * Create an alternate to clock_gettime(CLOCK_REALTIME, &time) for Mach. See
 * http://stackoverflow.com/questions/5167269/clock-gettime-alternative-in-mac-os-x
 * However, it appears that later versions of Mac OS X support clock_gettime(),
 * so this alternative may no longer be necessary for Mac OS X.
 */
#ifdef MACH
#include <mach/mach_time.h>

#define MACH_NANO (+1.0E-9)
#define MACH_GIGA UINT64_C(1000000000)

static double mach_timebase = 0.0;
static uint64_t mach_timestart = 0;

struct timespec getTime(void) {
  // be more careful in a multithreaded environement
  if (!mach_timestart) {
    mach_timebase_info_data_t tb = { 0, 1 }; // Initialize tb.numer and tb.denom
    mach_timebase_info(&tb);
    mach_timebase = tb.numer;
    mach_timebase /= tb.denom;
    mach_timestart = mach_absolute_time();
  }
  struct timespec t;
  double diff = (mach_absolute_time() - mach_timestart) * mach_timebase;
  t.tv_sec = diff * MACH_NANO;
  t.tv_nsec = diff - (t.tv_sec * MACH_GIGA);
  return t;
}
#else

#if defined(_WIN32) || defined(_WIN64)
 //see https://stackoverflow.com/questions/5404277/porting-clock-gettime-to-windows/5404467#5404467
#define NOMINMAX // Prevent Windows from getting confused about std::min vs. min, etc.
#include <windows.h>

int clock_gettime(int, struct timespec* tv)
{
  static int initialized = 0;
  static LARGE_INTEGER freq, startCount;
  static struct timespec tv_start;
  LARGE_INTEGER curCount;
  time_t sec_part;
  long nsec_part;

  if (!initialized) {
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&startCount);
    timespec_get(&tv_start, TIME_UTC);
    initialized = 1;
  }

  QueryPerformanceCounter(&curCount);

  curCount.QuadPart -= startCount.QuadPart;
  sec_part = curCount.QuadPart / freq.QuadPart;
  nsec_part = (long)((curCount.QuadPart - (sec_part * freq.QuadPart))
    * 1000000000UL / freq.QuadPart);

  tv->tv_sec = tv_start.tv_sec + sec_part;
  tv->tv_nsec = tv_start.tv_nsec + nsec_part;
  if (tv->tv_nsec >= 1000000000UL) {
    tv->tv_sec += 1;
    tv->tv_nsec -= 1000000000UL;
  }
  return 0;
}

#define CLOCK_REALTIME 0

#endif

struct timespec getTime(void) {
  struct timespec time;
  clock_gettime(CLOCK_REALTIME, &time);
  return time;
}
#endif


// The test case generated in main() creates artificially clustered data around 1600 random points 
// in n-space.  Each artificial cluster has 4000 points.  When DBSCAN is run on this data the
// checkClusters function should print out the following result:
// Cluster Count = 1600 Max = 4000  Min = 4000 Average = 4000
// This is followed by a quick check to see that all values were returned
// This test has been run with in64_t, float, and double for the dkey_t type and dimensions 1 to 5.
// Note: the randomInterval class used in this fuction is defined in KdTree.h

int main()
{
  std::cout << "Createing artificailly clustered data...." << std::endl;
  const int numDimensions = 3;
  const int clusterSpan = 1000;
  const int numClusters = 1600;
  const int numPointsPer = 4000;
  const dkey_t searchRad = (dkey_t)(clusterSpan * sqrt(numDimensions) / sqrt(3));
  auto riLarge = RandomInterval(-10000000000L, 10000000000L, 1);
  int64_t* clusterCenters = new int64_t[numClusters*numDimensions];
  for (size_t i = 0; i < numClusters; i++) {
    for (size_t k = 0; k < numDimensions; k++)
      clusterCenters[i * numDimensions + k] = riLarge();
  }

  // calculate the total number of points and create a coordiante array to hold them
  size_t numPoints = numClusters * numPointsPer;
  auto coordinates = new std::vector<dkey_t>[numPoints];
  // create a list that will be used later  check that all the right values have been returned.
  auto riSmall = RandomInterval(-clusterSpan, clusterSpan, 2);
  size_t idx = 0;
  bool* indices = new bool[numPoints];
  for (size_t j = 0; j < numPointsPer; j++) {
    for (size_t i = 0; i < numClusters; i++) {
      for (size_t k = 0; k < numDimensions; k++) {
        dkey_t tmp = riSmall();
        //std::cout << tmp << " ";
        coordinates[idx].push_back(clusterCenters[i * numDimensions + k] + tmp);
      }
      indices[idx] = false;
      idx++;
    }
  }

  std::cout << "Adding data to DBSCAN..." << std::endl;
  auto dbscan = new DBSCAN< dkey_t, dval_t, numDimensions>;
  // Add each pair points and values to the dbscan object
  for (size_t i = 0; i < numPoints; ++i) {
    dbscan->addPointWithValue(coordinates[i], i);
  }
  // get the search range to about cluster distance window
  std::vector<dkey_t> window(numDimensions);
  for (int i = 0; i < numDimensions; i++) {
    window[i] = searchRad;
  }

  std::cout << "Running DBSCAN clustering algorithm..." << std::endl;
  dbscan->setWindow(window);
  timespec startTime = getTime();
  dbscan->build();
  timespec endTime = getTime();

  double totalTime = (endTime.tv_sec - startTime.tv_sec) +
    1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));
  std::cout << "Total cluster time = " << std::fixed << std::setprecision(2) << totalTime << " seconds" << std::endl << std::endl;

  if (!dbscan->checkClusters(numPoints)) {
    exit(1);
  }

  std::cout << "Checking to see that all values ended up in some cluster" << std::endl;
  for (size_t i = 0; i < dbscan->getNumClusters(); i++) {
    for (dval_t val : *dbscan->getClusterValueList(i)) {
      indices[val] = true;
    }
  }

  bool passed = true;
  for (size_t i = 0; i < numPoints; i++) passed &= indices[i];
  if (passed) {
    std::cout << "Passed" << std::endl;
  } else {
    std::cout << "Failed to find all the values in the clusters." << std::endl;
  }

  dbscan->sortClustersBySize();

  delete[] indices;
  delete dbscan;
  exit(0);

}
