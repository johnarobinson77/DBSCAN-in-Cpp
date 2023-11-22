/**
 * Copyright (c) 2022 John Robinson.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */


// dbscan.cpp : This file contains the a test of the DBSCAN code in the 'main' function. 


#include "dbscan.h"

// set these to the type to be used for the test.  dkey_t for the point coordinate type
typedef int64_t dkey_t;
typedef int64_t dval_t;

double timeNow() {
  struct timespec now;
  if (timespec_get(&now, TIME_UTC) != TIME_UTC) printf("timespec_get failure\n");
  double dblTime = (double)now.tv_sec + (double)now.tv_nsec * 1e-9;
  return dblTime;
}

struct timespec getTime(void) {
  struct timespec time;
  if (timespec_get(&time, TIME_UTC) != TIME_UTC) printf("timespec_get failure\n");
  return time;
}


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
  double startTime = timeNow();
  dbscan->build();
  double endTime = timeNow();

  double totalTime = endTime - startTime;
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

  delete dbscan;

  std::cout << "Starting sorting test of a small 5 cluster case." << std::endl;
  auto dbscan1 = new DBSCAN<dkey_t, dval_t, 1>();

  std::vector<dkey_t> smallc{ 1, 101, 201, 301, 401, 2, 102, 202, 302, 3, 103, 203, 4, 104, 5 };
  for (size_t i = 0; i < smallc.size(); i++) {
    dbscan1->addPointWithValue(std::vector<dkey_t>{smallc[i]}, i);
  }
  std::vector<dkey_t> window1{ 2 };
  dbscan1->setWindow(window1);
  dbscan1->build();
  dbscan1->sortClustersBySize();
  std::cout << "Sorted cluster sized are:";
  for (size_t i = 0; i < dbscan1->getNumClusters(); i++) {
    std::cout << " " << dbscan1->getClusterSize(i);
  }
  std::cout << std::endl;

  delete[] indices;
  delete dbscan1;
  exit(0);

}
