## DBSCAN Clustering Implementation Using a K-d Tree

The code included here implements the popular DBSCAN clustering algorithm [1] in C++.  This implementation uses a K-d tree to accelerate the cluster-building process.  The k-d tree code is also presented here.  Using a k-d tree improves the clustering performance in two essential ways.  

First, the k-d tree code presented here is highly optimized to build and search the tree using multithreading, significantly reducing the time needed for window searches. 

Second, it supports removing items from the k-d tree after the items are added to a cluster.  Each clustered item is considered only once and does not have to be tagged as having been previously added to a cluster.  Also, the searches speed up as the k-d tree decreases in size.

Note that this is a very efficient implementation of a single-threaded DBSCAN algorithm.  There is a multithreaded version of the algorithm here: https://github.com/johnarobinson77/Multithreaded-DBSCAN-Clustering-in-Cpp.  But this is the better implementation if the use case needs to be single-threaded.

## Usage

This implementation assumes there is a set of objects that have a position in some n-space, geometric or otherwise.  The position is represented in a vector of K numbers and an object called a value here, which can be an object or pointer to an object.  Each position and value is presented as a pair to the DBSCAN object.  Also, a cluster window needs to be provided.  After running the clustering procedure, a list of clusters will be available where each cluster will be a list of the values presented to the DBSCAN object earlier. 

In the example below, that data is stored as a class in the Locations class.  The usage steps are as follows:



1. Create an instance of the DBSCAN class,
2. Loop through the positions and values, presenting them in pairs to the DBSCAN object
3. Call the build method with the build k-d tree and the search window
4. Call the checkCluster function to get the statistics of the clustering if desired.
5. Access the clusters list as needed.

This source code example below can be found in the main function in DBSCAN.cpp:


```c++
#include "dbscan.h"
	:
	:
 std::cout << "Adding data to DBSCAN..." << std::endl;
  auto dbscan = new DBSCAN< dkey_t, dval_t, numDimensions>;
  // Add each pair of points and values to the dbscan object
  for (size_t i = 0; i < numPoints; ++i) {
    dbscan->addPointWithValue(coordinates[i], i);
  }
  //Get the search range to about the cluster distance window
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
```



### Description of Uploaded Files

**DBSCAN.h** holds the DBSCAN class definition.  The DBSCAN algorithm is performed in the build function there. 

**KdTree.h**  contains the KdTree class.  This class is an API wrapper around the KdTreeExtras class and the KdNode class.  The KdTreeExtras class implements some special pick and search functions on the KdTree that are important to high performance.

**KdTreeKmanKnlogn.h** is where the KdNode class is defined.  The KdNode class implements the building of the KdTree using the Knlogn algorithm.  The algorithm description and original code can be found at [https://github.com/chezruss/kd-tree](https://github.com/chezruss/kd-tree). The version present here is an include file version that excludes the test code in the original version as well as a few modifications to make it compatible with the KdTree class.

**DBSCAN.cpp** provides an example and test case for the DBSCAN implementation. The test case is 1600 artificially clustered data with 400 points per cluster.  That data is presented to DBSCAN the results should therefore end up with a Clusters list that is 1600 in size, and each cluster should have 400 values.


## DBSCAN Class Methods


```c++
 // base constructor
  DBSCAN()

  // constructor with window
  DBSCAN(std::vector<K>& window)

  // addPointWithValuefuction adds a point or key and an associated value to the dbscan object
  size_t addPointWithValue(std::vector<K> point, V value)

  // setWindow sets the clustering widnow size
  void setWindow(std::vector<K>& window)

  // build builds the clusters.  All points and values must be added first.
  // It will return false if called a second time after dbscan was constructed.
  bool build()

  //sortClusterBySize sorts the cluster by the number of values in the cluster from largest to smallest
  void sortClustersBySize()

  // cluster setters and getters
  // getClusterSize returns the number of values in cluster clusterIdx
  size_t getClusterSize(size_t clusterIdx)

  // getClusterValueList returns the list of values in cluster clusterIdx
  std::list<V>* getClusterValueList(size_t clusterIdx) 

  // getClusterMaxCorner returns a vector of the maximum corner of the
  // bounding hypercube of cluster clusterIdx
  std::vector<K>* getClusterMaxCorner(size_t clusterIdx)

  // getClusterMinCorner returns a vector of the minimum corner of the
  // bounding hypercube of cluster clusterIdx
  std::vector<K>* getClusterMinCorner(size_t clusterIdx) 

  // getClusterTag returns the tag string set by setClusterTag
  std::string* getClusterTag(size_t clusterIdx)

  // cluster tag setter
  bool setClusterTag(size_t clusterIdx, std::string& tagIn)

  // getNumClusters returns the number of clusters
  size_t getNumClusters()

  // checCluster does some basic checks on the clusters and prints some statistics
  bool checkClusters(const size_t numLocations, const std::string* tag = nullptr)  

```



## Notes on Implementation

The DBSCAN clustering algorithm  works as follows:

1. Choose an arbitrary point from the dataset of items
2. Use it to seed a new cluster.
3. Loop on all items in the cluster
    1. Search for items in the dataset that are within a search window of the point
    2. Add any points returned from that search to the cluster and mark them in some way so that they are not included in future searches.
    3. When all the items in the current cluster have been searched against the dataset, exit the loop.  This means no more items in the dataset are within a search window of any of the points in this cluster.
4. Choose the next arbitrary point in the dataset that has not been added to a cluster and go back to step 2.  If there are no more points, the process is complete.

Using the k-d tree for the above algorithm helps in the following ways.



1. Searching the data set that is contained in the k-d tree is very efficient.  The k-d tree includes a function that returns all points in the tree, which is within a hypercube window, which is multithreaded, so it handles large trees very well.
2. There is a searchAndRemove method that removes the items found in the search from the k-d tree.  Other descriptions of DBSCAN algorithms, such as [https://en.wikipedia.org/wiki/DBSCAN](https://en.wikipedia.org/wiki/DBSCAN), talk about tagging each item added to a cluster so that it can be ignored if it is returned in a future dataset search.  But by deleting the item from the k-d tree, that item will never appear in a future search, so tagging is unnecessary.
3. 


## Differences from Original Paper

This implementation uses a hypercube search instead of a radial search to search for adjacent objects to add to a cluster.  An additional search method could be added to the KdTree class to test against the radial distance from the center instead of being inside the hypercube.  Because of the sum of squares calculation required for a radial distance test, this choice of search kernel will be slower.  Of course, it depends on the particular needs of the application.

Another difference is that this code does not explicitly implement the noise part of DBSCAN.  It would be easy to add a Min Cluster size parameter to the buildClusters() method.  Or a "noise" tag could be added to each cluster by looping over the clusters in the DBSCAN.clusters list.  A sort() method sorts clusters by size in the DBSCAN_Cluster class, which could aid in doing that tagging.

References

[1] Ester, Martin; Kriegel, Hans-Peter; Sander, Jörg; Xu, Xiaowei (1996). Simoudis, Evangelos; Han, Jiawei; Fayyad, Usama M. (eds.). A density-based algorithm for discovering clusters in large spatial databases with noise. Proceedings of the Second International Conference on Knowledge Discovery and Data Mining (KDD-96). AAAI Press. pp. 226–231. CiteSeerX 10.1.1.121.9220. ISBN 1-57735-004-9.

[2] Russell A. Brown, Building a Balanced k-d Tree in O(kn log n) Time, Journal of Computer Graphics Techniques (JCGT), vol. 4, no. 1, 50-68, 2015.

[3] Robert Sedgewick. Optimized Implementations, in Algorithms in C++, 173-174, Addison-Wesley, New York, 1992.
