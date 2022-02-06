#pragma once
/**
 * Copyright (c) 2022 John Robinson.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */


#ifndef DBSCAN_H
#define DBSCAN_H

#include <iostream>
#include <algorithm>
#include "KdTree.h"


template <typename K, typename V, size_t N>
class DBSCAN {

  // This is a shotened name of an std::pair that contains a tuple in the fist spot and a value in the second
  // WHen the KdTree search fucntions are called, a list of these pairs are returned.
  typedef std::pair<std::vector<K>, std::list<V>> retPair_t;

  // this struct maintains a bounding hyper rectangle
  struct bounds {
    std::vector<K> max;
    std::vector<K> min;
    void addToBounds(std::vector<K>& nv) {
      for (size_t i = 0; i < max.size(); i++) {
        max[i] = (nv[i] > max[i]) ? nv[i] : max[i];
        min[i] = (nv[i] < min[i]) ? nv[i] : min[i];
      }
    }
    void resetBounds() {
      for (size_t i = 0; i < max.size(); i++) {
        max[i] = std::numeric_limits<K>::min();
        min[i] = std::numeric_limits<K>::max();
      }
    }
    void resetBounds(int n) {
      max.resize(n);
      min.resize(n);
      resetBounds();
    }
    bounds(size_t n) {
      max.resize(n);
      min.resize(n);
      resetBounds();
    }
    bounds() = delete;
  };
  typedef bounds bounds_t;



  // Object of this structure contain the list of values of a single cluster
  struct Cluster {
    std::list<V>* values = nullptr;  // the calue list
    bounds_t clusterBounds{ N }; // the minimum corner of a bounding hyper rectangle
    std::string* tag = nullptr; // optional tag string

    // distructor
    ~Cluster() {
      if (values != nullptr) delete values;
      if (tag != nullptr) delete tag;
    }
  };

  // this is the list of clusters.
  std::vector<Cluster*>* clusters;
  // Pointer to the temporary KdTree that will be used to accelrated the clustering.
  KdTree< K, V, N>* kdTree;
  // a pointer to the current cluster being built.
  Cluster* current = nullptr;
  // this vector contains the distence to the edge of the search window in each dimention.
  std::vector<K> clusterWindowRadius;

  // this fuction creates an instence of a new Cluster object.
  bool newCluster() {
    auto cluster = new Cluster();
    cluster->values = new std::list<V>();
    cluster->tag = nullptr;
    current = cluster;
    clusters->push_back(cluster);
    return true;
  }

  bool addToCurrent(std::vector<K>& kdKey, std::list<V>& values) {
    current->values->splice(current->values->end(), values);
    current->clusterBounds.addToBounds(kdKey);
    return true;
  }

  bool buildCluster() {
    if (kdTree != nullptr && kdTree->buildTree(4) == false) return false;  // build a kdTree using 4 threads.
    std::vector<K> qLower(N);  //allocate 2 vectors for the upper and lower corners of the window
    std::vector<K> qUpper(N);
    auto keys = new std::list< std::vector<K> >();  // hold the list of the points to be searched for a cluster.
    retPair_t retPair; 
    std::list<retPair_t> retPairs;  // list of returned pairs from the kdTree Search
    while (kdTree->pickValue(retPair, 1, true)) {  // pick an initial cluster point until no more points
      newCluster();                 // create a new cluster
      keys->clear();
      addToCurrent(retPair.first, retPair.second); // add the picked Point to the cluster to the cluster
      keys->push_back(retPair.first);              // add the point to the keys list
      while (!keys->empty()) {                     // search around each point in the key list until empty
        std::vector<K> center = keys->front();
        keys->pop_front();
        for (size_t i = 0; i < N; i++) {           // build the search bounds
          qLower.at(i) = center[i] - clusterWindowRadius[i];
          qUpper.at(i) = center[i] + clusterWindowRadius[i];
        }
        retPairs.clear();
        kdTree->searchRegionAndRemove(retPairs, qLower, qUpper, 1); //search the tree for points within the window
        typename std::list< retPair_t >::iterator kit;
        for (kit = retPairs.begin(); kit != retPairs.end(); kit++) { // add the points to the keys list and cluster
          keys->push_back(kit->first);
          addToCurrent(kit->first, kit->second);
        }
      }
    }
    delete kdTree;
    kdTree = nullptr;
    delete keys;
    return true;
  }

public:
  // base constructor
  DBSCAN() {
    clusters = new std::vector<Cluster*>;
    kdTree = new KdTree<K, V, N>();
  }
  // constructor with window
  DBSCAN(std::vector<K>& window) {
    clusters = new std::vector<Cluster*>;
    kdTree = new KdTree<K, V, N>();
    setWindow(window);
  }

  ~DBSCAN() {
    for (size_t i = 0; i < clusters->size(); i++) {
      delete clusters->at(i);
    }
    delete clusters;
  }

  // addPoint fuction adds a point or key and an associated value to the dbscan object
  size_t addPointWithValue(std::vector<K> point, V value) {
    return kdTree->add(point, value);
  }

  // setWindow sets the clustering widnow size
  void setWindow(std::vector<K>& window) {
    clusterWindowRadius = window;
  }

  // build builds the clusters.  All points and vaules must be added first.
  // It will return false if called a second time after dbscan was constructed.
  bool build() {
    return buildCluster();
  }

  // sortClusterBySize sorts the cluster by number of values in the cluster 
  // from largest to smallest
  void sortClustersBySize() {
    std::sort(clusters->begin(), clusters->end(), 
      [](const Cluster* a, const Cluster* b) -> bool
      {
        return a->values->size() > b->values->size();
      });
  }

  // cluster setters and getters
  // getClusterSize returns the number of values in cluster clusterIdx
  size_t getClusterSize(size_t clusterIdx) {
    if (clusterIdx >= clusters->size()) return 0;
    return clusters->at(clusterIdx)->values->size();
  }
  // getClusterValueList returns the list of values in cluster clusterIdx
  std::list<V>* getClusterValueList(size_t clusterIdx) {
    if (clusterIdx >= clusters->size()) return nullptr;
    return (clusters->at(clusterIdx))->values;
  }
  // getClusterMaxCorner returns the a vector of the maximum corner of the
  // bounding hypercube of cluster clusterIdx
  std::vector<K>* getClusterMaxCorner(size_t clusterIdx) {
    if (clusterIdx >= clusters->size()) return nullptr;
    Cluster* cluster = clusters->at(clusterIdx);
    return &(cluster->clusterBounds.min);
  }
  // getClusterMinCorner returns the a vector of the minimum corner of the
  // bounding hypercube of cluster clusterIdx
  std::vector<K>* getClusterMinCorner(size_t clusterIdx) {
    if (clusterIdx >= clusters->size()) return nullptr;
    Cluster* cluster = clusters->at(clusterIdx);
    return &(cluster->clusterBounds.max);
  }
  // getClusterTag returns the tag string set by setClusterTag
  std::string* getClusterTag(size_t clusterIdx) {
    if (clusterIdx >= clusters->size()) return nullptr;
    return clusters[clusterIdx]->tag;
  }

  // cluster tag setter
  bool setClusterTag(size_t clusterIdx, std::string& tagIn) {
    clusters[clusterIdx]->tag = new std::string(tagIn);
  }


  // getNumClusters returns the number of clusters
  size_t getNumClusters() {
    return clusters->size();
  }

  // checCluster does some basic check on the clusters and print some statistics
  bool checkClusters(const size_t numLocations, const std::string* tag = nullptr) {
    size_t max = std::numeric_limits<size_t>::min();
    size_t min = std::numeric_limits<size_t>::max();
    size_t avgSize = 0;
    size_t count = 0;
    bool rb = true;
    for (size_t i = 0; i < clusters->size(); i++) {
      Cluster* c = clusters->at(i);
      //      if (tag == nullptr || c.tag != nullptr && c.hasTag(tag)) {
      if (tag == nullptr) {
        size_t t = c->values->size();
        if (t == 0) {
          std::cout << "Cluster " << i << " has 0 entries" << std::endl;
        }
        if (t > max) max = t;
        if (t < min) min = t;
        avgSize += t;
        count++;
      }
    }
    if (numLocations != -1 && avgSize != numLocations) {
      std::cout << "Number of locations in all clusters not equal input locations." << std::endl;
      rb = false;
    }
    avgSize = (count > 0) ? avgSize / count : 0;
    min = (count > 0) ? min : 0;
    max = (count > 0) ? max : 0;
    std::string t_tag = (tag == nullptr ? "" : *tag + " ");
    std::cout << "Cluster " << t_tag << "Count = " << count << " Max = " << max <<
      "  Min = " << min << " Average = " << avgSize << std::endl;
    return rb;
  }
};

#endif // DBSCAN_H
