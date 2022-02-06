/*
 * Copyright (c) 2022 John Robinson.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

//
//  KdTree.h
//  KdTree
//
//  Created by John Robinson on 11/18/21.
//

#ifndef KdTree_h
#define KdTree_h
#include <random>
#include "kdTreeKmapKNlogn.h"

// a slight rewrite of the Romdomer class from
// https://stackoverflow.com/questions/13445688/how-to-generate-a-random-number-in-c/53887645#53887645
class RandomInterval {
  // random seed by default
  std::mt19937 gen_;
  std::uniform_int_distribution<int64_t> dist_;

public:
  RandomInterval(int64_t min, int64_t max, unsigned int seed = std::random_device{}())
    : gen_{ seed }, dist_{ min, max } {}

  // if you want predictable numbers
  void SetSeed(unsigned int seed) { gen_.seed(seed); }

  int64_t operator()() { return dist_(gen_); }
};

static std::vector<int64_t> getPermutations( size_t size, size_t dimension){
  // Determine the maximum depth of the k-d tree, which is log2(size).
  unsigned int maxDepth = 1;
  int64_t lsize = size;
  while (lsize > 0) {
    maxDepth++;
    lsize >>= 1;
  }
  // The partition coordinate permutes n the order 0, 1, 2, 3, 0, 1, 2, 3, etc.
  // for e.g. 4-dimensional data.
  std::vector<int64_t> permutation(maxDepth);
  for (size_t i = 0; i < permutation.size(); ++i) {
    permutation.at(i) = i % dimension;
  }
  return permutation;
}



template<typename K, typename V, size_t N>
class KdNodeExtras {

  typedef std::pair<std::vector<K>, std::list<V>> retPair_t;

private:
  std::vector<int64_t> permutation;

public:

  KdNodeExtras(size_t size) {
    permutation = getPermutations(size, N);

  }

  bool searchRegion(KdNode<K, V, N>* root, std::list<retPair_t>& result,
    std::vector<K>& queryLower, std::vector<K>& queryUpper, size_t maximumSubmitDepth, size_t size) {

    // Ensure that each query lower bound <= the corresponding query upper bound.
    for (size_t i = 0; i < queryLower.size(); ++i) {
      if (queryLower[i] > queryUpper[i]) {
        K tmp = queryLower[i];
        queryLower[i] = queryUpper[i];
        queryUpper[i] = tmp;
      }
    }

    // Search the tree and return the resulting list of KdNodes.
    regionSearch(root, result, queryLower, queryUpper, maximumSubmitDepth, 0);
    // copy KdNode list result to return pair
    typename std::list< KdNode<K, V, N>* >::iterator it;  
    return true;
  }
  
  bool searchRegionAndRemove(KdNode<K, V, N>* root, std::list<retPair_t>& result,
                             std::vector<K>& queryLower, std::vector<K>& queryUpper, signed_size_t maximumSubmitDepth, size_t size) {

    // Ensure that each query lower bound <= the corresponding query upper bound.
    for (size_t i = 0; i < queryLower.size(); ++i) {
      if (queryLower[i] > queryUpper[i]) {
        K tmp = queryLower[i];
        queryLower[i] = queryUpper[i];
        queryUpper[i] = tmp;
      }
    }

    // Search the tree and return the resulting list of KdNodes.
    auto retCode = regionSearchAndRemove(root, result, queryLower, queryUpper, maximumSubmitDepth, 0l);
    if (retCode == 1) root = nullptr;
    // copy the KdNode date ot a pair of tuple, valulist pair;
    typename std::list< KdNode<K, V, N> >::iterator it;
    return true;
  }
  
private:
  static int64_t pickValue(KdNode<K, V, N>* myThis, retPair_t& returnKV,
                       const uint64_t selector, bool removePick, const int depth){

      // init the return result to 0
      int64_t returnResult = 0;

      bool goGtThan = (selector & 0x1) == 1;

      // if greater than or less than, continue the search on  the appropriate child node
      if ((!goGtThan || myThis->gtChild == nullptr) && myThis->ltChild != nullptr) {
          returnResult = pickValue(myThis->ltChild, returnKV, selector >> 1, removePick, depth+1);
          // if the node below declared itself dead, remove the link
          if (removePick && returnResult == -1) {
            delete myThis->ltChild;
            myThis->ltChild = nullptr;
          }
      } else if((goGtThan || myThis->ltChild == nullptr) && myThis->gtChild != nullptr) {
          returnResult = pickValue(myThis->gtChild, returnKV, selector >> 1, removePick, depth+1);
          // if the node below declared itself dead, remove the link
          if (removePick && returnResult == -1) {
            delete myThis->gtChild;
            myThis->gtChild = nullptr;
          }
      } else {// both child pointers must be null to get here so this is a leaf node
          if (myThis->values != nullptr) {
            std::vector<K> key(N);
            for(size_t i = 0; i < N; i++){
                key[i] = myThis->tuple[i];
            }
            returnKV = retPair_t(key, *myThis->values);
            if (removePick) {
                delete myThis->values;
                myThis->values = nullptr ;
                returnResult = -1;  //flag for possible node removal
              }  else {
                returnResult = 1;
              }
          }
      }
      // Now figure out what to return.  If something was found either here or below return 1 or
      // if this node is still active as indicated by non null pointer set returnResult to 1.
      //Otherwise return a 0.
      if (returnResult == -1 && (myThis->values != nullptr || myThis->ltChild != nullptr || myThis->gtChild != nullptr)) {
          returnResult = 1;
      }
      return returnResult;
  }

      /**
   * <p>
   * The {@code pickValue) Picks an arbitrary value from the tree.
   * </p>
   *
       * @param key - array where the kee will be put.
       * @param query - Array containing the search point
       * @param valueToRemove - value to remove at that point
       * @returns a boolean which is true if that value was removed.
   */
public:
  static bool pickValue(retPair_t& returnKV, KdNode<K, V, N>* root, 
                        int selectionBias, bool remove) {

    // descent selector
    uint64_t selector = 0;
    switch (selectionBias) {
      case 0 : // less than side of the tree
        selector = 0L;
        break;
      case 1: // greater than side of the tree.
        selector = 0x7FFFFFFFFFFFFFFFUL;
        break;
      case 2: // middle case
        selector = 0x2AAAAAAAAAAAAAAAUL;
        break;
      case 3: // make a random selector
      {
        auto ri = RandomInterval(0, std::numeric_limits<int64_t>::max());
        selector = ri();
      }
        break;
      default:
        std::cout << "Selection Bias " << selectionBias << " not avaiable" << std::endl;
        return true;
    }
    
    // search the tree to find the node to return and possibly.
    int64_t returnResult = pickValue(root, returnKV, selector, remove, 0);
    // A result of 0 indicates that no value was found to pick which is an internal error or empty tree
    // -1 or +1 indicate a value was removed.
    if (returnResult == 0) return false;

    return true;
  }

  
  /**
   * <p>
   * The {@code searchKdTree} method searches the k-d tree and finds the KdNodes
   * that lie within a cutoff distance from a query node in all k dimensions.
   * </p>
   *
   * @param result - ArrayList to the k-d nodes that lie the query hypercube will be added.
   * @param queryPlus - Array containing the lager search bound for each dimension
   * @param queryMinus - Array containing the smaller search bound for each dimension
   * @param maximumSubmitDepth - the maximum tree depth at which a thread may be launched
   * @param depth - the depth in the k-d tree
   * @return void
   */
private:
  void regionSearch(KdNode<K, V, N>* myThis, std::list<retPair_t>& result,
                    const std::vector<K>& queryLower, const std::vector<K>& queryUpper,
                    const signed_size_t maximumSubmitDepth, const size_t depth) {

    // Look up the primary coordinate index.
    unsigned int p = (unsigned int)permutation.at(depth);

    // the branchCode will be used later to select the actual branch configuration in the switch statement
    // below.  0 = no branch, 1 = < branch only, 2 = > branch only, 3 = both branches.
    int branchCode = 0;

    // Search the < branch of the k-d tree if the partition coordinate of the queryPlus is
    // <= the partition coordinate of the k-d node.  The < branch
    // must be searched when the cutoff distance equals the partition coordinate because the super
    // key may assign a point to either branch of the tree if the sorting or partition coordinate,
    // which forms the most significant portion of the super key, shows equality.
    if (queryLower[p] <= myThis->tuple[p]) {
      // but only search if the ltChild pointer is not null;
      if (myThis->ltChild != nullptr) branchCode = 1;
      // Search the > branch of the k-d tree if the partition coordinate of the queryPlus is
      // >= the partition coordinate of the k-d node.  The < branch
      // must be searched when the cutoff distance equals the partition coordinate because the super
      // key may assign a point to either branch of the tree if the sorting or partition coordinate,
      // which forms the most significant portion of the super key, shows equality.
      if (queryUpper[p] >= myThis->tuple[p]) {
        // but only if the gtChild pointer is not null;
        if (myThis->gtChild != nullptr) branchCode += 2;
        // while here check to see if the local tuple is inside the the hypercube.
        if (myThis->values != nullptr){
            // If the distance from the query node to the k-d node is within the cutoff distance
            // in all k dimensions, add the k-d node to a list.
            bool inside = true;
            for (size_t i = 0; i < queryUpper.size(); i++) {
              if ((queryUpper[i] < myThis->tuple[i]) || (queryLower[i] > myThis->tuple[i]) ) {
                  inside = false;
                  break;
              }
            }
            if (inside) {
              std::vector<K> tempTuple(myThis->tuple, myThis->tuple + queryLower.size());
              retPair_t tmpPair(tempTuple, *myThis->values);
              result.push_back(tmpPair);
            }
        }
      }
    } else { // will not decend the lt branch so lets check the gt.
      if (myThis->gtChild != nullptr && queryUpper[p] >= myThis->tuple[p]) branchCode = 2;
    }

    // now implenent the branching decided on earlier
    switch (branchCode) {
      case 0: // child pointer are both null so just return
        break;
      case 1: // only go down the less than branch
        regionSearch(myThis->ltChild, result, queryLower, queryUpper, maximumSubmitDepth, depth + 1);
        break;
      case 2: // only go down the greater than branch
        regionSearch(myThis->gtChild, result, queryLower, queryUpper, maximumSubmitDepth, depth + 1);
        break;
      case 3: // go down both branches
        if(depth <= maximumSubmitDepth) {
          // get a future and another list ready for a child thread
          std::future< void > searchFuture ;
          std::list<retPair_t> ltResult;
          // Search the < branch asynchronously with a child thread.
          searchFuture = std::async(std::launch::async, [&] {
                         return regionSearch(myThis->ltChild, ltResult, queryLower, queryUpper,
                                                              maximumSubmitDepth, depth + 1); });
          // Search the > branch  with the master thread.
          regionSearch(myThis->gtChild, result, queryLower, queryUpper, maximumSubmitDepth, depth + 1);

          // Get the result of searching the < branch with the child thread.
          try {
            searchFuture.get();
          } catch (std::exception const& e) {
            std::cout << "caught exception " << e.what() << std::endl;
          }
          result.splice(result.end(), ltResult);
        } else { // if below the maximum submit depth
          regionSearch(myThis->ltChild, result, queryLower, queryUpper, maximumSubmitDepth, depth + 1);
          regionSearch(myThis->gtChild, result, queryLower, queryUpper, maximumSubmitDepth, depth + 1);
        }
        break;
    }
    return;
  }

  int64_t regionSearchAndRemove(KdNode<K, V, N>* myThis, std::list<retPair_t>& result,
                    const std::vector<K>& queryLower, const std::vector<K>& queryUpper,
                    const signed_size_t maximumSubmitDepth, const signed_size_t depth) {

    // Look up the primary coordinate index.
    unsigned int p = (unsigned int)permutation[(unsigned int)depth];

    // the branchCode will be used later to select the actual branch configuration in the switch statement
    // below.  0 = no branch, 1 = < branch only, 2 = > branch only, 3 = both branches.
    int branchCode = 0;
    // the return codes indicate whether the status of the node that was just returned from
    // 0 = active node
    // 1 = the node and all nodes below are dead and it's can be removed from the tree.
    int64_t ltRetCode = 0;
    int64_t gtRetCode = 0;

    // Search the < branch of the k-d tree if the partition coordinate of the queryPlus is
    // <= the partition coordinate of the k-d node.  The < branch
    // must be searched when the cutoff distance equals the partition coordinate because the super
    // key may assign a point to either branch of the tree if the sorting or partition coordinate,
    // which forms the most significant portion of the super key, shows equality.
    if (queryLower[p] <= myThis->tuple[p]) {
      // but only search if the ltChild pointer is not null;
      if (myThis->ltChild != nullptr) branchCode = 1;
      // Search the > branch of the k-d tree if the partition coordinate of the queryPlus is
      // >= the partition coordinate of the k-d node.  The < branch
      // must be searched when the cutoff distance equals the partition coordinate because the super
      // key may assign a point to either branch of the tree if the sorting or partition coordinate,
      // which forms the most significant portion of the super key, shows equality.
      if (queryUpper[p] >= myThis->tuple[p]) {
        // but only if the gtChild pointer is not null;
        if (myThis->gtChild != nullptr) branchCode += 2;
        // while here check to see if the local tuple is inside the the hypercube.
        if (myThis->values != nullptr){
            // If the distance from the query node to the k-d node is within the query rect
            // in all k dimensions, add the k-d node to a list.
            bool inside = true;
            for (size_t i = 0; i < queryUpper.size(); i++) {
              if ((queryUpper[i] < myThis->tuple[i]) || (queryLower[i] > myThis->tuple[i]) ) {
                  inside = false;
                  break;
              }
            }
            if (inside) {
              std::vector<K> tempTuple(myThis->tuple, myThis->tuple + queryLower.size());
              retPair_t tmpPair(tempTuple, *myThis->values);
              result.push_back(tmpPair);
              delete myThis->values;
              myThis->values = nullptr;   // mark the node dead by nulling the pointer
            }
        }
      }
    } else { // will not decend the lt branch so lets check the gt.
      if (myThis->gtChild != nullptr && queryUpper[p] >= myThis->tuple[p]) branchCode = 2;
    }

    // now implenent the branching decided on earlier
    switch (branchCode) {
      case 0: // child pointer are both null so just return
        break;
      case 1: // only go down the less than branch
        ltRetCode = regionSearchAndRemove(myThis->ltChild, result, queryLower, queryUpper, maximumSubmitDepth, depth + 1);
        break;
      case 2: // only go down the greater than branch
        gtRetCode = regionSearchAndRemove(myThis->gtChild, result, queryLower, queryUpper, maximumSubmitDepth, depth + 1);
        break;
      case 3: // go down both branches
        if(depth <= maximumSubmitDepth) {
          // get a future and another list ready for a child thread
          std::future< int64_t > searchFuture ;
          std::list<retPair_t> ltResult;
          // Search the < branch asynchronously with a child thread.
          searchFuture = std::async(std::launch::async, [&] {
                         return regionSearchAndRemove(myThis->ltChild, ltResult, queryLower, queryUpper,
                                                              maximumSubmitDepth, depth + 1); });
          // Search the > branch  with the master thread.
          gtRetCode = regionSearchAndRemove(myThis->gtChild, result, queryLower, queryUpper, maximumSubmitDepth, depth + 1);

          // Get the result of searching the < branch with the child thread.
          try {
            ltRetCode = searchFuture.get();
          } catch (std::exception const& e) {
            std::cout << "caught exception " << e.what() << std::endl;
          }
          result.splice(result.end(), ltResult);
        } else { // if below the maximum submit depth
          ltRetCode = regionSearchAndRemove(myThis->ltChild, result, queryLower, queryUpper, maximumSubmitDepth, depth + 1);
          gtRetCode = regionSearchAndRemove(myThis->gtChild, result, queryLower, queryUpper, maximumSubmitDepth, depth + 1);
        }
        break;
    }
    if (ltRetCode == 1) {
      delete myThis->ltChild;
      myThis->ltChild = nullptr;
    }
    if (gtRetCode == 1) {
      delete myThis->gtChild;
      myThis->gtChild = nullptr;
    }
    // retun a dead node code if all the pointers in this node are null;
    if (myThis->values == nullptr && myThis->ltChild == nullptr && myThis->gtChild == nullptr) {
      return 1;
    }
    return 0;
  }


}; // KdNodeExtras

template<typename K, typename V, size_t N>
class KdTree {

  typedef std::pair<std::vector<K>, std::list<V>> retPair_t;

private:
  size_t numPoints = 0;
  std::vector<KdNode<K, V, N>*> kdNodes;
  KdNode<K, V, N>* root = nullptr;
  KdNodeExtras<K, V, N>* kdNodeExtras = nullptr;

  size_t calcMaximumSubmitDepth(size_t numThreads) {
    size_t n = 0, maximumSubmitDepth;
    if (numThreads > 0) {
      while (numThreads > 0) {
        n++;
        numThreads >>= 1L;
      }
      numThreads = 1LL << (n - 1);
    }
    else {
      numThreads = 0;
    }
    size_t childThreads = numThreads - 1LL;
    maximumSubmitDepth = -1;
    if (numThreads < 2) {
      maximumSubmitDepth = -1; // The sentinel value -1 specifies no child threads.
    }
    else if (numThreads == 2) {
      maximumSubmitDepth = 0;
    }
    else {
      maximumSubmitDepth = (long)floor(log((double)childThreads) / log(2.));
    }
    return maximumSubmitDepth;
  }

public:
  // main constructor
  KdTree<K, V, N>() {
   }

  // add points and values to the KdTree.
  size_t add(std::vector<K>& tuple, V value) {
    if (tuple.size() != N) {
      return 0;
    }
    auto kp = new KdNode<K, V, N>(tuple, value);
    kdNodes.push_back(kp);
    numPoints = kdNodes.size();
    return numPoints;
  }

  bool searchRegion(std::list<retPair_t>& retPair, std::vector<K>& queryLower, std::vector<K>& queryUpper, size_t numThreads = 1) {
      // if the tree is not built yet, build it
      if (root == nullptr) {
          buildTree(numThreads);
      }
      size_t maximumSubmitDepth = calcMaximumSubmitDepth(numThreads);
      return kdNodeExtras->searchRegion(root, retPair, queryLower, queryUpper, maximumSubmitDepth, kdNodes.size());
  }

  bool searchRegion(std::list<V>& retVal, std::vector<K>& queryLower, std::vector<K>& queryUpper, size_t numThreads = 1) {
      // if the tree is not built yet, build it
      if (root == nullptr) {
          buildTree(numThreads);
      }
      std::list<retPair_t> retPair;
      size_t maximumSubmitDepth = calcMaximumSubmitDepth(numThreads);
      bool retFlag = kdNodeExtras->searchRegion(root, retPair, queryLower, queryUpper, maximumSubmitDepth, kdNodes.size());
      for (const retPair_t & pp : retPair) {
        retVal.splice(retVal.begin(), (std::list<V>)pp.second);
      }
      return retFlag;
  }

  bool searchRegionAndRemove(std::list<retPair_t>& retPair, std::vector<K>& queryLower, std::vector<K>& queryUpper, size_t numThreads = 1) {
    // if the tree is not built yet, build it
    if (root == nullptr) { buildTree(1); }
    
    signed_size_t maximumSubmitDepth = calcMaximumSubmitDepth(numThreads);
    kdNodeExtras->searchRegionAndRemove(root, retPair, queryLower, queryUpper, maximumSubmitDepth, kdNodes.size());
    return true;
  }

  bool buildTree(size_t numThreads) {
    size_t maximumSubmitDepth = calcMaximumSubmitDepth(numThreads);
    root = KdNode<K, V, N>::createKdTree(kdNodes, numThreads, maximumSubmitDepth);
    kdNodeExtras = new KdNodeExtras<K, V, N>(numPoints);
    return true;
  }

  bool pickValue(std::pair<std::vector<K>, std::list<V>>& returnKV, int selectionBias, bool remove)  {
  // if the tree is not built yet, build it
    if (root == nullptr) {
        buildTree(1);
        // if root is still null; return a null
        if (root == nullptr) return false;
    }
    return KdNodeExtras<K, V, N>::pickValue(returnKV, root, selectionBias, remove);
  }

  ~KdTree() {
    root->deleteKdTree();
    delete kdNodeExtras;
  }

}; // KdTree
