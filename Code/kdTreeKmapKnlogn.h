#pragma once
/*
 * Copyright (c) 2015, 2021 Russell A. Brown
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

 /*
  * The k-d tree was described by Jon Bentley in "Multidimensional Binary Search Trees
  * Used for Associative Searching", CACM 18(9): 509-517, 1975.  For k dimensions and
  * n elements of data, a balanced k-d tree is built in O(kn log n) + O((k-1)n log n)
  * time by first sorting the data in each of k dimensions, then building the k-d tree
  * in a manner that preserves the order of the k sorts while recursively partitioning
  * the data at each level of the k-d tree. No further sorting is necessary.
  *
  * Gnu g++ compilation options are: -lm -O3 -std=c++11 -DTEST_KD_TREE
  *
  * Optional compilation options are:
  *
  * -D K_DIMENSIONALITY=n - The k-dimensionality of the k-d tree (default 3).
  * -D INSERTION_SORT_CUTOFF=n - A cutoff for switching from merge sort to insertion sort
  *                              in the KdNode::mergeSort* functions (default 15)
  * -D NO_SUPER_KEY - Do not compare super-keys in the KdNode::regionSearch function.
  * -D MACH - Use a Mach equivalent to the clock_gettime(CLOCK_REALTIME, &time) function
  *           but this option appears to no longer be necessary.
  */

/* 
  Modifications from the original file made by John Robinson  2022.
  1. Seperated the class definigions from the test code in main and put them in this .h file
  2. Changed the method of removing duplicate nodes to adding to the value list in each node
      rather than createing a linklist of KdNodes.
  3. Put ifdefs around the code that prints out the execution time.
  4. Changed the uses of the signed_size_t type to size_t where appropriate.
  5. Added friend class reference to KdNodeExtras where other sort and access functions are defined.
  */

#include <limits>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <list>
#include <iostream>
#include <iomanip>
#include <exception>
#include <future>

/* A cutoff for switching from merge sort to insertion sort in the KdNode::mergeSort* functions */
#ifndef INSERTION_SORT_CUTOFF
#define INSERTION_SORT_CUTOFF 15
#endif

/*
 * This type is the signed equivalent of size_t and might be equivalent to intmax_t
 */
typedef int64_t signed_size_t;


/* A forward reference to NearestNeighborHeap */
template <typename, typename, size_t>
class NearestNeighborHeap;

/* A forward reference to KdNodeExtras */
template <typename, typename, size_t>
class KdNodeExtras;

/* One node of a k-d tree where K is key type, V is value type, and N is tuple dimensionality */
template <typename K, typename V, size_t N>
class KdNode {
public:
  K tuple[N];
private:
  KdNode<K, V, N>* ltChild;
  KdNode<K, V, N>* gtChild;
  std::list<V>* values;

public:
  KdNode(V const value) { // Pass non-primitive types as 'V const&'
    for (size_t i = 0; i < N; i++) this->tuple[i] = (K)0;
    ltChild = gtChild = nullptr; // redundant
    values = new std::list<V>();
    values->push_back(value);
  }
  KdNode(std::vector<K>& tuple, V const value) { // Pass non-primitive types as 'V const&'
    for (size_t i = 0; i < N; i++) this->tuple[i] = tuple[i];
    ltChild = gtChild = nullptr; // redundant
    values = new std::list<V>();
    values->push_back(value);
  }

public:
  ~KdNode() {
    // Delete each KdNode from the duplicates list.
    delete values;
  }

public:
  K const* getTuple() {
    return this->tuple;
  }

  friend class KdNodeExtras<K, V, N>;

  /*
   * The superKeyCompare function compares two K arrays in all k dimensions,
   * and uses the sorting or partition coordinate as the most significant dimension.
   *
   * Calling parameters:
   *
   * a - a K*
   * b - a K*
   * p - the most significant dimension
   *
   * returns a K result of comparing two K arrays
   */
private:
  inline
    static K superKeyCompare(K const* a, K const* b, size_t p) {
    // Typically, this first calculation of diff will be non-zero and bypass the 'for' loop.
    K diff = a[p] - b[p];
    for (size_t i = 1; diff == 0 && i < N; i++) {
      size_t r = i + p;
      // A fast alternative to the modulus operator for (i + p) < 2 * N.
      r = (r < N) ? r : r - N;
      diff = a[r] - b[r];
    }
    return diff;
  }

  /*
   * The following four merge sort functions are adapted from the mergesort function that is shown
   * on p. 166 of Robert Sedgewick's "Algorithms in C++", Addison-Wesley, Reading, MA, 1992.
   * That elegant implementation of the merge sort algorithm eliminates the requirement to test
   * whether the upper and lower halves of an auxiliary array have become exhausted during the
   * merge operation that copies from the auxiliary array to a result array.  This elimination is
   * made possible by inverting the order of the upper half of the auxiliary array and by accessing
   * elements of the upper half of the auxiliary array from highest address to lowest address while
   * accessing elements of the lower half of the auxiliary array from lowest address to highest
   * address.
   *
   * The following four merge sort functions also implement two suggestions from p. 275 of Robert
   * Sedgewick's and Kevin Wayne's "Algorithms 4th Edition", Addison-Wesley, New York, 2011.  The
   * first suggestion is to replace merge sort with insertion sort when the size of the array to
   * sort falls below a threshold.  The second suggestion is to avoid unnecessary copying to the
   * auxiliary array prior to the merge step of the algorithm by implementing two versions of
   * merge sort and by applying some "recursive trickery" to arrange that the required result is
   * returned in an auxiliary array by one version and in a result array by the other version.
   * The following four merge sort methods build upon this suggestion and return their result in
   * either ascending or descending order, as discussed on pp. 173-174 of Robert Sedgewick's
   * "Algorithms in C++", Addison-Wesley, Reading, MA, 1992.
   *
   * During multi-threaded execution, the upper and lower halves of the result array may be filled
   * from the auxiliary array (or vice versa) simultaneously by two threads.  The lower half of the
   * result array is filled by accessing elements of the upper half of the auxiliary array from highest
   * address to lowest address while accessing elements of the lower half of the auxiliary array from
   * lowest address to highest address, as explained above for elimination of the test for exhaustion.
   * The upper half of the result array is filled by addressing elements from the upper half of the
   * auxiliary array from lowest address to highest address while accessing the elements from the lower
   * half of the auxiliary array from highest address to lowest address.  Note: for the upper half
   * of the result array, there is no requirement to test for exhaustion provided that the upper half
   * of the result array never comprises more elements than the lower half of the result array.  This
   * provision is satisfied by computing the median address of the result array as shown below for
   * all four merge sort methods.
   *
   *
   * The mergeSortReferenceAscending function recursively subdivides the reference array then
   * merges the elements in ascending order and leaves the result in the reference array.
   *
   * Calling parameters:
   *
   * reference - a KdNode** array to sort via its (x, y, z, w...) tuples array
   * temporary - a KdNode** temporary array from which to copy sorted results;
   *             this array must be as large as the reference array
   * low - the start index of the region of the reference array to sort
   * high - the end index of the region of the reference array to sort
   * p - the sorting partition (x, y, z, w...)
   n  * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortReferenceAscending(KdNode<K, V, N>** reference, KdNode<K, V, N>** temporary,
    size_t low, size_t high, size_t p,
    signed_size_t maximumSubmitDepth, size_t depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the reference array with a child thread at as many levels of subdivision as possible.
      // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
      // Is a child thread available to subdivide the lower half of the reference array?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, recursively subdivide the lower half of the reference array with the current
        // thread and return the result in the temporary array in ascending order.
        mergeSortTemporaryAscending(reference, temporary, low, mid, p, maximumSubmitDepth, depth + 1);

        // Then recursively subdivide the upper half of the reference array with the current
        // thread and return the result in the temporary array in descending order.
        mergeSortTemporaryDescending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Compare the results in the temporary array in ascending order and merge them into
        // the reference array in ascending order.
        for (size_t i = low, j = high, k = low; k <= high; ++k) {
          reference[k] =
            (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p) < 0) ? temporary[i++] : temporary[j--];
        }

      }
      else {

        // Yes, a child thread is available, so recursively subdivide the lower half of the reference
        // array with a child thread and return the result in the temporary array in ascending order.
        std::future<void> sortFuture = std::async(std::launch::async, mergeSortTemporaryAscending, reference, temporary,
          low, mid, p, maximumSubmitDepth, depth + 1);

        // And simultaneously, recursively subdivide the upper half of the reference array with
        // the current thread and return the result in the temporary array in descending order.
        mergeSortTemporaryDescending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          sortFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }

        // Compare the results in the temporary array in ascending order with a child thread
        // and merge them into the lower half of the reference array in ascending order.
        std::future<void> mergeFuture =
          std::async(std::launch::async, [&] {
          for (size_t i = low, j = high, k = low; k <= mid; ++k) {
            reference[k] =
              (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p) <= 0)
              ? temporary[i++] : temporary[j--];
          }
            });

        // And simultaneously compare the results in the temporary array in descending order with the
        // current thread and merge them into the upper half of the reference array in ascending order.
        for (size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
          reference[k] =
            (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p) > 0) ? temporary[i--] : temporary[j++];
        }

        // Wait for the child thread to finish execution.
        try {
          mergeFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
      }

    }
    else {

      // Here is Jon Benley's implementation of insertion sort from "Programming Pearls", pp. 115-116,
      // Addison-Wesley, 1999, that sorts in ascending order and leaves the result in the reference array.
      for (size_t i = low + 1; i <= high; ++i) {
        KdNode<K, V, N>* tmp = reference[i];
        size_t j;
        for (j = i; j > low && superKeyCompare(reference[j - 1]->tuple, tmp->tuple, p) > 0; --j) {
          reference[j] = reference[j - 1];
        }
        reference[j] = tmp;
      }
    }
  }

  /*
   * The mergeSortReferenceDecending function recursively subdivides the reference array then
   * merges the elements in descending order and leaves the result in the reference array.
   *
   * Calling parameters:
   *
   * reference - a KdNode** array to sort via its (x, y, z, w...) tuples array
   * temporary - a KdNode** temporary array from which to copy sorted results;
   *             this array must be as large as the reference array
   * low - the start index of the region of the reference array to sort
   * high - the end index of the region of the reference array to sort
   * p - the sorting partition (x, y, z, w...)
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortReferenceDescending(KdNode<K, V, N>** reference, KdNode<K, V, N>** temporary,
    size_t low, size_t high, size_t p,
    signed_size_t maximumSubmitDepth, size_t depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the reference array with a child thread at as many levels of subdivision as possible.
      // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
      // Is a child thread available to subdivide the lower half of the reference array?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, recursively subdivide the lower half of the reference array with the current
        // thread and return the result in the temporary array in descending order.
        mergeSortTemporaryDescending(reference, temporary, low, mid, p, maximumSubmitDepth, depth + 1);

        // Then recursively subdivide the upper half of the reference array with the current
        // thread and return the result in the temporary array in ascending order.
        mergeSortTemporaryAscending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Compare the results in the temporary array in ascending order and merge them into
        // the reference array in descending order.
        for (size_t i = low, j = high, k = low; k <= high; ++k) {
          reference[k] =
            (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p) > 0) ? temporary[i++] : temporary[j--];
        }

      }
      else {

        // Yes, a child thread is available, so recursively subdivide the lower half of the reference
        // array with a child thread and return the result in the temporary array in descending order.
        std::future<void> sortFuture = std::async(std::launch::async, mergeSortTemporaryDescending, reference, temporary,
          low, mid, p, maximumSubmitDepth, depth + 1);

        // And simultaneously, recursively subdivide the upper half of the reference array with
        // the current thread and return the result in the temporary array in ascending order.
        mergeSortTemporaryAscending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          sortFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }

        // Compare the results in the temporary array in ascending order with a child thread
        // and merge them into the lower half of the reference array in descending order.
        std::future<void> mergeFuture =
          std::async(std::launch::async, [&] {
          for (size_t i = low, j = high, k = low; k <= mid; ++k) {
            reference[k] =
              (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p) >= 0)
              ? temporary[i++] : temporary[j--];
          }
            });

        // And simultaneously compare the results in the temporary array in descending order with the
        // current thread and merge them into the upper half of the reference array in descending order.
        for (size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
          reference[k] =
            (superKeyCompare(temporary[i]->tuple, temporary[j]->tuple, p) < 0) ? temporary[i--] : temporary[j++];
        }

        // Wait for the child thread to finish execution.
        try {
          mergeFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
      }

    }
    else {

      // Here is Jon Benley's implementation of insertion sort from "Programming Pearls", pp. 115-116,
      // Addison-Wesley, 1999, that sorts in descending order and leaves the result in the reference array.
      for (size_t i = low + 1; i <= high; ++i) {
        KdNode<K, V, N>* tmp = reference[i];
        size_t j;
        for (j = i; j > low && superKeyCompare(reference[j - 1]->tuple, tmp->tuple, p) < 0; --j) {
          reference[j] = reference[j - 1];
        }
        reference[j] = tmp;
      }
    }
  }

  /*
   * The mergeSortTemporaryAscending function recursively subdivides the reference array then
   * merges the elements in ascending order and leaves the result in the temporary array.
   *
   * Calling parameters:
   *
   * reference - a KdNode** array to sort via its (x, y, z, w...) tuples array
   * temporary - a KdNode** temporary array from which to copy sorted results;
   *             this array must be as large as the reference array
   * low - the start index of the region of the reference array to sort
   * high - the end index of the region of the reference array to sort
   * p - the sorting partition (x, y, z, w...)
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortTemporaryAscending(KdNode<K, V, N>** reference, KdNode<K, V, N>** temporary,
    size_t low, size_t high, size_t p,
    signed_size_t maximumSubmitDepth, size_t depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the reference array with a child thread at as many levels of subdivision as possible.
      // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
      // Is a child thread available to subdivide the lower half of the reference array?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, recursively subdivide the lower half of the reference array with the current
        // thread and return the result in the reference array in ascending order.
        mergeSortReferenceAscending(reference, temporary, low, mid, p, maximumSubmitDepth, depth + 1);

        // Then recursively subdivide the upper half of the reference array with the current
        // thread and return the result in the reference array in descending order.
        mergeSortReferenceDescending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Compare the results in the reference array in ascending order and merge them into
        // the temporary array in ascending order.
        for (size_t i = low, j = high, k = low; k <= high; ++k) {
          temporary[k] =
            (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p) < 0) ? reference[i++] : reference[j--];
        }

      }
      else {

        // Yes, a child thread is available, so recursively subdivide the lower half of the reference
        // array with a child thread and return the result in the reference array in ascending order.
        std::future<void> sortFuture = std::async(std::launch::async, mergeSortReferenceAscending, reference, temporary,
          low, mid, p, maximumSubmitDepth, depth + 1);

        // And simultaneously, recursively subdivide the upper half of the reference array with
        // the current thread and return the result in the reference array in descending order.
        mergeSortReferenceDescending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          sortFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }

        // Compare the results in the reference array in ascending order with a child thread
        // and merge them into the lower half of the temporary array in ascending order.
        std::future<void> mergeFuture =
          std::async(std::launch::async, [&] {
          for (size_t i = low, j = high, k = low; k <= mid; ++k) {
            temporary[k] =
              (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p) <= 0)
              ? reference[i++] : reference[j--];
          }
            });

        // And simultaneously compare the results in the reference array in descending order with the
        // current thread and merge them into the upper half of the temporary array in ascending order.
        for (size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
          temporary[k] =
            (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p) > 0) ? reference[i--] : reference[j++];
        }

        // Wait for the child thread to finish execution.
        try {
          mergeFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
      }

    }
    else {

      // Here is John Robinson's implementation of insertion sort that sorts in ascending order
      // and leaves the result in the temporary array.
      temporary[high] = reference[high];
      signed_size_t i;
      signed_size_t j; // MUST be signed because it can decrement to -1
      for (j = high - 1; j >= (signed_size_t)low; --j) {
        for (i = j; i < (signed_size_t)high; ++i) {
          if (superKeyCompare(reference[j]->tuple, temporary[i + 1]->tuple, p) > 0) {
            temporary[i] = temporary[i + 1];
          }
          else {
            break;
          }
        }
        temporary[i] = reference[j];
      }
    }
  }

  /*
   * The mergeSortTemporaryDecending function recursively subdivides the reference array
   * then merges the elements in descending order and leaves the result in the reference array.
   *
   * Calling parameters:
   *
   * reference - a KdNode** array to sort via its (x, y, z, w...) tuples array
   * temporary - a KdNode** temporary array from which to copy sorted results;
   *             this array must be as large as the reference array
   * low - the start index of the region of the reference array to sort
   * high - the end index of the region of the reference array to sort
   * p - the sorting partition (x, y, z, w...)
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the tree depth
   */
private:
  static void mergeSortTemporaryDescending(KdNode<K, V, N>** reference, KdNode<K, V, N>** temporary,
    size_t low, size_t high, size_t p,
    signed_size_t maximumSubmitDepth, size_t depth) {

    if (high - low > INSERTION_SORT_CUTOFF) {

      // Avoid overflow when calculating the median.
      size_t const mid = low + ((high - low) >> 1);

      // Subdivide the lower half of the reference array with a child thread at as many levels of subdivision as possible.
      // Create the child threads as high in the subdivision hierarchy as possible for greater utilization.
      // Is a child thread available to subdivide the lower half of the reference array?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, recursively subdivide the lower half of the reference array with the current
        // thread and return the result in the reference array in descending order.
        mergeSortReferenceDescending(reference, temporary, low, mid, p, maximumSubmitDepth, depth + 1);

        // Then recursively subdivide the upper half of the reference array with the current
        // thread and return the result in the reference array in ascending order.
        mergeSortReferenceAscending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Compare the results in the reference array in ascending order and merge them into
        // the temporary array in descending order.
        for (size_t i = low, j = high, k = low; k <= high; ++k) {
          temporary[k] =
            (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p) > 0) ? reference[i++] : reference[j--];
        }

      }
      else {

        // Yes, a child thread is available, so recursively subdivide the lower half of the reference
        // array with a child thread and return the result in the reference array in descending order.
        std::future<void> sortFuture = std::async(std::launch::async, mergeSortReferenceDescending, reference, temporary,
          low, mid, p, maximumSubmitDepth, depth + 1);

        // And simultaneously, recursively subdivide the upper half of the reference array with
        // the current thread and return the result in the reference array in ascending order.
        mergeSortReferenceAscending(reference, temporary, mid + 1, high, p, maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          sortFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }

        // Compare the results in the reference array in ascending order with a child thread
        // and merge them into the lower half of the temporary array in descending order.
        std::future<void> mergeFuture =
          std::async(std::launch::async, [&] {
          for (size_t i = low, j = high, k = low; k <= mid; ++k) {
            temporary[k] =
              (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p) >= 0)
              ? reference[i++] : reference[j--];
          }
            });

        // And simultaneously compare the results in the reference array in descending order with the
        // current thread and merge them into the upper half of the temporary array in descending order.
        for (size_t i = mid, j = mid + 1, k = high; k > mid; --k) {
          temporary[k] =
            (superKeyCompare(reference[i]->tuple, reference[j]->tuple, p) < 0) ? reference[i--] : reference[j++];
        }

        // Wait for the child thread to finish execution.
        try {
          mergeFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
      }

    }
    else {

      // Here is John Robinson's implementation of insertion sort that sorts in descending order
      // and leaves the result in the temporary array.
      temporary[high] = reference[high];
      size_t i;
      size_t j; // MUST be signed because it can decrement to -1
      for (j = high - 1; j >= low; --j) {
        for (i = j; i < high; ++i) {
          if (superKeyCompare(reference[j]->tuple, temporary[i + 1]->tuple, p) < 0) {
            temporary[i] = temporary[i + 1];
          }
          else {
            break;
          }
        }
        temporary[i] = reference[j];
      }
    }
  }

  /*
   * The removeDuplicates function checks the validity of the merge sort and
   * removes duplicates from the kdNodes array.
   *
   * Calling parameters:
   *
   * kdNodes - a KdNode** array that has been sorted via merge sort according to (x,y,z,w...) tuples
   * i - the leading dimension for the super key
   *
   * returns the end index of the reference array following removal of duplicate elements
   */
private:
  inline
    static size_t removeDuplicates(KdNode<K, V, N>** kdNodes, size_t i, size_t size) {
    size_t end = 0;
    for (size_t j = 1; j < size; ++j) {
      K compare = superKeyCompare(kdNodes[j]->tuple, kdNodes[end]->tuple, i);
      if (compare < 0) {
        std::cout << "merge sort failure: superKeyCompare(kdNodes[" << j << "], kdNodes["
          << end << "], " << i << ") = " << compare << std::endl;
        exit(1);
      }
      else if (compare > 0) {
        // Keep the jth element of the kdNodes array.
        kdNodes[++end] = kdNodes[j];
      }
      else {
        // append vales to the end of fte kpet kdNode and delete the now unused node
        kdNodes[end]->values->splice(kdNodes[end]->values->end(), *kdNodes[j]->values);
        delete kdNodes[j];
      }
    }
    return end;
  }

  /*
   * The buildKdTree function builds a k-d tree by recursively partitioning
   * the reference arrays and adding KdNodes to the tree.  These arrays
   * are permuted cyclically for successive levels of the tree in
   * order that sorting occur in the order x, y, z, w...
   *
   * Calling parameters:
   *
   * reference - a KdNode*** array to recursively sort via its (x, y, z, w...) tuples array
   * temporary - a KdNode*** temporary array from which to copy sorted results;
   * start - start element of the reference array
   * end - end element of the reference array
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the depth in the tree
   *
   * returns: a KdNode pointer to the root of the k-d tree
   */
private:
  static KdNode<K, V, N>* buildKdTree(KdNode<K, V, N>*** references,
    std::vector< std::vector<size_t> > const& permutation,
    size_t start, size_t end,
    signed_size_t maximumSubmitDepth, size_t depth) {

    KdNode<K, V, N>* node = nullptr;

    // The partition permutes as x, y, z, w... and specifies the most significant key.
    size_t p = permutation.at(depth).at(permutation.at(0).size() - 1);

    // Obtain the reference array that corresponds to the most significant key.
    KdNode<K, V, N>** reference = references[permutation.at(depth).at(N)];

    if (end == start) {

      // Only one reference was passed to this function, so add it to the tree.
      node = reference[end];

    }
    else if (end == start + 1) {

      // Two references were passed to this function in sorted order, so store the start
      // element at this level of the tree and store the end element as the > child.
      node = reference[start];
      node->gtChild = reference[end];

    }
    else if (end == start + 2) {

      // Three references were passed to this function in sorted order, so
      // store the median element at this level of the tree, store the start
      // element as the < child and store the end element as the > child.
      node = reference[start + 1];
      node->ltChild = reference[start];
      node->gtChild = reference[end];

    }
    else if (end > start + 2) {

      // Four or more references were passed to this function, so the
      // median element of the reference array is chosen as the tuple
      // about which the other reference arrays will be partitioned
      // Avoid overflow when computing the median.
      size_t median = start + ((end - start) / 2);

      // Store the median element of the reference array in a new KdNode.
      node = reference[median];

      // Build both branches with child threads at as many levels of the tree
      // as possible.  Create the child threads as high in the tree as possible.
      // Are child threads available to build both branches of the tree?
      if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

        // No, child threads are not available, so one thread will be used.
        // Initialize startIndex=1 so that the 'for' loop that partitions the
        // reference arrays will partition a number of arrays equal to N.
        size_t startIndex = 1;

        // If depth < N-1, copy references[permut[N]] to references[permut[0]]
        // where permut is the permutation vector for this level of the tree.
        // Sort the two halves of references[permut[0]] with p+1 as the most
        // significant key of the super key. Use as the temporary array
        // references[permut[1]] because that array is not used for partitioning.
        // Partition a number of reference arrays equal to the tree depth because
        // those reference arrays are already sorted.
        if (depth < N - 1) {
          startIndex = N - depth;
          KdNode<K, V, N>** dst = references[permutation.at(depth).at(0)];
          KdNode<K, V, N>** tmp = references[permutation.at(depth).at(1)];
          for (size_t i = start; i <= end; ++i) {
            dst[i] = reference[i];
          }
          // Sort the lower half of references[permut[0]] with the current thread.
          mergeSortReferenceAscending(dst, tmp, start, median - 1, p + 1, maximumSubmitDepth, depth);
          // Sort the upper half of references[permut[0]] with the current thread.
          mergeSortReferenceAscending(dst, tmp, median + 1, end, p + 1, maximumSubmitDepth, depth);
        }

        // Partition the reference arrays specified by 'startIndex' in
        // a priori sorted order by comparing super keys.  Store the
        // result from references[permut[i]]] in references[permut[i-1]]
        // where permut is the permutation vector for this level of the
        // tree, thus permuting the reference arrays. Skip the element
        // of references[permut[i]] that equals the tuple that is stored
        // in the new KdNode.
        K* tuple = node->tuple;
        for (size_t i = startIndex; i < N; ++i) {
          // Specify the source and destination reference arrays.
          KdNode<K, V, N>** src = references[permutation.at(depth).at(i)];
          KdNode<K, V, N>** dst = references[permutation.at(depth).at(i - 1)];

          // Fill the lower and upper halves of one reference array
          // in ascending order with the current thread.
          for (size_t j = start, lower = start - 1, upper = median; j <= end; ++j) {
            KdNode<K, V, N>* src_j = src[j];
            K compare = superKeyCompare(src_j->tuple, tuple, p);
            if (compare < 0) {
              dst[++lower] = src_j;
            }
            else if (compare > 0) {
              dst[++upper] = src_j;
            }
          }
        }

        // Recursively build the < branch of the tree with the current thread.
        node->ltChild = buildKdTree(references, permutation, start, median - 1,
          maximumSubmitDepth, depth + 1);

        // Then recursively build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(references, permutation, median + 1, end,
          maximumSubmitDepth, depth + 1);

      }
      else {

        // Yes, child threads are available, so two threads will be used.
        // Initialize endIndex=0 so that the 'for' loop that partitions the
        // reference arrays will partition a number of arrays equal to N.
        size_t startIndex = 1;

        // If depth < N-1, copy references[permut[N]] to references[permut[0]]
        // where permut is the permutation vector for this level of the tree.
        // Sort the two halves of references[permut[0]] with p+1 as the most
        // significant key of the super key. Use as the temporary array
        // references[permut[1]] because that array is not used for partitioning.
        // Partition a number of reference arrays equal to the tree depth because
        // those reference arrays are already sorted.
        if (depth < N - 1) {
          startIndex = N - depth;
          KdNode<K, V, N>** dst = references[permutation.at(depth).at(0)];
          KdNode<K, V, N>** tmp = references[permutation.at(depth).at(1)];
          // Copy and sort the lower half of references[permut[0]] with a child thread.
          std::future<void> copyFuture =
            std::async(std::launch::async, [&] {
            for (size_t i = start; i <= median - 1L; ++i) {
              dst[i] = reference[i];
            }
            mergeSortReferenceAscending(dst, tmp, start, median - 1, p + 1, maximumSubmitDepth, depth);
              });

          // Copy and sort the upper half of references[permut[0]] with the current thread.
          for (size_t i = median + 1; i <= end; ++i) {
            dst[i] = reference[i];
          }
          mergeSortReferenceAscending(dst, tmp, median + 1, end, p + 1, maximumSubmitDepth, depth);

          // Wait for the child thread to finish execution.
          try {
            copyFuture.get();
          }
          catch (std::exception const& e) {
            std::cout << "caught exception " << e.what() << std::endl;
          }
        }

        // Create a copy of the node->tuple array so that the current thread
        // and the child thread do not contend for read access to this array.
        K* tuple = node->tuple;
        K* point = new K[N];
        for (size_t i = 0; i < N; ++i) {
          point[i] = tuple[i];
        }

        // Partition the reference arrays specified by 'startIndex' in
        // a priori sorted order by comparing super keys.  Store the
        // result from references[permut[i]]] in references[permut[i-1]]
        // where permut is the permutation vector for this level of the
        // tree, thus permuting the reference arrays. Skip the element
        // of references[permut[i]] that equals the tuple that is stored
        // in the new KdNode.
        for (size_t i = startIndex; i < N; ++i) {
          // Specify the source and destination reference arrays.
          KdNode<K, V, N>** src = references[permutation.at(depth).at(i)];
          KdNode<K, V, N>** dst = references[permutation.at(depth).at(i - 1)];

          // Two threads may be used to partition the reference arrays, analogous to
          // the use of two threads to merge the results for the merge sort algorithm.
          // Fill one reference array in ascending order with a child thread.
          std::future<void> partitionFuture =
            std::async(std::launch::async, [&] {
            for (size_t lower = start - 1, upper = median, j = start; j <= median; ++j) {
              KdNode<K, V, N>* src_j = src[j];
              K compare = superKeyCompare(src_j->tuple, point, p);
              if (compare < 0) {
                dst[++lower] = src_j;
              }
              else if (compare > 0) {
                dst[++upper] = src_j;
              }
            }
              });

          // Simultaneously fill the same reference array in descending order with the current thread.
          for (size_t lower = median, upper = end + 1, k = end; k > median; --k) {
            KdNode<K, V, N>* src_k = src[k];
            K compare = superKeyCompare(src_k->tuple, tuple, p);
            if (compare < 0) {
              dst[--lower] = src_k;
            }
            else if (compare > 0) {
              dst[--upper] = src_k;
            }
          }

          // Wait for the child thread to finish execution.
          try {
            partitionFuture.get();
          }
          catch (std::exception const& e) {
            std::cout << "caught exception " << e.what() << std::endl;
          }
        }

        // Delete the point array.
        delete[] point;

        // Recursively build the < branch of the tree with a child thread.
        // The recursive call to buildKdTree must be placed in a lambda
        // expression because buildKdTree is a template not a function.
        std::future<KdNode<K, V, N>*> buildFuture =
          std::async(std::launch::async, [&] {
          return buildKdTree(references, permutation, start, median - 1,
            maximumSubmitDepth, depth + 1);
            });

        // And simultaneously build the > branch of the tree with the current thread.
        node->gtChild = buildKdTree(references, permutation, median + 1, end,
          maximumSubmitDepth, depth + 1);

        // Wait for the child thread to finish execution.
        try {
          node->ltChild = buildFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
      }

    }
    else if (end < start) {

      // This is an illegal condition that should never occur, so test for it last.
      std::cout << "error has occurred at depth = " << depth << " : end = " << end
        << "  <  start = " << start << std::endl;
      exit(1);

    }

    // Return the pointer to the root of the k-d tree.
    return node;
  }

  /*
   * The verifyKdTree function walks the k-d tree and checks that the
   * children of a node are in the correct branch of that node.
   *
   * Calling parameters:
   *
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the depth in the k-d tree
   *
   * returns: a count of the number of kdNodes in the k-d tree
   */
private:
  size_t verifyKdTree(signed_size_t maximumSubmitDepth, size_t depth) {

    size_t count = 1;

    // The partition cycles as x, y, z, w...
    size_t p = depth % N;

    if (ltChild != nullptr) {
      if (ltChild->tuple[p] > tuple[p]) {
        std::cout << "child is > node!" << std::endl;
        exit(1);
      }
      if (superKeyCompare(ltChild->tuple, tuple, p) >= 0) {
        std::cout << "child is >= node!" << std::endl;
        exit(1);
      }
    }
    if (gtChild != nullptr) {
      if (gtChild->tuple[p] < tuple[p]) {
        std::cout << "child is < node!" << std::endl;
        exit(1);
      }
      if (superKeyCompare(gtChild->tuple, tuple, p) <= 0) {
        std::cout << "child is <= node" << std::endl;
        exit(1);
      }
    }

    // Verify the < branch with a child thread at as many levels of the tree as possible.
    // Create the child thread as high in the tree as possible for greater utilization.

    // Is a child thread available to build the < branch?
    if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

      // No, so verify the < branch with the current thread.
      if (ltChild != nullptr) {
        count += ltChild->verifyKdTree(maximumSubmitDepth, depth + 1);
      }

      // Then verify the > branch with the current thread.
      if (gtChild != nullptr) {
        count += gtChild->verifyKdTree(maximumSubmitDepth, depth + 1);
      }
    }
    else {

      // Yes, so verify the < branch with a child thread. Note that a
      // lambda is required to instantiate the verifyKdTree template.
      std::future<size_t> verifyFuture;
      if (ltChild != nullptr) {
        verifyFuture =
          std::async(std::launch::async, [&] {
          return ltChild->verifyKdTree(maximumSubmitDepth, depth + 1);
            });
      }

      // And simultaneously verify the > branch with the current thread.
      size_t gtCount = 0;
      if (gtChild != nullptr) {
        gtCount = gtChild->verifyKdTree(maximumSubmitDepth, depth + 1);
      }

      // Wait for the child thread to finish execution.
      size_t ltCount = 0;
      if (ltChild != nullptr) {
        try {
          ltCount = verifyFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
      }
      count += ltCount + gtCount;
    }

    return count;
  }

  /*
   * The swap function swaps two elements in a vector<size_t>.
   *
   * Calling parameters:
   *
   * a - the vector
   * i - the index of the first element
   * j - the index of the second element
   */
private:
  inline
    static void swap(std::vector<size_t>& a, size_t i, size_t j) {
    size_t t = a.at(i);
    a.at(i) = a.at(j);
    a.at(j) = t;
  }

  /*
   * The createKdTree function performs the necessary initialization then calls the buildKdTree function.
   *
   * Calling parameters:
   *
   * kdNodes - a vector<KdNode*> wherein each KdNode contains a (x,y,z,w...) tuple
   * numDimensions - the number of dimensions
   * numThreads - the number of threads
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   *
   * returns: a KdNode pointer to the root of the k-d tree
   */
public:
  static KdNode<K, V, N>* createKdTree(std::vector<KdNode<K, V, N>*>& kdNodes,
    size_t numThreads, signed_size_t maximumSubmitDepth) {

#ifdef TEST_KD_TREE
    struct timespec startTime, endTime;
#endif //TEST_KD_TREE

    // Create the references arrays including one additional array for use in building the k-d tree.
    KdNode<K, V, N>*** references = new KdNode<K, V, N>**[N + 1];

    // The first references array is the .data() array of the kdNodes vector.
    references[0] = kdNodes.data();

    // Allocate the remaining references arrays.
    for (size_t i = 1; i < N + 1; ++i) {
      references[i] = new KdNode<K, V, N>*[kdNodes.size()];
    }

    // Sort the first reference array using multiple threads. Importantly,
    // for compatibility with the 'permutation' vector initialized below,
    // use the first dimension (0) as the leading key of the super key.
#ifdef TEST_KD_TREE
    startTime = getTime();
#endif //TEST_KD_TREE
    mergeSortReferenceAscending(references[0], references[N], 0, kdNodes.size() - 1,
      0, maximumSubmitDepth, 0);
#ifdef TEST_KD_TREE
    endTime = getTime();
    double sortTime = (endTime.tv_sec - startTime.tv_sec) +
      1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));
#endif//TEST_KD_TREE
    // Remove references to duplicate coordinates via one pass through the first reference array.
#ifdef TEST_KD_TREE
    startTime = getTime();
#endif //TEST_KD_TREE
    size_t end = removeDuplicates(references[0], 0, kdNodes.size());
#ifdef TEST_KD_TREE
    endTime = getTime();
    double removeTime = (endTime.tv_sec - startTime.tv_sec) +
      1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));
#endif//TEST_KD_TREE

#ifdef TEST_KD_TREE
    // Start the timer to time building the k-d tree.
    startTime = getTime();
#endif //TEST_KD_TREE

    // Determine the maximum depth of the k-d tree, which is log2( kdNodes.size() ).
    size_t size = kdNodes.size();
    size_t maxDepth = 1;
    while (size > 0) {
      ++maxDepth;
      size >>= 1;
    }

    // It is unnecessary to compute either the permutation of the reference array or
    // the partition coordinate upon each recursive call of the buildKdTree function
    // because both depend only on the depth of recursion, so they may be pre-computed.
    // Create and initialize an 'indices' vector for the permutation calculation.
    // Because this vector is initialized with 0, 1, 2, 3, 0, 1, 2, 3, etc. (for
    // e.g. 4-dimensional data), the leading key of the super key will be 0 at the
    // first level of the nascent tree, consistent with having sorted the reference
    // array above using 0 as the leading key of the super key.
    std::vector<size_t> indices(N + 2);
    for (size_t i = 0; i < indices.size() - 1; ++i) {
      indices[i] = i;
    }

    // Create a 'permutation' vector from the 'indices' vector to specify permutation
    // of the reference arrays and of the partition coordinate.
    std::vector< std::vector<size_t> > permutation(maxDepth, std::vector<size_t>(N + 2));

    // Fill the permutation vector by calculating the permutation of the indices vector
    // and the the partition coordinate of the tuple at each depth in the tree.
    for (size_t i = 0; i < permutation.size(); ++i) {
      // The last entry of the indices vector contains the partition coordinate.
      indices.at(N + 1) = i % N;
      // Swap the first and second to the last elements of the indices vector.
      swap(indices, 0, N);
      // Copy the indices vector to one row of the permutation vector.
      permutation.at(i) = indices;
      // Swap the third and second to the last elements of the indices vector.
      swap(indices, N - 1, N);
    }

    // Build the k-d tree with multiple threads if possible.
    KdNode<K, V, N>* root = buildKdTree(references, permutation, 0, end,
      maximumSubmitDepth, 0);
#ifdef TEST_KD_TREE
    endTime = getTime();
    double kdTime = (endTime.tv_sec - startTime.tv_sec) +
      1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));
#endif //TEST_KD_TREE
    // Verify the k-d tree and report the number of kdNodes.
#ifdef TEST_KD_TREE
    startTime = getTime();
#endif //TEST_KD_TREE
    size_t numberOfNodes;
    numberOfNodes = root->verifyKdTree(maximumSubmitDepth, 0);
#ifdef TEST_KD_TREE
    endTime = getTime();
    double verifyTime = (endTime.tv_sec - startTime.tv_sec) +
      1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));
    std::cout << "Number of nodes = " << numberOfNodes << std::endl << std::endl;

    std::cout << "totalTime = " << std::fixed << std::setprecision(2) << (sortTime + removeTime + kdTime + verifyTime)
      << "  sortTime = " << sortTime << "  removeTime = " << removeTime
      << "  kdTime = " << kdTime << "  verifyTime = " << verifyTime << std::endl << std::endl;
#endif //TEST_KD_TREE
    // Delete all but the first of the references arrays.
    for (size_t i = 1; i < N + 1; ++i) {
      delete[] references[i];
    }
    delete[] references;

    // Return the pointer to the root of the k-d tree.
    return root;
  }

  /*
   * Walk the k-d tree to delete each KdNode.
   */
public:
  void deleteKdTree() {

    // Delete the < sub-tree.
    if (ltChild != nullptr) {
      ltChild->deleteKdTree();
    }
    // Delete the > sub-tree.
    if (gtChild != nullptr) {
      gtChild->deleteKdTree();
    }
    // Delete the current KdNode.
    delete this;
  }

  /*
   * The insideBounds function determines whether KdNode::tuple lies inside the
   * hyper-rectangle defined by the query lower and upper bound vectors.
   *
   * Calling parameters:
   *
   * queryLower - the query lower bound vector
   * queryUpper - the query upper bound vector
   * enable - a vector that specifies the dimensions on which to test for insidedness
   *
   * return true if inside, false if outside
   */
private:
  bool insideBounds(std::vector<K> const& queryLower, std::vector<K> const& queryUpper,
    std::vector<bool> const& enable) {
    bool inside = true;
    for (size_t i = 0; i < queryLower.size(); ++i) {
      if (enable[i] && (queryLower[i] > tuple[i] || queryUpper[i] < tuple[i])) {
        inside = false;
        break;
      }
    }
    return inside;
  }

  /*
   * The regionSearch function searches the k-d tree to find the KdNodes that
   * lie within a hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * queryLower - the query lower bound vector
   * queryUpper - the query upper bound vector
   * permutation - vector that specifies permutation of the partition coordinate
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * depth - the depth in the k-d tree
   * enable - a vector that specifies the dimensions on which to prune the region search
   *
   * return a list that contains the KdNodes that lie within the cutoff distance of the query node
   */
private:
  std::list<KdNode<K, V, N>*> regionSearch(std::vector<K> const& queryLower, std::vector<K> const& queryUpper,
    std::vector<size_t> const& permutation,
    signed_size_t maximumSubmitDepth, size_t depth,
    std::vector<bool> const& enable) {

    // The partition cycles as x, y, z, w...
    size_t p = permutation[depth];

    // If the KdNode is within the query hyper-rectangle for each of the k dimensions,
    // add the KdNode to the list of KdNodes that lie inside the hyper-cube. The
    // following loop is equivalent to the IN_REGION pseudo-Algol code proposed
    // by Jon Bentley in his CACM article.
    std::list<KdNode<K, V, N>*> result;
    if (insideBounds(queryLower, queryUpper, enable)) {
      result.push_back(this);
    }

    // Determine whether to search the < and > branches of the k-d tree. Although
    // the superKeyCompare function can produce a different result for the == case
    // than does comparing only the leading keys of the super-keys, that result
    // will avoid unnecessary searching of a sub-tree (at the expense of a more
    // precise super-key comparison) but the unnecessary search/ appears not to
    // change the outcome of this recursive regionSearch function.
#ifdef NO_SUPER_KEY
    bool searchLT = ltChild != nullptr && (tuple[p] >= queryLower[p] || !enable[p]);
    bool searchGT = gtChild != nullptr && (tuple[p] <= queryUpper[p] || !enable[p]);;
#else
    bool searchLT = ltChild != nullptr && (superKeyCompare(tuple, queryLower.data(), p) >= 0 || !enable[p]);
    bool searchGT = gtChild != nullptr && (superKeyCompare(tuple, queryUpper.data(), p) <= 0 || !enable[p]);
#endif

    // Do both branches require searching and is a child thread available?
    if (searchLT && searchGT && maximumSubmitDepth >= 0 && depth <= maximumSubmitDepth) {

      // Yes, both branches of the tree require searching and a child thread is available,
      // so prepare to search the < branch with a child thread.
      std::future< std::list<KdNode<K, V, N>*> > searchFuture;

      // Search the < branch?
      if (searchLT) {

        // Yes, search the < branch asynchronously with a child thread.
        searchFuture =
          std::async(std::launch::async, [&] {
          return ltChild->regionSearch(queryLower, queryUpper, permutation,
            maximumSubmitDepth, depth + 1, enable);
            });

        // Search the > branch?
        std::list<KdNode<K, V, N>*> gtResult;
        if (searchGT) {

          // Yes, search the > branch  with the master thread.
          gtResult = gtChild->regionSearch(queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
        }

        // Get the result of searching the < branch with the child thread.
        std::list<KdNode<K, V, N>*> ltResult;
        try {
          ltResult = searchFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }

        // Append the results of searching the < and > branches to the result (if any) for this KdNode.
        result.splice(result.end(), ltResult);
        result.splice(result.end(), gtResult);

      }
      else {

        // No, don't search the < branch. Search the > branch?
        std::list<KdNode<K, V, N>*> gtResult;
        if (searchGT) {

          // Yes, search the > branch  with the master thread.
          gtResult = gtChild->regionSearch(queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
        }

        // Append the result of searching the > branch to the result (if any) for this KdNode.
        result.splice(result.end(), gtResult);
      }

    }
    else {

      // No, both branches do not require searching. Search the < branch with the master thread?
      if (searchLT) {
        auto ltResult = ltChild->regionSearch(queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
        result.splice(result.end(), ltResult); // Can't substitute regionSearch(...) for ltResult.
      }

      // Search the > branch with the master thread?
      if (searchGT) {
        auto gtResult = gtChild->regionSearch(queryLower, queryUpper, permutation, maximumSubmitDepth, depth + 1, enable);
        result.splice(result.end(), gtResult); // Can't substitute regionSearch(...) for gtResult.
      }

    }

    return result;
  }

  /*
   * The searchRegion function searches the k-d tree to find the KdNodes that
   * lie within a hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * queryLower - the query lower bound vector
   * queryUpper - the query upper bound vector
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * size - the number of points in the coordinates vector after removal of duplicates
   *
   * return a list of KdNodes that lie within the query hyper-rectangle
   */
public:
  std::list<KdNode<K, V, N>*> searchRegion(std::vector<K>& queryLower, std::vector<K>& queryUpper,
    signed_size_t maximumSubmitDepth, size_t size) {

    // Determine the maximum depth of the k-d tree, which is log2(size).
    size_t maxDepth = 1;
    while (size > 0) {
      ++maxDepth;
      size >>= 1;
    }

    // It is unnecessary to compute the partition coordinate upon each recursive call
    // of the regionSearch function because that coordinate depends only on the depth
    // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
    // The partition coordinate permutes n the order 0, 1, 2, 3, 0, 1, 2, 3, etc.
    // for e.g. 4-dimensional data.
    std::vector<size_t> permutation(maxDepth);
    for (size_t i = 0; i < permutation.size(); ++i) {
      permutation.at(i) = i % queryLower.size();
    }

    // Ensure that each query lower bound <= the corresponding query upper bound.
    for (size_t i = 0; i < queryLower.size(); ++i) {
      if (queryLower[i] > queryUpper[i]) {
        K tmp = queryLower[i];
        queryLower[i] = queryUpper[i];
        queryUpper[i] = tmp;
      }
    }

    // Search the tree over all dimensions and return the resulting list of KdNodes.
    std::vector<bool> enable(queryLower.size(), true);
    return regionSearch(queryLower, queryUpper, permutation, maximumSubmitDepth, 0, enable);
  }

  /*
   * The searchRegion function searches the k-d tree to find the KdNodes that
   * lie within a hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * queryLower - the query lower bound vector
   * queryUpper - the query upper bound vector
   * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
   * size - the number of points in the coordinates vector after removal of duplicates
   * enable - a vector that specifies the dimensions on which to test for insidedness
   *          and prune the region search
   *
   * return a list of KdNodes that lie within the query hyper-rectangle
   */
public:
  std::list<KdNode<K, V, N>*> searchRegion(std::vector<K>& queryLower, std::vector<K>& queryUpper,
    signed_size_t maximumSubmitDepth, size_t size,
    std::vector<bool> const& enable) {

    // Determine the maximum depth of the k-d tree, which is log2(size).
    size_t maxDepth = 1;
    while (size > 0) {
      ++maxDepth;
      size >>= 1;
    }

    // It is unnecessary to compute the partition coordinate upon each recursive call
    // of the regionSearch function because that coordinate depends only on the depth
    // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
    // The partition coordinate permutes n the order 0, 1, 2, 3, 0, 1, 2, 3, etc.
    // for e.g. 4-dimensional data.
    std::vector<size_t> permutation(maxDepth);
    for (size_t i = 0; i < permutation.size(); ++i) {
      permutation.at(i) = i % queryLower.size();
    }

    // Ensure that each query lower bound <= the corresponding query upper bound.
    for (size_t i = 0; i < queryLower.size(); ++i) {
      if (queryLower[i] > queryUpper[i]) {
        K tmp = queryLower[i];
        queryLower[i] = queryUpper[i];
        queryUpper[i] = tmp;
      }
    }

    // Search the tree over the enabled dimensions and return the resulting list of KdNodes.
    return regionSearch(queryLower, queryUpper, permutation, maximumSubmitDepth, 0, enable);
  }

  /*
   * Walk the k-d tree and append to a list each KdNode that lies inside
   * the hyper-rectangle defined by the query lower and upper bounds.
   *
   * Calling parameters:
   *
   * queryLower - the query lower bound vector
   * queryUpper - the query upper bound vector
   *
   * return a list of KdNodes that lie within the query hyper-rectangle.
   */
public:
  std::list<KdNode<K, V, N>*> bruteRegion(std::vector<K> const& queryLower, std::vector<K> const& queryUpper) {

    // Append the KdNode to the list if it lies inside the query bounds.
    std::list<KdNode<K, V, N>*> result;
    std::vector<bool> enable(queryLower.size(), true);
    if (insideBounds(queryLower, queryUpper, enable)) {
      result.push_back(this);
    }
    // Visit the < sub-tree.
    if (ltChild != nullptr) {
      result.splice(result.end(), ltChild->bruteRegion(queryLower, queryUpper));
    }
    // Visit the > sub-tree.
    if (gtChild != nullptr) {
      result.splice(result.end(), gtChild->bruteRegion(queryLower, queryUpper));
    }

    return result;
  }

  /*
   * Search the k-d tree for all possible M nearest geometric neighbors by adding them
   * to the NearestNeighborHeap.  Exclude a branch of the tree wherein it is guaranteed
   * that all nodes in that branch are farther way than the current farthest node stored
   * in the heap.
   *
   * Calling parameters:
   *
   * heap - an instance of NearestNeighborHeap
   * permutation - vector that specifies permutation of the partition coordinate
   * depth - depth in the k-d tree
   */
private:
  void nearestNeighbors(NearestNeighborHeap<K, V, N>& heap, std::vector<size_t> const& permutation,
    size_t depth) {

    // The partition permutes as x, y, z, w...
    size_t p = permutation.at(depth);

    // If query[p] < tuple[p], descend the < branch to the bottom of the tree before adding a point to the
    // heap, which increases the probability that closer nodes to the query point will get added earlier,
    // thus reducing the likelihood of adding more distant points that get kicked out of the heap later.
    if (heap.query[p] < tuple[p]) {
      if (ltChild != nullptr) {  // If not at the bottom of the tree, descend the < branch unconditionally.
        ltChild->nearestNeighbors(heap, permutation, depth + 1);
      }
      // If the current node is closer to the query point than the farthest item in the heap, or if this
      // component of the array is not part of the nearest neighbor search, or if the heap is not full,
      // descend the > branch and then attempt to add the node to the heap.
      double dist = static_cast<double>(tuple[p] - heap.query[p]); // May result in loss of precision.
      if (dist * dist <= heap.curMaxDist() || !heap.enable[p] || !heap.heapFull()) {
        if (gtChild != nullptr) { // If not at the bottom of the tree, descend the > branch
          gtChild->nearestNeighbors(heap, permutation, depth + 1);
        }
        heap.add(this);  // Attempt to add the current KdNode to the heap.
      }
    }
    // If query[p] > tuple[p], descend the > branch to the bottom of the tree before adding a point to the
    // heap, which increases the probability that closer nodes to the query point will get added earlier,
    // thus reducing the likelihood of adding more distant points that get kicked out of the heap later.
    else if (heap.query[p] > tuple[p]) {
      if (gtChild != nullptr) {  // If not at the bottom of the tree, descend the > branch unconditionally.
        gtChild->nearestNeighbors(heap, permutation, depth + 1);
      }
      // If the current node is closer to the query point than the farthest item in the heap, or if this
      // component of the array is not part of the nearest neighbor search, or if the heap is not full,
      // descend the < branch and then attempt to add the node to the heap.
      double dist = static_cast<double>(tuple[p] - heap.query[p]); // May result in loss of precision.
      if (dist * dist <= heap.curMaxDist() || !heap.enable[p] || !heap.heapFull()) {
        if (ltChild != nullptr) {
          ltChild->nearestNeighbors(heap, permutation, depth + 1);
        }
        heap.add(this);  // Attempt to add the current node to the heap.
      }
    }
    // Because query[p] == tuple[p], the probability of finding nearest neighbors is equal for both branches
    // of the tree, so descend both branches and then attempt to add the current node to the heap.
    else {
      if (ltChild != nullptr) {
        ltChild->nearestNeighbors(heap, permutation, depth + 1);
      }
      if (gtChild != nullptr) {
        gtChild->nearestNeighbors(heap, permutation, depth + 1);
      }
      heap.add(this);  // Attempt to add the current node to the heap.
    }
  }

  /*
   * Find M nearest neighbors to the query vector and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * query - the query vector
   * numNeighbors - the number M of nearest neighbors to find
   * size - the number of points in the coordinates vector after removal of duplicates
   *
   */
public:
  std::list< std::pair<double, KdNode<K, V, N>*> > findNearestNeighbors(std::vector<K> const& query,
    size_t numNeighbors,
    size_t size) {

    // Determine the maximum depth of the k-d tree, which is log2(size).
    size_t maxDepth = 1;
    while (size > 0) {
      ++maxDepth;
      size >>= 1;
    }

    // It is unnecessary to compute the partition coordinate upon each recursive call
    // of the nearestNeighbors function because that coordinate depends only on the depth
    // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
    // The partition coordinate permutes n the order 0, 1, 2, 3, 0, 1, 2, 3, etc.
    // for e.g. 4-dimensional data.
    size_t numDimensions = query.size();
    std::vector<size_t> permutation(maxDepth);
    for (size_t i = 0; i < permutation.size(); ++i) {
      permutation.at(i) = i % numDimensions;
    }

    // Create the heap and search the k-d tree for nearest neighbors.
    NearestNeighborHeap<K, V, N> heap(query, numNeighbors);
    nearestNeighbors(heap, permutation, 0);

    // Empty the heap by successively removing the top of the heap and appending it to a list.
    // Then reverse the list so that the results are ordered by increasing distance to the query.
    std::list< std::pair<double, KdNode<K, V, N>*> > result;
    for (size_t i = 0; i < numNeighbors; ++i) {
      result.push_back(heap.removeTop());
    }
    result.reverse();
    return result;
  }

  /*
   * Find M nearest neighbors to the query vector and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * query - the query vector
   * numNeighbors - the number M of nearest neighbors to find
   * size - the number of points in the coordinates vector after removal of duplicates
   * enable - a vector that specifies the dimensions on which to test distance
   *
   */
public:
  std::list< std::pair<double, KdNode<K, V, N>*> > findNearestNeighbors(std::vector<K> const& query,
    size_t numNeighbors,
    size_t size,
    std::vector<bool> const& enable) {

    // Determine the maximum depth of the k-d tree, which is log2(size).
    size_t maxDepth = 1;
    while (size > 0) {
      ++maxDepth;
      size >>= 1;
    }

    // It is unnecessary to compute the partition coordinate upon each recursive call
    // of the nearestNeighbors function because that coordinate depends only on the depth
    // of recursion, so it may be pre-computed and stored in the 'permutation' vector.
    // The partition coordinate permutes n the order 0, 1, 2, 3, 0, 1, 2, 3, etc.
    // for e.g. 4-dimensional data.
    size_t numDimensions = query.size();
    std::vector<size_t> permutation(maxDepth);
    for (size_t i = 0; i < permutation.size(); ++i) {
      permutation.at(i) = i % numDimensions;
    }

    // Create the heap and search the k-d tree for nearest neighbors.
    NearestNeighborHeap<K, V, N> heap(query, numNeighbors, enable);
    nearestNeighbors(heap, permutation, 0);

    // Empty the heap by successively removing the top of the heap and appending it to a list.
    // Then reverse the list so that the results are ordered by increasing distance to the query.
    std::list< std::pair<double, KdNode<K, V, N>*> > result;
    for (size_t i = 0; i < numNeighbors; ++i) {
      result.push_back(heap.removeTop());
    }
    result.reverse();
    return result;
  }

  /*
   * Walk the k-d tree and attempt to add each KdNode to the NearestNeighborHeap.
   *
   * Calling parameter:
   *
   * heap - an instance of NearestNeighborHeap
   */
private:
  void allNeighbors(NearestNeighborHeap<K, V, N>& heap) {

    // Visit the < sub-tree.
    if (ltChild != nullptr) {
      ltChild->allNeighbors(heap);
    }
    // Visit the > sub-tree.
    if (gtChild != nullptr) {
      gtChild->allNeighbors(heap);
    }
    // Attempt to add the current KdNode to the heap.
    heap.add(this);
  }

  /*
   * Find M nearest neighbors to the query vector via brute force and return them as a list ordered by increasing distance.
   *
   * Calling parameters:
   *
   * query - the query vector
   * numNeighbors - the number M of nearest neighbors to find
   *
   */
public:
  std::list< std::pair<double, KdNode<K, V, N>*> > bruteNearestNeighbors(std::vector<K> const& query, size_t numNeighbors) {

    // Create the heap, walk the k-d tree, and attempt to add each KdNode to the heap.
    NearestNeighborHeap<K, V, N> heap(query, numNeighbors);
    allNeighbors(heap);

    // Empty the heap by successively removing the top of the heap and appending it to a list.
    // Then reverse the list so that the results are ordered by increasing distance to the query.
    std::list< std::pair<double, KdNode<K, V, N>*> > result;
    for (size_t i = 0; i < numNeighbors; ++i) {
      result.push_back(heap.removeTop());
    }
    result.reverse();
    return result;
  }

  /*
   * The printTuple function prints one tuple.
   *
   * Calling parameters:
   *
   * tuple - the tuple to print
   * dim - the number of dimensions
   */
public:
  static void printTuple(K const* tuple)
  {
    std::cout << "(" << tuple[0] << ",";
    for (size_t i = 1; i < N - 1; ++i) {
      std::cout << tuple[i] << ",";
      std::cout << tuple[N - 1] << ")";
    }
  }

  /*
   * The printTuple function prints one tuple.
   *
   * Calling parameter:
   *
   * tuple - the tuple as a vector
   */
public:
  static void printTuple(std::vector<K> const& tuple)
  {
    std::cout << "(" << tuple[0] << ",";
    for (size_t i = 1; i < tuple.size() - 1; ++i) {
      std::cout << tuple[i] << ",";
      std::cout << tuple[tuple.size() - 1] << ")";
    }
  }

  /*
   * The printKdTree function prints the k-d tree "sideways" with the root at the ltChild.
   *
   * Calling parameters:
   *
   * dim - the number of dimensions
   * depth - the depth in the k-d tree
   */
public:
  void printKdTree(size_t depth)
  {
    if (gtChild != nullptr) {
      gtChild->printKdTree(depth + 1);
    }
    for (size_t i = 0; i < depth; ++i) std::cout << "       ";
    printTuple(tuple);
    std::cout << std::endl;
    if (ltChild != nullptr) {
      ltChild->printKdTree(depth + 1);
    }
  }
}; // class KdNode

/*
 * The NearestNeighborHeap Class implements a std::fixed length heap of both containing both a KdNode and euclidean distance
 * from the tuple in the node to a query point.  When a KdNode is added to the heap it is unconditionally placed in
 * the heap until the heap is full.  After the heap is full, a KdNode is added to the heap only if the calculated
 * distance from the query point to the tuple is less than the farthest KdNode currently in the heap; and in that
 * case, the current farthest KdNode and distance are removed from the heap to make room for it.
 *
 * The heap is maintained in two corresponding std::fixed length vectors, one for the KdNodes and one for the distance to
 * the query point.  These vectors are stored in order of decreasing distance.  When a new node is added, regardless
 * of whether or not the heap is full, a heap sort is done to place the new KdNode in the proper order in the heap.
 * Hence, the farthest KdNode is always at the top of the heap (index 1).
 *
 * For a discussion of heap sort and a priority queue implemented via a heap, see Section 2.4 "Priority Queues"
 * pp. 308-335 in "Algorithms Fourth Edition" by Robert Sedgewick and Kevin Wayne, Addison-Wesley, 2011.
 */
template <typename K, typename V, size_t N>
class NearestNeighborHeap {
public:
  std::vector<K> query; // query point for which nearest neighbors will be found
  std::vector<bool> enable;
private:
  size_t reqDepth; // requested number of nearest neighbors
  std::vector<KdNode<K, V, N>* > nodes; // vector of KdNodes that are the nearest neighbors
  std::vector<double> dists; // vector of squared distances
  size_t curDepth; // number of nearest nodes/distances on the heap

  /*
   * Constructor that enables distance test on all dimensions
   *
   * Calling parameters:
   *
   * query - a vector that defines the query point
   * numNeighbors - the number of nearest neighbors desired
   */
public:
  NearestNeighborHeap(std::vector<K> const& query, size_t numNeighbors) {
    this->nodes.resize(numNeighbors + 1, nullptr); // heap of KdNode* (address 0 is unused)
    this->dists.resize(numNeighbors + 1, 0); // corresponding heap of distances (initialized to 0)
    this->reqDepth = numNeighbors;
    this->curDepth = 0;
    this->query = query;
    this->enable.assign(query.size(), true);
  }

  /*
   * Constructor that enables distance test on only specified dimensions
   *
   * Calling parameters:
   *
   * query - a vector that defines the query point
   * numNeighbors - the number of nearest neighbors desired
   * enable - a vector that specifies the dimensions on which to test distance
   */
public:
  NearestNeighborHeap(std::vector<K> const& query, size_t numNeighbors, std::vector<bool> const& enable) {
    this->nodes.resize(numNeighbors + 1, nullptr); // heap of KdNode* (address 0 is unused)
    this->dists.resize(numNeighbors + 1, 0); // corresponding heap of distances (initialized to 0)
    this->reqDepth = numNeighbors;
    this->curDepth = 0;
    this->query = query;
    this->enable = enable;
  }

  /*
   * Swap two elements in the heap.
   *
   * Calling parameters:
   *
   * i - the index of the first element
   * j - the index of the second element
   */
private:
  void swap(size_t i, size_t j) {
    double tempDist = dists[i];
    KdNode<K, V, N>* tempNode = nodes[i];
    dists[i] = dists[j];
    nodes[i] = nodes[j];
    dists[j] = tempDist;
    nodes[j] = tempNode;
  }

  /*
   * Allow an element to rise upward through the heap.
   *
   * Calling parameter:
   *
   * k - the index of the element
   */
private:
  void rise(size_t k) {
    while (k > 1 && dists[k / 2] < dists[k]) {
      swap(k / 2, k);
      k = k / 2;
    }
  }
  /*
   * Allow an element to fall downward through the heap.
   *
   * Calling parameter:
   *
   * k - the index of the element
   */
private:
  void fall(size_t k) {
    while (2 * k <= curDepth) {
      size_t j = 2 * k;
      if (j < curDepth && dists[j] < dists[j + 1]) {
        ++j;
      }
      if (dists[k] >= dists[j]) {
        break;
      }
      swap(k, j);
      k = j;
    }
  }

  /*
   * Remove the top element of the heap and re-order the remaining elements.
   *
   * return a pair that contains a pointer to the top KdNode and the distance to that KdNode
   */
public:
  std::pair<double, KdNode<K, V, N>*> removeTop() {
    std::pair<double, KdNode<K, V, N>*> returnPair = std::make_pair(dists[1], nodes[1]);
    swap(1, curDepth--);
    nodes[curDepth + 1] = nullptr;
    fall(1);
    return returnPair;
  }

  /*
   * Add a new KdNode to the NearestNeighborHeap if its distance to the
   * query point is less than curMaxDistance.
   *
   * Calling parameter:
   *
   * node - KdNode to potentially add to the heap
   */
public:
  void add(KdNode<K, V, N>* node) {
    // Find the distance by subtracting the query from the tuple and
    // calculating the sum of the squared distances. Note that conversion
    // from type K to double may result in loss of precision but avoids
    // the possibility of integer overflow.
    double dist2 = 0.0;
    for (size_t i = 0; i < query.size(); ++i) {
      // Add the squared coordinate distance only if the dimension is enabled.
      if (enable[i]) {
        K comp = node->tuple[i] - query[i];
        dist2 += static_cast<double>(comp) * static_cast<double>(comp);
      }
    }
    // If the queue is not full, add the point to the bottom of the heap unconditionally and let it rise;
    if (!heapFull()) {
      dists[++curDepth] = dist2;
      nodes[curDepth] = node;
      rise(curDepth);
    }
    // otherwise, if the point is closer than the top of the heap, overwrite the top and let it fall.
    else if (dist2 < curMaxDist()) {
      dists[1] = dist2;
      nodes[1] = node;
      fall(1);
    }
    return;
  }

  /* Return the current maximum distance, i.e., dists[1] */
public:
  double curMaxDist() {
    return dists[1];
  }

  /* Return true if the heap is full. */
public:
  bool heapFull() {
    return curDepth >= reqDepth;
  }
}; // class NearestNeighborHeap
