#include "stdafx.h"
#include "LPCons.h"
#include <algorithm>

// Perform binary search to obtain the index
int GetIndex(const vector<unsigned> &array, const unsigned &neighborId) {
  unsigned len = array.size();
  if (len == 0)
    return -1;
  if (neighborId < array[0])
	return -1;
  if (neighborId > array[len-1])
	return -1;

  unsigned left = 0;
  unsigned right = len-1;
  while (left+1 < right) {
    unsigned mid = static_cast<unsigned> ((left+right)/2);
	if (array[mid] == neighborId)
	  return mid;
	if (array[mid] > neighborId)
	  right = mid;
	if (array[mid] < neighborId)
	  left = mid;
  }
  if (array[left] == neighborId)
	return left;
  if (array[right] == neighborId)
	return right;
  return -1;
}

void LPCons::Compute(
  vector<unsigned> *E1_corIds,
  vector<unsigned> *E2_corIds,
  vector<unsigned> *rowCorIds) {
  // the graphs have to be symmetric
 unsigned numV_s = G_s.size();
 for (unsigned sId = 0; sId < numV_s; ++sId)
   sort(G_s[sId].begin(), G_s[sId].end());

 unsigned numV_t = G_t.size();
 for (unsigned tId = 0; tId < numV_t; ++tId)
   sort(G_t[tId].begin(), G_t[tId].end());

 // Count the sparse structure of the correspondence matrix
 vector<unsigned> offs;
 offs.resize(numV_s+1);
 offs[0] = 0;
 for (unsigned i = 0; i < numV_s; ++i)
   offs[i+1] = offs[i] + Cost[i].size();
 unsigned numCorres  = offs[numV_s];

  // number of neighbors
 vector<unsigned> valance;
 valance.resize(G_s.size());
 for (unsigned vId = 0; vId < G_s.size(); ++vId)
   valance[vId] = G_s[vId].size();

 vector<unsigned> offs_varX;
 offs_varX.resize(numCorres+1);
 offs_varX[0] = 0;
 for (unsigned vId = 0; vId < numV_s; ++vId) {
   for (unsigned corId = offs[vId]; corId < offs[vId+1]; ++corId) {
     offs_varX[corId+1] = offs_varX[corId] + valance[vId];
   }
 }

 rowCorIds->resize(offs_varX[numCorres]);
 unsigned off = 0;
 for (unsigned vId = 0; vId < numV_s; ++vId) {
   for (unsigned corId = offs[vId]; corId < offs[vId+1]; ++corId) {
     unsigned vId_t = Cost[vId][corId-offs[vId]];
     unsigned corId_dense = vId_t*numV_s + vId;
     for (unsigned i = 0; i < G_s[vId].size(); ++i) {
       (*rowCorIds)[off] = corId_dense;
       off++;
     }
   }
 }


 unsigned numAlloc = 8192;
 E1_corIds->resize(numAlloc);
 E2_corIds->resize(numAlloc);

 unsigned num = 0;
 for (unsigned vId = 0; vId < numV_s; ++vId) {
   for (unsigned corId = offs[vId]; corId < offs[vId+1]; ++corId) {
     unsigned vId_t = Cost[vId][corId-offs[vId]];
     for (unsigned i = 0; i < G_s[vId].size(); ++i) {
       unsigned rowId1 = offs_varX[corId] + i;

       unsigned nId = G_s[vId][i];
       if (vId >= nId)
           continue;
       
       int j = GetIndex(G_s[nId], vId);

       for (unsigned corId2 = offs[nId]; corId2 < offs[nId+1]; ++corId2) {
         unsigned nId_t = Cost[nId][corId2-offs[nId]];
         if (GetIndex(G_t[vId_t], nId_t) >= 0) {
           unsigned rowId2 = offs_varX[corId2] + j;
           if (num >= numAlloc) {
             numAlloc = numAlloc*2;
             E1_corIds->resize(numAlloc);
             E2_corIds->resize(numAlloc);
           }
           (*E1_corIds)[num] = rowId1;
           (*E2_corIds)[num] = rowId2;
           num++;
         }
       }
     }
   }
 }
 E1_corIds->resize(num);
 E2_corIds->resize(num);
}