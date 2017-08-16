/*---
function [E_rowIds, F_rowIds, rowCorIds] = unproject(...
 * edges_S, edges_T, cost_ST, [nS, nT])
Input: 

Output

--*/

#include "mex.h"
#include "LPCons.h"
#include <vector>
#include <algorithm>
using namespace std;

// The three input matrices are
// 1) The mesh vertices
// 2) The mesh faces
// 3) The camera parameters

// The output parameters
// 1) The intersecting pixels

void mexFunction(
     int nargout,
     mxArray *output[],
     int nargin,
     const mxArray *input[]) {
  /* check argument */
  if (nargin != 4) {
    mexErrMsgTxt("Four input arguments required.");
  }
  if (nargout != 3) {
    mexErrMsgTxt("Three output output arguments."); 
  }

  double *data_edges_S = (double*)mxGetData(input[0]);
  unsigned numEdges_S = static_cast<unsigned> (mxGetN(input[0]));
  double *data_edges_T = (double*)mxGetData(input[1]);
  unsigned numEdges_T = static_cast<unsigned> (mxGetN(input[1]));
  double *data_cost = (double*)mxGetData(input[2]);
  unsigned numCorres = static_cast<unsigned> (mxGetN(input[2]));
  double *data = (double*)mxGetData(input[3]);
  unsigned numV_S = static_cast<unsigned> (data[0]);
  unsigned numV_T = static_cast<unsigned> (data[1]);
  
  LPCons lp_cons;
  
  lp_cons.G_s.resize(numV_S);
  for (unsigned eId = 0; eId < numEdges_S; ++eId) {
      unsigned v1Id = static_cast<unsigned> (data_edges_S[2*eId]-0.5);
      unsigned v2Id = static_cast<unsigned> (data_edges_S[2*eId+1]-0.5);
      if (v1Id < v2Id) {
          lp_cons.G_s[v1Id].push_back(v2Id);
          lp_cons.G_s[v2Id].push_back(v1Id);
      }
  }
  lp_cons.G_t.resize(numV_T);
  for (unsigned eId = 0; eId < numEdges_T; ++eId) {
      unsigned v1Id = static_cast<unsigned> (data_edges_T[2*eId]-0.5);
      unsigned v2Id = static_cast<unsigned> (data_edges_T[2*eId+1]-0.5);
      if (v1Id < v2Id) {
          lp_cons.G_t[v1Id].push_back(v2Id);
          lp_cons.G_t[v2Id].push_back(v1Id);
      }
  }
  
  lp_cons.Cost.resize(numV_S);
  for (unsigned cId = 0; cId < numCorres; ++cId) {
      unsigned sId = static_cast<unsigned> (data_cost[2*cId]-0.5);
      unsigned tId = static_cast<unsigned> (data_cost[2*cId+1]-0.5);
      lp_cons.Cost[sId].push_back(tId);
  }
  
  vector<unsigned> E_rowIds, F_rowIds, rowCorIds;
  lp_cons.Compute(&E_rowIds, &F_rowIds, &rowCorIds);
  
  output[0] = mxCreateDoubleMatrix(1, E_rowIds.size(), mxREAL);
  data = mxGetPr(output[0]);
  for (unsigned i = 0; i < E_rowIds.size(); ++i)
      data[i] = E_rowIds[i]+1;
  
  output[1] = mxCreateDoubleMatrix(1, F_rowIds.size(), mxREAL);
  data = mxGetPr(output[1]);
  for (unsigned i = 0; i < F_rowIds.size(); ++i)
      data[i] = F_rowIds[i]+1;
  
  output[2] = mxCreateDoubleMatrix(1, rowCorIds.size(), mxREAL);
  data = mxGetPr(output[2]);
  for (unsigned i = 0; i < rowCorIds.size(); ++i)
      data[i] = rowCorIds[i]+1;
}


