#ifndef LPCons_h
#define LPCons_h

#include <vector>
using namespace std;

class LPCons {
 public:
  LPCons() {
  }
  ~LPCons() {
  }

  void Compute(vector<unsigned> *E1_corIds,
	  vector<unsigned> *E2_corIds,
	  vector<unsigned> *rowCorIds);


  vector<vector<unsigned>> G_s;
  vector<vector<unsigned>> G_t;
  vector<vector<unsigned>> Cost;

  unsigned numV_s;
  unsigned numV_t;
};

#endif
