
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "test/ExtractHitsTracks/interface/Data.h"
#include <vector>

//class TTree;

class BuildGraph {

 public:
 BuildGraph() {}  

 void select_hits(std::vector<float> *sim_pt){ 
/* Returns the positions for which the hits have accepted pt values 
 * In this function also need something to remove multiple hits in a layer for any given particle 
 *
 */
      //apply pt cut 
      std::vector<size_t> accepted_hit_indices;

      //should be a config variable instead
      int pt_cut = 2;  
      std::vector<size_t> y((*sim_pt).size());

      std::iota(y.begin(), y.end(), 0);
      std::copy_if(y.begin(), y.end(), 
      std::ostream_iterator<size_t>(std::cout, " "), 
      [&](size_t i) { return (*sim_pt)[i] > pt_cut; });
      
      //std::cout<<"The first index after applying pt cut is"<<accepted_hit_indices[0]<<std::endl; 
     
  }

 
};
