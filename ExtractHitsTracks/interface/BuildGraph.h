
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "test/ExtractHitsTracks/interface/Data.h"
#include <vector>
#include <math.h> 
#include <cmath> 
#include "test/ExtractHitsTracks/interface/OutputGraph.h"






template <class T>

std::vector<int> find_all_indices_in_vec(std::vector<T> vec, float num){ 
/* Will return the indicies of all the values num in a vector vec
 */ 


  std::vector<int> B; 
  typename std::vector<T>::iterator it = vec.begin();
  while ((it = std::find_if(it, vec.end(), [num](int x){return x == num; })) != vec.end())
  {
    B.push_back(std::distance(vec.begin(), it));
    it++;
  }
return B; 
}


template <class S>
std::vector<S> subset_vec(std::vector<S> vec, std::vector<int> indices){ 
/*returns a vector vec subset by a vector of indices 
 */

  std::vector<S> result(indices.size(), 0); 
  std::transform(indices.begin(), indices.end(), result.begin(), [vec](size_t pos) {return vec[pos];});
  return result; 
}



 
class BuildGraph {
  
  public:
  BuildGraph() {}  
  OutputGraph outgraph; 
 // std::vector<int> vec_hit_index1, vec_hit_index2; 
 
  float calc_r (float x, float y) { return sqrt(pow(x, 2)+ pow(y, 2));}
  
  //OutputGraph build_graph(Ntuple *nt){ 
   // std::vector<float> accepted_indices = select_hits(nt);  
    
  std::vector<int> select_hits(Ntuple *nt){
  /* Performs pt cut 
   * If there are multiple hits in a layer for a given particle, it selects the one with smallest sim_dxy_sig (should be changed in the long run) 
   */
  
  
    // pt cut 
    
    std::vector<int> accepted_hit_indices;
    std::vector<int> new_accepted_hit_indices;  
    //should be a config variable instead
    int pt_cut = 2;  
  
    auto it = std::find_if(std::begin(nt->sim_pt_), std::end(nt->sim_pt_), [pt_cut](int i){return i > pt_cut;});
    while (it != std::end(nt->sim_pt_)) {
      accepted_hit_indices.emplace_back(std::distance(std::begin(nt->sim_pt_), it));
      it = std::find_if(std::next(it), std::end(nt->sim_pt_), [pt_cut](int i){return i > pt_cut;});
    }
  
  
    // remove duplicate hits in a layer for a given particle 
  
    std::vector<int> pt_cut_particle_id = subset_vec(nt->particle_id_,accepted_hit_indices); 
    // find unqiue particle ids 
    std::unordered_set<int> s(pt_cut_particle_id.begin(), pt_cut_particle_id.end());
    std::vector<int> unique_particle_ids; 
    unique_particle_ids.assign(s.begin(), s.end());  
    for(auto particle: pt_cut_particle_id){
      // find the positions in the list of that particle 
          
      //contains the indices of the partilce in the particle cut list 
      std::vector<int> positions = find_all_indices_in_vec(pt_cut_particle_id, particle); 
      
      //layer ids for a given particle 
      auto particle_layer_ids = subset_vec(nt->layer_id_,positions); 
      std::unordered_set<int> s(particle_layer_ids.begin(), particle_layer_ids.end());
      std::vector<int> unique_particle_layer_ids; 
      unique_particle_layer_ids.assign(s.begin(), s.end());  
      //check if the layer id is repeated 
      std::map<int, int> CountMap;
      for (auto layer = unique_particle_layer_ids.begin(); layer!= unique_particle_layer_ids.end(); layer++)
  	CountMap[*layer]++; 
      
      // for each layer for the particle 
      for (auto count = CountMap.begin(); count!=CountMap.end(); count++){
  	// if the layer id is repeated, find the sim_dxy_sig and return the index of the smallest values for this 
  	if (count->second > 1){
  	  int repeated_layer = count->first; 
  	  std::vector<int> indices_repeated_layers = find_all_indices_in_vec(particle_layer_ids, repeated_layer); 
  	  // the position in the data array will be the index for the particle plus the index for the layer 
            std::for_each(indices_repeated_layers.begin(), indices_repeated_layers.end(), [positions](int &d) { d+=positions.front();});
            
  	  std::vector<float> sim_dxy_sig_multiple_layer_hit = subset_vec(nt->sim_dxy_sig_,indices_repeated_layers); 
  	  //std::vector<int>::iterator smallest_sim_dxy_sig = std::min_element(sim_dxy_sig_for_multiple_layer_hit.begin(), sim_dxy_sig_for_multiple_layer_hit.end());
  	  float smallest_sim_dxy_sig = *std::min_element(sim_dxy_sig_multiple_layer_hit.begin(), sim_dxy_sig_multiple_layer_hit.end());
      	  auto it = std::find(sim_dxy_sig_multiple_layer_hit.begin(), sim_dxy_sig_multiple_layer_hit.end(), smallest_sim_dxy_sig);
  	  int index_min_sim_dxy_sig = std::distance(std::begin(sim_dxy_sig_multiple_layer_hit), it); 
  	  new_accepted_hit_indices.emplace_back(positions.front() + index_min_sim_dxy_sig); 	    	    
   	}     
     
      }//layer id counts for a given particle
   
    }// particle id 
     
   // populate the output graph 
   std::vector<float> r; 
   for (size_t i = 0; i < nt->x_.size(); i++){
     r.push_back(calc_r(nt->x_.at(i), nt->y_.at(i)));   
   }
   // scale r and fill to graph  
   std::transform(r.begin(), r.end(), r.begin(), [](float &c){ return c/1000; });
   outgraph.X.r = subset_vec(r, new_accepted_hit_indices); 
   // scale phi and fill 
   std::vector<float> sim_phi = nt->sim_phi_; 
   std::transform(sim_phi.begin(), sim_phi.end(), sim_phi.begin(), [](float &c){ return c/3.14; });
   outgraph.X.phi = subset_vec(sim_phi, new_accepted_hit_indices); 


   std::vector<float> z = nt->z_; 
   std::transform(z.begin(), z.end(), z.begin(), [](float &c){ return c/1000; });
   outgraph.X.z = subset_vec(z, new_accepted_hit_indices); 

   //particle level - for inference this will be removed! 
 //  outgraph.particle_id = subset_vec(nt->particle_id_, new_accepted_hit_indices); 
  // outgraph.sim_pt = subset_vec(nt->sim_pt_, new_accepted_hit_indices); 
   //outgraph.sim_eta = subset_vec(nt->sim_eta_, new_accepted_hit_indices); 

   return new_accepted_hit_indices; 
  }
  
 void  create_hit_pairs(Ntuple *nt, std::vector<int> accepted_indices){ 
    // create the layer combos 
    std::vector<std::pair<int, int>> layer_combos; 
    for(int i =1; i!=28; i++){
    // all layers can connect to the next except one endcap to another   
      if(i!= 16){layer_combos.push_back(std::pair(i, i+1));}
    }
    // also allow any of the inner barrel layers to connect to first layer in either endcap 
    for(int i = 1; i!=5; i++){
      layer_combos.push_back(std::pair(i, 5)); 
      layer_combos.push_back(std::pair(i, 17)); 
    }
    // create two list of indices containing all the hits in the two compatible layers 
    for(auto layer_combo: layer_combos){
      std::vector<int> hit1_indices = find_all_indices_in_vec(nt->layer_id_, layer_combo.first);     
      std::vector<int> hit2_indices = find_all_indices_in_vec(nt->layer_id_, layer_combo.second);     
      select_segments(nt, hit1_indices, hit2_indices); 
    }

    //return std::pair(hit1_indices, hit2_indices); 

  }

  void select_segments(Ntuple *nt, std::vector<int> hit1_indices, std::vector<int> hit2_indices){ 
    // find all possible combinations of the two lists 
    std::vector<std::pair<int, int>> index_accepted_segment; 
    for(auto hit1_index: hit1_indices){ 
      for(auto hit2_index: hit2_indices){ 
        float dphi = calc_dphi(nt->sim_phi_.at(hit1_index), nt->sim_phi_.at(hit2_index));
	float dz = nt->z_.at(hit2_index) - nt->z_.at(hit1_index); 
	float r1 = sqrt(pow(nt->x_.at(hit1_index),2) + pow(nt->y_.at(hit1_index),2)); 
	float r2 = sqrt(pow(nt->x_.at(hit2_index),2) + pow(nt->y_.at(hit2_index),2)); 
	float dr = r2 - r1; 
	float eta1 = calc_eta(r1, nt->z_.at(hit1_index)); 
	float eta2 = calc_eta(r2, nt->z_.at(hit2_index)); 
	float deta = eta2 - eta1; 
	float dR = sqrt(pow(deta, 2) + pow(dphi,2)); 
	float phi_slope = dphi/dr; 
	float z0 = nt->z_.at(hit1_index) - r1*dz/dr; 


	// these parameters should go in the python config file 
	float phi_slope_max = 0.06; 
	float z0_max = 27; 
        if ((z0 < z0_max) & (phi_slope < phi_slope_max)) { 
	  index_accepted_segment.push_back(std::pair(hit1_index, hit2_index)); 
   
	  outgraph.edge_attr.dr.push_back(dr/1000); 
          outgraph.edge_attr.dphi.push_back(dphi/3.14); 
	  outgraph.edge_attr.dz.push_back(dz/1000); 
          outgraph.edge_attr.dR.push_back(dR); 
 	  outgraph.edge_index.seg_start.push_back(hit1_index); 
	  outgraph.edge_index.seg_end.push_back(hit2_index);  
	
	
        } 
      }
    }
}		


  // probably root functions for this somewhere
  float calc_dphi(float phi1, float phi2){ 
    float dphi = phi2 - phi1; 
    if (dphi > 3.14){ 
	dphi -= 2*3.14; 
     }  

    if (dphi <  -3.14){ 
	dphi += 2*3.14; 
     }
   return dphi; 
  }  

  //probably also a root function 
  float calc_eta(float r, float z){
    float theta = atan2(r, z); 
  return -log(atan(theta/2)); 
  } 
   
 
//  void build_graph(*Ntuple nt, std::vector<std::pair<int,int>> accepted_segment_indices){
    // 
    // feature scaling  

  //}
  };




