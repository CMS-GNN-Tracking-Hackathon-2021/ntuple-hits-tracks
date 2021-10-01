struct Xvec { 
  std::vector<float> r, phi, z; 
}; 

struct EdgeAttributes{ 
  std::vector<float> dr, dphi, dz, dR; 
}; 


struct EdgeIndex{ 
  std::vector<int> seg_start, seg_end; 
}; 

class OutputGraph{ 
  public: 
    Xvec X; 
    EdgeAttributes edge_attr; 
    EdgeIndex edge_index;  
   // std::vector<float> sim_pt, sim_eta; 
    //std::vector<int> particle_id; 
    

}; 
