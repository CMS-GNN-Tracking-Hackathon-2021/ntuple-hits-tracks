#include<fstream> 
#include<iostream> 
#include "interface/OutGraph.h"

int main(){ 
  std::ifstream file; 
  file.open("built_graph.csv", std::ios::in); 
  OutGraph obj; 
  file.read((char*)&obj, sizeof(obj)); 
  file.close(); 
}  
