#include "hierarchicalMaterials.h"

hexagonal_unit_cell::hexagonal_unit_cell(int i, int j, int a, set<string> * all_cells){
  this->i = i;
  this->j = j;
  this->a = a;
  ostringstream stringStream;
  stringStream << i << "_" << j;
  this->id = stringStream.str();

  //if(all_cells->find(this->id) != all_cells->end()){
  //  ostringstream stringStream;
  // stringStream << "Unit cell at (" << i << ", " << j << ") already exists";
  //  throw invalid_argument(stringStream.str());
  //}
  //else{
    all_cells->insert(this->id);
  //}

  this->nodes.push_back(new node(2*i*a, 2*j*a));
  this->nodes.push_back(new node((2*i+1)*a, (2*j+1)*a));
  this->interior_bonds.push_back(new bond(this->nodes[0], this->nodes[1], "I"));
}

bool hexagonal_unit_cell::has_neighbor(hexagonal_unit_cell *neigh){
  for(auto const& nb: this->neighbors) {
    if (nb->id == neigh->id){
      return true;
    };
  }
  return false;
}

void hexagonal_unit_cell::add_neighbor(hexagonal_unit_cell *neigh, string neigh_type){
  /*if(this->has_neighbor(neigh) && neigh->has_neighbor(this)){
    ostringstream stringStream;
    stringStream << "[huc::add_neighbor] Unit cells " << this->id << ", " << neigh->id << " are already neighbors";
    throw invalid_argument(stringStream.str());
  }
  if(this->has_neighbor(neigh) || neigh->has_neighbor(this)){
    ostringstream stringStream;
    stringStream << "[huc::add_neighbor] Unit cells " << this->id << ", " << neigh->id << " are partial neighbors";
    throw invalid_argument(stringStream.str());
  }*/
  vector <int> bond_n1, bond_n2;
  if (neigh_type == "PR" || neigh_type == "R"){
    bond_n1.push_back(0);
    bond_n1.push_back(1);
    bond_n1.push_back(1);
    bond_n2.push_back(0);
    bond_n2.push_back(1);
    bond_n2.push_back(0);
  }
  else if (neigh_type == "PU" || neigh_type == "U" || neigh_type == "UR" || neigh_type == "PUR" || neigh_type == "C"){
    bond_n1.push_back(1);
    bond_n2.push_back(0);
  }
  else{
    throw invalid_argument("[huc::add_neighbor] Unknown neighboring cell type");
  }
  int i = 0;
  for (i = 0; i < bond_n1.size(); i++){
    bond *b = new bond(this->nodes[bond_n1[i]], neigh->nodes[bond_n2[i]], neigh_type);
    this->exterior_bonds.push_back(b);
    neigh->exterior_bonds.push_back(b);
    this->neighbors.push_back(neigh);
    neigh->neighbors.push_back(this);
  }

}