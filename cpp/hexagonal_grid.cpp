#include "hierarchicalMaterials.h"

hexagonal_grid::hexagonal_grid(int nx, int ny, int level, int magnification, int l0){
  this->nx = nx;
  this->ny = ny;
  this->level = level;
  this->magnification = magnification;
  this->l0 = l0;
  this->a = l0*(pow(magnification,level));
  this->all_cells = new set<string>;
  int i, j;
  cout << "Hexagonal grid: " << nx << ", " << ny << ", " << level << ", " << this->l0 << ", " << this->magnification <<", " << this->a << endl;
  for(i = 0; i<nx; i++){
    for(j=0; j<ny; j++){
      hexagonal_unit_cell *hc = new hexagonal_unit_cell(i, j, this->a, this->all_cells);
      this->master_grid[(nx*j + i)*this->a] = hc;
    }
  }

  for(i = 0; i<nx; i++){
    for(j=0; j<ny; j++){
      int cell_id = (i + j*this->nx)*this->a;
      vector <int> offset_i = {1,0,1};
      vector <int> offset_j = {0,1,1};
      vector <string> n_type = {"R", "U", "UR"};
      for (int k = 0; k < 3; k++){
        int o_i = offset_i[k], o_j = offset_j[k];
        string neigh_type = n_type[k];
        if (o_i + i == this->nx && o_j + j == this->ny)
          neigh_type = "C";
        else if (o_i + i == this->nx && o_j == 1)
          neigh_type = "PUR";
        else if (o_i + i == this->nx)
          neigh_type = "PR";
        else if (o_j + j == this->ny)
          neigh_type = "PU";
        int n_i = (o_i + i) % this->nx;
        int n_j = (o_j + j) % this->ny;
        int neigh_id = (n_i + n_j*this->nx)*this->a;
        if (neigh_type == "R" || neigh_type == "U" || neigh_type == "UR"){
          this->master_grid[cell_id]->add_neighbor(this->master_grid[neigh_id], neigh_type);
        }
      }
    }
  }
}

map<string, node *> * hexagonal_grid::node_map(){
  map<string, node *> * m = new map<string, node *> ();
  for(auto const& huc : this->master_grid){
    for(auto const& n :huc.second->nodes){
      if(m->find(n->id) != m->end()){
        throw invalid_argument("[hexagonal_grid::node_map] Duplicate Node");
      }
      (*m)[n->id] = n;
    }
  }
  return m;
}

map<string, bond *> * hexagonal_grid::bond_map(){
  map<string, bond *> * m = new map<string, bond *> ();
  for(auto const& huc : this->master_grid){
    for(auto const& b :huc.second->interior_bonds){
      (*m)[b->id] = b;
    }
    for(auto const& b :huc.second->exterior_bonds){
      (*m)[b->id] = b;
    }
  }
  return m;
}

coordinates_2d hexagonal_grid::node_coordinates(){
  coordinates_2d c;
  for(auto const& huc : this->master_grid){
    for(auto const& n :huc.second->nodes){
      c.append(n->coordinates());
    }
  }
  return c;
}