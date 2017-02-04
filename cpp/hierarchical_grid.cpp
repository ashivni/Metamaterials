#include "hierarchicalMaterials.h"

hierarchical_grid::hierarchical_grid(int nx, int ny, int levels, int l0, int magnification, bool notch, bool damage, double damage_fraction){
  this->nx = nx;
  this->ny = ny;
  this->levels = levels;
  this->l0 = l0;
  this->magnification = magnification;
  this->notch = notch;
  this->damage = damage;
  this->damage_fraction = damage_fraction;
  this->n_step = 0;
  this->a = l0 * (pow(magnification,levels));
  this->ly = this->a * (ny * 2 - 1) * (pow(3.0,0.5)) * 0.5;
  this->lx = this->a * (nx * 2 - 1) * 0.5;
  this->ly_ind = this->a * (ny * 2 - 1);
  this->lx_ind = this->a * (nx * 2 - 1);
  this->outline = new hexagonal_grid(nx, ny, levels, magnification, l0);
  this->level_nodes = new map<int, map<string, node *> *>;
  this->level_bonds = new map<int, map<string, bond *> *>;
  this->level_broken_bonds = new map<int, map<string, bond *> *>;
  this->level_is_build = new map<int, bool>;

  (*(this->level_is_build))[this->levels] = false;
  (*(this->level_nodes))[this->levels] = this->outline->node_map();
  (*(this->level_bonds))[this->levels] = this->outline->bond_map();
  (*(this->level_broken_bonds))[this->levels] = new map<string, bond *>;

  int l = this->levels - 1;

  while(l >= 0){
    map<string, node *> * nodes = new map<string, node *>;
    map<string, bond *> * bonds = new map<string, bond *>;

    for(auto const& level_bond : *((*(this->level_bonds))[l+1])){
      nodes_and_bonds nb = level_bond.second->refine(this->magnification);
      for(auto const& refined_node_pts : nb.nodes){
        if(nodes->find(refined_node_pts.first) == nodes->end()){
          (*nodes)[refined_node_pts.first] = refined_node_pts.second;
        }
      }
      int count = 0;
      for(auto const& refined_bond_pts : nb.bonds){
        if(bonds->find(refined_bond_pts.first) == bonds->end()){
          (*bonds)[refined_bond_pts.first] = refined_bond_pts.second;
        }
        node *n1 = (*nodes)[refined_bond_pts.second->n1->id];
        node *n2 = (*nodes)[refined_bond_pts.second->n2->id];
        if(!n1->bonded_to(n2))  n1->neighbors->insert(std::make_pair(n2->get_id(),n2));
        if(!n2->bonded_to(n1))  n2->neighbors->insert(std::make_pair(n1->get_id(),n1));
      }
    }
    (*(this->level_nodes))[l] = nodes;
    (*(this->level_bonds))[l] = bonds;

    for(auto const& b: *((*(this->level_bonds))[l])){
      (b.second)->n1 = (*((*(this->level_nodes))[l]))[(b.second)->n1->id];
      (b.second)->n2 = (*((*(this->level_nodes))[l]))[(b.second)->n2->id];
    }

    set<string> duplicate_bonds;
    for(auto const& b: *((*(this->level_bonds))[l])){
      set<std::string>::iterator it = duplicate_bonds.find((b.second)->id);
      if (it == duplicate_bonds.end()){
        string dup_id = (b.second)->n2->id + (b.second)->n1->id;
        if ((*((*(this->level_bonds))[l])).find(dup_id) !=  (*((*(this->level_bonds))[l])).end()){
          duplicate_bonds.insert(dup_id);
        }
      }
    }

    for(auto const& db: duplicate_bonds){
      (*((*(this->level_bonds))[l])).erase(db);
    }
    l-=1;
  }
};
