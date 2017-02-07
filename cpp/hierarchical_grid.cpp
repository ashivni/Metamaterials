#include "hierarchicalMaterials.h"

hierarchical_grid::hierarchical_grid(int nx, int ny, int levels, int l0, int magnification, bool notch, bool damage, double damage_fraction, double EY){
  this->nx = nx;
  this->ny = ny;
  this->levels = levels;
  this->l0 = l0;
  this->magnification = magnification;
  this->notch = notch;
  this->damage = damage;
  this->damage_fraction = damage_fraction;
  this->EY = EY;
  this->n_step = 0;
  this->a = l0 * (pow(magnification,levels));
  this->ly = this->a * (ny * 2 - 1) * (pow(3.0,0.5)) * 0.5;
  this->lx = this->a * (nx * 2 - 1) * 0.5;
  this->ly_ind = this->a * (ny * 2 - 1);
  this->lx_ind = this->a * (nx * 2 - 1);
  this->outline = new hexagonal_grid(nx, ny, levels, magnification, l0);
  this->level_nodes = new map<int, my_node_map *>;
  this->level_bonds = new map<int, map<string, bond *> *>;
  this->level_broken_bonds = new map<int, map<string, bond *> *>;
  this->level_is_build = new map<int, bool>;

  (*(this->level_is_build))[this->levels] = false;
  (*(this->level_nodes))[this->levels] = this->outline->node_map();
  (*(this->level_bonds))[this->levels] = this->outline->bond_map();
  (*(this->level_broken_bonds))[this->levels] = new map<string, bond *>;

  int l = this->levels - 1;

  while(l >= 0){
    my_node_map * nodes = new my_node_map();
    map<string, bond *> * bonds = new map<string, bond *>;

    for(auto const& level_bond : *((*(this->level_bonds))[l+1])){
      nodes_and_bonds nb = level_bond.second->refine(this->magnification);
      for(auto const& refined_node_pts : nb.nodes){
        if(!nodes->has(refined_node_pts.first)){
          (*nodes).insert(refined_node_pts.first, refined_node_pts.second);
        }
      }
      int count = 0;
      for(auto const& refined_bond_pts : nb.bonds){
        //if(bonds->find(refined_bond_pts.first) == bonds->end()){
          (*bonds)[refined_bond_pts.first] = refined_bond_pts.second;
        //}
        node *n1 = (*nodes).value(refined_bond_pts.second->n1->id);
        node *n2 = (*nodes).value(refined_bond_pts.second->n2->id);
        //if(!n1->bonded_to(n2))  n1->neighbors->insert(std::make_pair(n2->get_id(),n2));
        //if(!n2->bonded_to(n1))  n2->neighbors->insert(std::make_pair(n1->get_id(),n1));
        n1->neighbors->insert(n2->get_id(), n2);
        n2->neighbors->insert(n1->get_id(), n1);
      }
    }
    (*(this->level_nodes))[l] = nodes;
    (*(this->level_bonds))[l] = bonds;

    for(auto const& b: *((*(this->level_bonds))[l])){
      (b.second)->n1 = (*((*(this->level_nodes))[l])).value((b.second)->n1->id);
      (b.second)->n2 = (*((*(this->level_nodes))[l])).value((b.second)->n2->id);
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


void hierarchical_grid::build_eqns(){
  int a = this->a;
  int ly = this->ly;
  int l = this->levels;
  double EY = this->EY;
  while(l >=0){
    if ( !(*(this->level_is_build))[l] ){
      int i_min = numeric_limits<int>::max(), j_min = numeric_limits<int>::max();
      int i_max = numeric_limits<int>::min(), j_max = numeric_limits<int>::min();
      int counter = 0;
      vector <string> * interior_node_ids = new vector <string>;
      vector <string> * exterior_up_node_ids = new vector <string>;
      vector <string> * exterior_dn_node_ids = new vector <string>;

      for(auto const & nm : *(*(this->level_nodes))[l]){
        counter += 1;
        int x = nm.second->i, y = nm.second->j;
        i_min = i_min < x ? i_min : x;
        i_max = i_max > x ? i_max : x;
        j_min = j_min < y ? j_min : y;
        j_max = j_max > y ? j_max : y;

        if(y < this->ly_ind && y > 0){
          interior_node_ids->push_back(nm.first);
        }
        else if(y >= this->ly_ind){
          exterior_up_node_ids->push_back(nm.first);
        }
        else{
          exterior_dn_node_ids->push_back(nm.first);
        }
      }
      sort(interior_node_ids->begin(), interior_node_ids->end());
      sort(exterior_up_node_ids->begin(), exterior_up_node_ids->end());
      sort(exterior_dn_node_ids->begin(), exterior_dn_node_ids->end());

      int N_interior = interior_node_ids->size();
      int N_exterior_up = exterior_up_node_ids->size();
      int N_exterior_dn = exterior_dn_node_ids->size();

      map <string, int> *interior_node_indices = new map <string, int>;
      map <string, int> *exterior_up_node_indices = new map <string, int>;
      map <string, int> *exterior_dn_node_indices = new map <string, int>;

      for(int k = 0; k < N_interior; k++){
        counter += 1;
        (*interior_node_indices)[(*interior_node_ids)[k]] = k;
      }

      for(int k = 0; k < N_exterior_up; k++){
        counter += 1;
        (*exterior_up_node_indices)[(*exterior_up_node_ids)[k]] = k;
      }

      for(int k = 0; k < N_exterior_dn; k++){
        counter += 1;
        (*exterior_dn_node_indices)[(*exterior_dn_node_ids)[k]] = k;
      }
      vector <int> * Li = new vector <int>, * Lj = new vector <int>;
      vector <double> * L_data = new vector<double>;
      vector <int> * Bi = new vector <int>, *Bj = new vector <int>;
      vector <double> * B_data = new vector<double>;
      vector <int> * Vi = new vector <int>, * Vj = new vector <int>;
      vector <double> * V_data = new vector<double>;
      vector <int> * Ci = new vector <int>, * Cj = new vector <int>;
      vector <double> *C_data = new vector<double>;

      for(int k = 0; k < N_interior; k++){
        counter += 1;
        int node_index = k;
        string node_id = (*interior_node_ids)[k];
        int node_x = (*((*(this->level_nodes))[l]))[node_id]->i;
        int node_y = (*((*(this->level_nodes))[l]))[node_id]->j;
        my_map <string, node *> *neighs = (*((*(this->level_nodes))[l]))[node_id]->neighbors;
        int diag = 0;
        double load = 0;

        for(auto const & neigh_pair: *neighs){
          counter += 1;
          diag += 1;
          if(find(interior_node_ids->begin(), interior_node_ids->end(), neigh_pair.first) != interior_node_ids->end()){
            int neigh_ind = (*interior_node_indices)[neigh_pair.first];
            Li->push_back(k);
            Lj->push_back(neigh_ind);
            L_data->push_back(-1.0);
          }
          else if(find(exterior_up_node_ids->begin(), exterior_up_node_ids->end(), neigh_pair.first) != exterior_up_node_ids->end()){
            load += ly_ind * EY;
          }
          else if(find(exterior_dn_node_ids->begin(), exterior_dn_node_ids->end(), neigh_pair.first) != exterior_dn_node_ids->end()){
            load += 0;
          }
          else{
            throw invalid_argument("Unknown node id");
          }

          int neigh_index;
          if(interior_node_indices->find(neigh_pair.first) != interior_node_indices->end()){
            neigh_index = (*interior_node_indices)[neigh_pair.first];
          }
          else if(exterior_dn_node_indices->find(neigh_pair.first) != exterior_dn_node_indices->end()){
            neigh_index = (*exterior_dn_node_indices)[neigh_pair.first] + N_interior;
          }
          else if(exterior_up_node_indices->find(neigh_pair.first) != exterior_up_node_indices->end()){
            neigh_index = (*exterior_up_node_indices)[neigh_pair.first] + N_interior + N_exterior_dn;
          }

          Ci->push_back(node_index);
          Cj->push_back(neigh_index);
          C_data->push_back(1.0);
        }

        if(abs(load)>1E-3){
          Bi->push_back(k);
          Bj->push_back(0);
          B_data->push_back(load);
        }

        Li->push_back(k);
        Lj->push_back(k);
        L_data->push_back(diag);

        Vi->push_back(node_index);
        Vj->push_back(node_index);
        V_data->push_back(-1.0);
      }

      for(int k = 0; k<N_exterior_up; k++){
        counter += 1;
        int node_index = k + N_interior + N_exterior_dn;
        string node_id = (*exterior_up_node_ids)[k];
        int y = (*((*(this->level_nodes))[l]))[node_id]->j;

        Vi->push_back(node_index);
        Vj->push_back(node_index);
        V_data->push_back(y * EY);
      }

      for(int k = 0; k<N_exterior_dn; k++){
        counter += 1;
        int node_index = k + N_interior ;
        string node_id = (*exterior_dn_node_ids)[k];
        my_map<string, node *> neighs = *(*((*(this->level_nodes))[l]))[node_id]->neighbors;
        for(auto const & neigh_pair: neighs){
          counter += 1;
          int neigh_index;
          if(interior_node_indices->find(neigh_pair.first) != interior_node_indices->end()){
            neigh_index = (*interior_node_indices)[neigh_pair.first];
          }
          else if(exterior_dn_node_indices->find(neigh_pair.first) != exterior_dn_node_indices->end()){
            neigh_index = (*exterior_dn_node_indices)[neigh_pair.first] + N_interior;
          }
          else if(exterior_up_node_indices->find(neigh_pair.first) != exterior_up_node_indices->end()){
            neigh_index = (*exterior_up_node_indices)[neigh_pair.first] + N_interior + N_exterior_dn;
          }
          Ci->push_back(node_index);
          Cj->push_back(neigh_index);
          C_data->push_back(1.0);
        }
      }

      for(int k = 0; k<N_exterior_up; k++){
        counter += 1;
        string node_id = (*exterior_up_node_ids)[k];
        int node_index = k + N_interior + N_exterior_dn;
        my_node_map neighs = *(*((*(this->level_nodes))[l]))[node_id]->neighbors;
        for(auto const & neigh_pair: neighs){
          counter += 1;
          int neigh_index;
          if(interior_node_indices->find(neigh_pair.first) != interior_node_indices->end()){
            neigh_index = (*interior_node_indices)[neigh_pair.first];
          }
          else if(exterior_dn_node_indices->find(neigh_pair.first) != exterior_dn_node_indices->end()){
            neigh_index = (*exterior_dn_node_indices)[neigh_pair.first] + N_interior;
          }
          else if(exterior_up_node_indices->find(neigh_pair.first) != exterior_up_node_indices->end()){
            neigh_index = (*exterior_up_node_indices)[neigh_pair.first] + N_interior + N_exterior_dn;
          }
          Ci->push_back(node_index);
          Cj->push_back(neigh_index);
          C_data->push_back(1.0);
        }
      }
      (*(this->level_is_build))[l] = true;
      cout << counter << endl;
    }
    l -= 1;
  }

}

