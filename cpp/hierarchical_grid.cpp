#include "hierarchicalMaterials.h"

hierarchical_grid::hierarchical_grid(int nx, int ny, int levels, int l0, int magnification, bool notch, bool damage, double damage_fraction, double EY, double notch_len){
  this->nx = nx;
  this->ny = ny;
  this->levels = levels;
  this->l0 = l0;
  this->magnification = magnification;
  this->notch = notch;
  this->damage = damage;
  this->damage_fraction = damage_fraction;
  this->EY = EY;
  this->notch_len = notch_len;
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
  this->level_is_solved = new map<int, bool>;
  this->level_not_broken = new map<int, bool>;
  this->level_C = new map<int, CPPSparse *>;
  this->level_B = new map<int, CPPSparse *>;
  this->level_L = new map<int, CPPSparse *>;
  this->level_V = new map<int, CPPSparse *>;
  this->level_curr = new map<int, CPPSparse *>;
  this->level_curr_dn = new map<int, double>;
  this->level_stress = new map<int, vector<double> *>;
  this->level_strain = new map<int, vector<double> *>;
  this->level_voltage = new map<int, vector<double> *>;
  this->level_current = new map<int, vector<double> *>;
  this->level_scaled_stress = new map<int, vector<double> *>;

  this->level_interior_node_indices = new  map<int, map <string, int> *>;
  this->level_exterior_up_node_indices = new  map<int, map <string, int> *>;
  this->level_exterior_dn_node_indices = new  map<int, map <string, int> *>;
  this->level_interior_node_ids = new  map<int, vector<string> *>;
  this->level_exterior_up_node_ids = new  map<int, vector <string> *>;
  this->level_exterior_dn_node_ids = new  map<int, vector <string> *>;
  this->level_int_vars = new map<int, map<string, int> *>;

  (*(this->level_is_build))[this->levels] = false;
  (*(this->level_is_solved))[this->levels] = false;
  (*(this->level_not_broken))[this->levels] = true;
  (*(this->level_nodes))[this->levels] = this->outline->node_map();
  (*(this->level_bonds))[this->levels] = this->outline->bond_map();
  (*(this->level_broken_bonds))[this->levels] = new map<string, bond *>;

  (*(this->level_stress))[this->levels] = new vector<double>;
  (*(this->level_strain))[this->levels] = new vector<double>;
  (*(this->level_voltage))[this->levels] = new vector<double>;
  (*(this->level_current))[this->levels] = new vector<double>;
  (*(this->level_scaled_stress))[this->levels] = new vector<double>;


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
        //if(bonds->find(refined_bond_pts.first) == bonds->end()){
          (*bonds)[refined_bond_pts.first] = refined_bond_pts.second;
        //}
        node *n1 = (*nodes)[refined_bond_pts.second->n1->id];
        node *n2 = (*nodes)[refined_bond_pts.second->n2->id];
        //if(!n1->bonded_to(n2))  n1->neighbors->insert(std::make_pair(n2->get_id(),n2));
        //if(!n2->bonded_to(n1))  n2->neighbors->insert(std::make_pair(n1->get_id(),n1));
        (*(n1->neighbors))[n2->get_id()] = n2;
        (*(n2->neighbors))[n1->get_id()] = n1;
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
      if (it == duplicate_bonds.end() && (b.second)->n1->id > (b.second)->n2->id){
        string dup_id = (b.second)->n2->id + "_" + (b.second)->n1->id;
        if ((*((*(this->level_bonds))[l])).find(dup_id) !=  (*((*(this->level_bonds))[l])).end()){
          duplicate_bonds.insert(dup_id);
        }
      }
    }

    for(auto const& db: duplicate_bonds){
      (*((*(this->level_bonds))[l])).erase(db);
    }
    (*(this->level_is_build))[l] = false;
    (*(this->level_is_solved))[l] = false;
    (*(this->level_not_broken))[l] = true;
    (*(this->level_broken_bonds))[l] = new map<string, bond *>;
    
    (*(this->level_stress))[l] = new vector<double>;
    (*(this->level_strain))[l] = new vector<double>;
    (*(this->level_voltage))[l] = new vector<double>;
    (*(this->level_current))[l] = new vector<double>;
    (*(this->level_scaled_stress))[l] = new vector<double>;


    l-=1;
  }

  if(this->notch){
    int l = this->levels;
    this->level_notched_bonds = new map<int, map<string, bond *> *>;
    this->level_notch_len = new map<int, double>;

    while(l >= 0){
      int i_min = numeric_limits<int>::max(), j_min = numeric_limits<int>::max();
      int i_max = numeric_limits<int>::min(), j_max = numeric_limits<int>::min();

      for(auto const & nm : *(*(this->level_nodes))[l]){
        int x = nm.second->i, y = nm.second->j;
        i_min = i_min < x ? i_min : x;
        i_max = i_max > x ? i_max : x;
        j_min = j_min < y ? j_min : y;
        j_max = j_max > y ? j_max : y;
      }
      int notch_len_ind = round((this->notch_len)*(i_max - i_min)) + (this->l0)*pow(this->magnification,(this->levels));
      double j_mid = 0.5*(j_min + j_max);
      if (abs(round(j_mid) - j_mid) < 1E-3){
        j_mid += 0.25;
      }

      map<string, bond *> * notched_bonds = new map<string, bond *>;
      for(auto const & nm : *(*(this->level_nodes))[l]){
        int node_x = nm.second->i, node_y = nm.second->j;
        map<string, node *> neigh_map;
        for(auto const & neigh: *(nm.second->neighbors)){
          neigh_map[neigh.first]  = neigh.second;
        }
        for(auto const & neigh: neigh_map){
          int neigh_x = neigh.second->i, neigh_y = neigh.second->j;
          if(node_x - i_min <= notch_len_ind && (node_y - j_mid)*(neigh_y - j_mid) < 0){
            nm.second->debond_from(neigh.second);
            ostringstream s1, s2;
            s1 << nm.first << "_" << neigh.first;
            string b12_id = s1.str();
            s2 << neigh.first << "_" << nm.first;
            string b21_id = s2.str();
            if((*(this->level_bonds))[l]->find(b12_id) != (*(this->level_bonds))[l]->end()){
              (*notched_bonds)[b12_id] = (*(*(this->level_bonds))[l])[b12_id];
              (*(this->level_bonds))[l]->erase(b12_id);
            }
            if((*(this->level_bonds))[l]->find(b21_id) != (*(this->level_bonds))[l]->end()){
              (*notched_bonds)[b21_id] = (*(*(this->level_bonds))[l])[b21_id];
              (*(this->level_bonds))[l]->erase(b21_id);
            }
          }
        }
      }
      (*(this->level_notched_bonds))[l] = notched_bonds;
      (*(this->level_notch_len))[l] = notch_len_ind*0.5;
      l -= 1;
    }
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
        map <string, node *> *neighs = (*((*(this->level_nodes))[l]))[node_id]->neighbors;
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
        map<string, node *> neighs = *(*((*(this->level_nodes))[l]))[node_id]->neighbors;
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
        map<string, node *> neighs = *(*((*(this->level_nodes))[l]))[node_id]->neighbors;
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
      int n_nodes = N_interior + N_exterior_dn + N_exterior_up;
      (*(this->level_C))[l] = new CPPSparse(n_nodes, n_nodes, C_data->size(), Ci, Cj, C_data);
      (*(this->level_L))[l] = new CPPSparse(N_interior, N_interior, L_data->size(), Li, Lj, L_data);
      (*(this->level_B))[l] = new CPPSparse(N_interior, 1, B_data->size(), Bi, Bj, B_data);
      (*(this->level_V))[l] = new CPPSparse(n_nodes, n_nodes, V_data->size(), Vi, Vj, V_data);
      (*(this->level_interior_node_indices))[l] = interior_node_indices;
      (*(this->level_exterior_dn_node_indices))[l] = exterior_dn_node_indices;
      (*(this->level_exterior_up_node_indices))[l] = exterior_up_node_indices;
      (*(this->level_interior_node_ids))[l] = interior_node_ids;
      (*(this->level_exterior_dn_node_ids))[l] = exterior_dn_node_ids;
      (*(this->level_exterior_up_node_ids))[l] = exterior_up_node_ids;
      map<string, int> * int_vars = new map<string, int>;
      (*int_vars)["N_interior"] = N_interior;
      cout << "Level, N_vars: " << l << ", " << N_interior << endl;
      (*int_vars)["N_exterior_dn"] = N_exterior_dn;
      (*int_vars)["N_exterior_up"] = N_exterior_up;
      (*int_vars)["i_range"] = i_max - i_min + 1;
      (*int_vars)["j_range"] = j_max - j_min + 1;
      (*int_vars)["i_offset"] = i_min;
      (*int_vars)["j_offset"] = j_min;
      (*(this->level_int_vars))[l] = int_vars;

      (*(this->level_is_build))[l] = true;
    }
    l -= 1;
  }
}

void hierarchical_grid::solve(){
  this->build_eqns();
  int l = this->levels;
  double EY = this->EY;

  while(l >= 0){
    if(!(*(this->level_is_solved))[l]){
      CPPSparse * L = (*(this->level_L))[l];
      CPPSparse * B = (*(this->level_B))[l];
      CPPSparse * C = (*(this->level_C))[l];
      CPPSparse * V = (*(this->level_V))[l];

      int N_interior = (*((*(this->level_int_vars))[l]))["N_interior"];
      int N_exterior_up = (*((*(this->level_int_vars))[l]))["N_exterior_up"];
      int N_exterior_dn = (*((*(this->level_int_vars))[l]))["N_exterior_dn"];

      double *temp = new double [B->cs_c->n];
      vector <double> B_vec;
      for(int k = 0; k < B->cs_c->n; k++){
        temp[k] = 0.0;
      }
      for(int k = 0; k < B->cs_c->nzmax; k++){
        temp[B->cs_c->i[k]] = B->cs_c->x[k];
      }
      for(int k = 0; k < B->cs_c->n; k++){
        B_vec.push_back(temp[k]);
      }

      vector < double > interior_vol = (*L).solve(B_vec);
      for(int k = 0; k < interior_vol.size(); k++){
        V->cs_c->x[k] = interior_vol[k];
      }

      CPPSparse curr = ((*V)*(*C) - (*C)*(*V));

      double curr_up = 0.0, curr_dn = 0.0;
      for(int k = 0; k < curr.cs_c->nzmax; k++){
        int row = curr.cs_c->i[k];
        if( row >= N_interior && row < N_interior + N_exterior_dn){
          curr_up += curr.cs_c->x[k];
        }
        if( row >= N_interior + N_exterior_dn && row < N_interior + N_exterior_dn + N_exterior_up){
          curr_dn += curr.cs_c->x[k];
        }
      }

      if( abs(curr_up) < 1E-2){
        (*(this->level_not_broken))[l] = false;
      }
      if(abs(curr_up + curr_dn) > 1E-2){
        throw runtime_error("Net current flux");
      }
      (*(this->level_curr))[l] = new CPPSparse(curr);
      (*(this->level_curr_dn))[l] = curr_dn;
      (*(this->level_is_solved))[l] = true;
    }
    l -= 1;
  }
}

void hierarchical_grid::step(){
  cout << "Building: " << endl;
  this->build_eqns();
  cout << "Solving: " << endl;
  this->solve();
  int l = this->levels;
  double EY = this->EY;
  int lx = this->lx, ly = this->ly;

  while(l >= 0){
    if((*(this->level_not_broken))[l]){
      CPPSparse * L = (*(this->level_L))[l];
      CPPSparse * B = (*(this->level_B))[l];
      CPPSparse * C = (*(this->level_C))[l];
      CPPSparse * V = (*(this->level_V))[l];
      CPPSparse * curr = (*(this->level_curr))[l];
      double curr_dn = (*(this->level_curr_dn))[l];

      int N_interior = (*((*(this->level_int_vars))[l]))["N_interior"];
      int N_exterior_up = (*((*(this->level_int_vars))[l]))["N_exterior_up"];
      int N_exterior_dn = (*((*(this->level_int_vars))[l]))["N_exterior_dn"];
      int N  = N_interior + N_exterior_dn + N_exterior_up;

      map <string, int> * interior_node_indices = (*(this->level_interior_node_indices))[l];
      map <string, int> * exterior_up_node_indices = (*(this->level_exterior_up_node_indices))[l];
      map <string, int> * exterior_dn_node_indices = (*(this->level_exterior_dn_node_indices))[l];

      map <string, node *> * level_nodes = (*(this->level_nodes))[l];
      map <string, bond *> * level_bonds = (*(this->level_bonds))[l];
      map <string, bond *> * level_broken_bonds = (*(this->level_broken_bonds))[l];

      vector <string> * interior_node_ids = (*(this->level_interior_node_ids))[l];
      vector <string> * exterior_up_node_ids = (*(this->level_exterior_up_node_ids))[l];
      vector <string> * exterior_dn_node_ids = (*(this->level_exterior_dn_node_ids))[l];

      double max_cur = -numeric_limits<double>::infinity();
      string n1_id, n2_id, b_id;
      for(auto const & lb: (*((*(this->level_bonds))[l])) ){
        string n1_b_id = lb.second->n1->id;
        string n2_b_id = lb.second->n2->id;
        int n1_ind = -1, n2_ind = -1;
        if(interior_node_indices->find(n1_b_id) != interior_node_indices->end()){
          n1_ind = (*interior_node_indices)[n1_b_id];
        }
        else if(exterior_dn_node_indices->find(n1_b_id) != exterior_dn_node_indices->end()){
          n1_ind = (*exterior_dn_node_indices)[n1_b_id] + N_interior;
        }
        else{
          n1_ind = (*exterior_up_node_indices)[n1_b_id] + N_interior + N_exterior_dn;
        }

        if(interior_node_indices->find(n2_b_id) != interior_node_indices->end()){
          n2_ind = (*interior_node_indices)[n2_b_id];
        }
        else if(exterior_dn_node_indices->find(n2_b_id) != exterior_dn_node_indices->end()){
          n2_ind = (*exterior_dn_node_indices)[n2_b_id] + N_interior;
        }
        else{
          n2_ind = (*exterior_up_node_indices)[n2_b_id] + N_interior + N_exterior_dn;
        }

        if( abs(curr->get(n1_ind, n2_ind)) > max_cur){
          max_cur = abs(curr->get(n1_ind, n2_ind));
          n1_id = n1_b_id;
          n2_id = n2_b_id;
          b_id = lb.first;
        }
      }

      if (abs(max_cur - curr->max_data()) > 1E-3){
        throw range_error("Max not consistent.");
      }
      // Debond
      (*level_nodes)[n1_id]->debond_from((*level_nodes)[n2_id]);
      ostringstream s1, s2;
      s1 << n1_id << "_" << n2_id;
      string b12_id = s1.str();
      s2 << n2_id << "_" << n1_id;
      string b21_id = s2.str();
      if(level_bonds->find(b12_id) != level_bonds->end()){
        (*level_broken_bonds)[b12_id] = (*level_bonds)[b12_id];
        level_bonds->erase(b12_id);
      }
     // if (l == 0) throw invalid_argument("HERE");

      if(level_bonds->find(b21_id) != level_bonds->end()){
        (*level_broken_bonds)[b21_id] = (*level_bonds)[b21_id];
        level_bonds->erase(b21_id);
      }

      // Update matrices
      vector <string> break_node_id = {n1_id, n2_id};
      vector <string> break_neigh_id = {n2_id, n1_id};
      
      for(int k = 0; k < 2; k++){
        string node_id = break_node_id[k], neigh_id = break_neigh_id[k];
        int neigh_ind, node_ind;
        
        if(interior_node_indices->find(node_id) != interior_node_indices->end()){
          node_ind = (*interior_node_indices)[node_id];
        }
        else if(exterior_dn_node_indices->find(node_id) != exterior_dn_node_indices->end()){
          node_ind = (*exterior_dn_node_indices)[node_id] + N_interior;
        }
        else if(exterior_up_node_indices->find(node_id) != exterior_up_node_indices->end()){
          node_ind = (*exterior_up_node_indices)[node_id] + N_interior + N_exterior_dn;
        }

        if(interior_node_indices->find(neigh_id) != interior_node_indices->end()){
          neigh_ind = (*interior_node_indices)[neigh_id];
        }
        else if(exterior_dn_node_indices->find(neigh_id) != exterior_dn_node_indices->end()){
          neigh_ind = (*exterior_dn_node_indices)[neigh_id] + N_interior;
        }
        else if(exterior_up_node_indices->find(neigh_id) != exterior_up_node_indices->end()){
          neigh_ind = (*exterior_up_node_indices)[neigh_id] + N_interior + N_exterior_dn;
        }

        if(interior_node_indices->find(node_id) != interior_node_indices->end()){
          double load = 0.0;
          int diag = 1;
          if(find(interior_node_ids->begin(), interior_node_ids->end(), neigh_id) != interior_node_ids->end()){
            L->update(node_ind, neigh_ind, L->get(node_ind, neigh_ind) + 1);
          }
          else if(find(exterior_up_node_ids->begin(), exterior_up_node_ids->end(), neigh_id) != exterior_up_node_ids->end()){
            load += this->ly_ind*EY;
          }
          else if(find(exterior_dn_node_ids->begin(), exterior_dn_node_ids->end(), neigh_id) != exterior_dn_node_ids->end()){
            load += 0;
          }
          else{
            throw invalid_argument("Unknown node id");
          }

          if(load!=0.0){
            B->update(node_ind, 0, B->get(node_ind, 0) - load);
          }

          L->update(node_ind, node_ind, L->get(node_ind, node_ind) - diag);
        }
        C->update(node_ind, neigh_ind, C->get(node_ind, neigh_ind) - 1);
      }
      (*(this->level_is_solved))[l] = false;
      (*(this->level_voltage))[l]->push_back(EY/max_cur);
      (*(this->level_current))[l]->push_back(curr_dn/max_cur);
      (*(this->level_scaled_stress))[l]->push_back(curr_dn*ly/(max_cur*(*(this->level_bonds))[l]->size()));
      (*(this->level_stress))[l]->push_back(curr_dn/(lx*max_cur));
      (*(this->level_strain))[l]->push_back(EY/(max_cur));
    }
  l-=1;
  }
  this->n_step += 1;
}

bool hierarchical_grid::is_broken(){
  bool is_broken = true;
  for(auto const & t : *(this->level_not_broken)){
    if(t.second){
      is_broken = false;
      break;
    }
  }
  return is_broken;
}

void hierarchical_grid::simulate_fracture(int max_step = 0){
  if (max_step == 0){
    while(!this->is_broken()){
      this->step();
    }
  }
  else{
    while(!this->is_broken() && this->n_step < max_step){
      this->step();
    }
  }
}

void hierarchical_grid::dump(string pref){
  int l = this->levels;
  while(l >=0){
    ostringstream s;
    s << pref << "_" << l << "_stats.csv";
    string fName = s.str();

    ofstream outFile;
    outFile.open(fName,std::ios_base::out);

    outFile << "Voltage,Current,Scaled_stress,Stress,Strain" << endl;
    for(int k = 0; k < (*((*(this->level_voltage))[l])).size(); k++){
      outFile << (*((*(this->level_voltage))[l]))[k] << ",";
      outFile << (*((*(this->level_current))[l]))[k] << ",";
      outFile << (*((*(this->level_scaled_stress))[l]))[k] << ",";
      outFile << (*((*(this->level_stress))[l]))[k] << ",";
      outFile << (*((*(this->level_strain))[l]))[k] << endl;
    }
    outFile.close();
    l -= 1;
  }
}