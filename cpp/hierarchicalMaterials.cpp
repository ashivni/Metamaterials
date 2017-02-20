#include "hierarchicalMaterials.h"

void check_bond(node * n1, node * n2){
  if (n1->bonded_to(n2)) {
    cout << n1->get_id() << " bonded to " << n2->get_id() << endl;
  }
  else{
    cout << n1->get_id() << " not bonded to " << n2->get_id() << endl;
  }

  if (n2->bonded_to(n1)) {
    cout << n2->get_id() << " bonded to " << n1->get_id() << endl;
  }
  else{
    cout << n2->get_id() << " not bonded to " << n1->get_id() << endl;
  }
}

void check_connectivity(node *n1, node * n2){
  if(n1->is_connected_to(n2)){
    cout << n1->get_id() << " connected to " << n2->get_id() << endl;
  }
  else{
    cout << n1->get_id() << " not connected to " << n2->get_id() << endl;
  }
  if(n2->is_connected_to(n1)){
    cout << n2->get_id() << " connected to " << n1->get_id() << endl;
  }
  else{
    cout << n2->get_id() << " not connected to " << n1->get_id() << endl;
  }
}

int main(){
std::cout.setf( std::ios_base::unitbuf );

/*
int m = 20, n = 30, nnz = 0;
vector<int> i, j;
vector <double> A;
srand(1);
for (int k = 0; k < m; k++){
  for (int l = 0; l < n; l++){
    if (rand()%100 > 90){
      nnz += 1;
      i.push_back(k);
      j.push_back(l);
      A.push_back((rand()%100)*1.0);
    }
  }
}
CPPSparse cs = CPPSparse(m,n,nnz,&i,&j,&A);
for (int k = 0; k < m; k++){
  for (int l = 0; l < n; l++){
    int d_index = cs.data_index(k, l);
    int d_index_loc = -1;
    for(int t = 0; t < i.size(); t++){
      if (i[t] == k && j[t] == l){
        d_index_loc = t;
      }
    }

    if (d_index == -1 && d_index_loc == -1){}
    else if (d_index == -1){cout << "error"<< endl;}
    else if (d_index_loc == -1){cout << "error" << endl;}
    else if (cs.cs_c->x[d_index] != A[d_index_loc]) {cout << "error" << endl;}
  }
}

int m = 3, n = 3, nnz = 4;
vector<int> i = {0,1,2,0}, j = {0,1,2,2};
vector <double> A = {1, 1, 1, 2.0};
CPPSparse cs = CPPSparse(m,n,nnz,&i,&j,&A);
vector <double> b = {1, 2, 3};
cs.CPPSprint();
for(auto const & t: cs.sp_index_t){
  cout << t.first << ":" << t.second << endl;
}
cs.update(0,2,0.0);
cs.CPPSprint();
vector <double> x = cs.solve(b);
cout << "x: " << "\t";
for(auto const & pp: x){
  cout << pp << "\t";
}
cout << endl;
return 1;
*/
/*
linear_system ls = linear_system(m,n,nnz,&i,&j,&A,&b);

x = ls.solve();
cout << "x: " << "\t";
for(auto const & pp: x){
  cout << pp << "\t";
}
cout << endl;

return 1;
*/

int nx=17, ny=10, levels=1, l0=1, magnification=6;

  hierarchical_grid hg(nx, ny, levels, l0, magnification,true, false, 0.0, 1.0);
  hg.build_eqns();
  //hg.solve();
  //hg.step();
  //hg.step();
  cout << '\a';
  //hg.simulate_fracture(0);
  cout << "Steps to break: " << hg.n_step << endl;
  hg.dump();


  int l = hg.levels;
  ofstream fOut;
  fOut.open("C_Out.txt");
  while(l >= 0){
    vector <string> node_ids;
    map <string, node *> nm = *(*(hg.level_nodes))[l];
    for(auto const & np: nm){
      node_ids.push_back(np.second->id);
    }
    sort(node_ids.begin(),node_ids.end());

    for(auto const & ni:node_ids){
      vector <string> neigh_ids;
      for(auto const & neigh_map: *(nm[ni]->neighbors)){
        neigh_ids.push_back(neigh_map.second->id);
      }
      sort(neigh_ids.begin(), neigh_ids.end());
      fOut << l << "\t" << ni;
      for (auto const & neigh: neigh_ids){
        fOut << "\t" << neigh;
      }
      fOut << endl;
    }

    vector <string> bond_ids;
    map <string, bond *> bm = *(*(hg.level_bonds))[l];
    for(auto const & bp: bm){
      bond_ids.push_back(bp.second->id);
    }
    sort(bond_ids.begin(),bond_ids.end());

    for(auto const & bi:bond_ids){
      fOut << l << "\t" << bi << endl;
    }
    l--;
  }
  fOut.close();

/*
  hexagonal_grid hg (nx, ny, levels, magnification, l0);
  ofstream fOut;
  fOut.open("C_Out.txt");
  fOut << "nx: " << hg.nx << "\t";
  fOut << "ny: " << hg.ny << "\t";
  fOut << "level: " << hg.level << "\t";
  fOut << "mag: " << hg.magnification << "\t";
  fOut << "l0: " << hg.l0 << "\t";
  fOut << "a: " << hg.a << endl ;

  vector <string> node_ids;
  map <string, node *> nm = *(hg.node_map());
  for(auto const & np: nm){
    node_ids.push_back(np.second->id);
  }
  sort(node_ids.begin(),node_ids.end());

  for(auto const & ni:node_ids){
    vector <string> neigh_ids;
    for(auto const & neigh_map: nm[ni]->neighbors){
      neigh_ids.push_back(neigh_map.second->id);
    }
    sort(neigh_ids.begin(), neigh_ids.end());

    fOut << ni;
    for (auto const & neigh: neigh_ids){
      fOut << "\t" << neigh;
    }
    fOut << endl;
  }

  fOut.close();
*/
/*
  hexagonal_unit_cell huc (1,1,4);
  ofstream fOut;
  fOut.open("C_Out.txt");
  fOut << "i: " << huc.i << "\t";
  fOut << "j: " << huc.j << "\t";
  fOut << "a: " << huc.a << endl;

  for(auto const & np: huc.nodes){
    fOut << "Node:\t" << np->id << "\t";
    for(auto const & mn: np->neighbors){
      fOut << mn.second->id << "\t";
    }
    fOut << endl;
  }

  for(auto const & bp: huc.interior_bonds){
    fOut << "IB:\t" << bp->id << "\t";
    fOut << endl;
  }
  for(auto const & bp: huc.exterior_bonds){
    fOut << "EB:\t" << bp->id << "\t";
    fOut << endl;
  }
  fOut.close();
*/
 }
