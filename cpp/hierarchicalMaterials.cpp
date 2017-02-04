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
int nx = 5, ny=5, levels=2, l0=1, magnification=6;

  hierarchical_grid hg(nx, ny, levels, l0, magnification,false, false, 0.0);

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
      if (l==0){
      fOut << l << "\t" << ni;
      for (auto const & neigh: neigh_ids){
        fOut << "\t" << neigh;
      }
      fOut << endl;
      }
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