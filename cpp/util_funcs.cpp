#include "hierarchicalMaterials.h"

vector < vector <double> > stress(hierarchical_grid hg, int l){
  int i_min = numeric_limits<int>::max(), j_min = numeric_limits<int>::max();
  int i_max = numeric_limits<int>::min(), j_max = numeric_limits<int>::min();

  for(auto const & nm : *(*(hg.level_nodes))[l]){
    int x = nm.second->i, y = nm.second->j;
    i_min = i_min < x ? i_min : x;
    i_max = i_max > x ? i_max : x;
    j_min = j_min < y ? j_min : y;
    j_max = j_max > y ? j_max : y;
  }

  int notch_len_ind = round((hg.notch_len)*(i_max - i_min));
  double j_mid = 0.5*(j_min + j_max);
  if (abs(round(j_mid) - j_mid) < 1E-3){
    j_mid += 0.25;
  }

  vector <double> dist, curr;

  for(auto const & lb: (*((*(hg.level_bonds))[l])) ){
    int n1_x = lb.second->n1->i, n1_y = lb.second->n1->j;
    int n2_x = lb.second->n2->i, n2_y = lb.second->n2->j;
    if(n1_x - i_min >= notch_len_ind && n2_x - i_min > notch_len_ind && (n1_y - j_mid)*(n2_y - j_mid) < 0){
      dist.push_back(0.25*(n1_x + n2_x -2*notch_len_ind));
      curr.push_back(abs(hg.bond_current(lb.second,l)));
    }
  }

  // Insertion sort
  for (int k = 0; k < curr.size(); k++){
    int min_pos = k;
    for (int m = k; m < curr.size(); m++){
      min_pos = dist[m] < dist[min_pos] ? m : min_pos;
    }
    double temp;
    temp = dist[min_pos];
    dist[min_pos] = dist[k];
    dist[k] = temp;
    temp = curr[min_pos];
    curr[min_pos] = curr[k];
    curr[k] = temp;

  }

  vector<vector<double>> res;
  res.push_back(dist);
  res.push_back(curr);

  return res;
}