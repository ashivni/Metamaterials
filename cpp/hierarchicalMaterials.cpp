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
  node n1(1,2), n2(2,3), n3(3,4);
  cout << n1.get_id() << endl;
  cout << n1.coordinates()[0] << ", " << n1.coordinates()[1] << endl;

  check_bond(&n1, &n2);
  n1.add_neighbor(&n2);
  check_bond(&n1,&n2);
  n1.debond_from(&n2);
  check_bond(&n1,&n2);
  n2.add_neighbor(&n1);
  check_bond(&n2,&n1);
  n2.add_neighbor(&n3);
  check_connectivity(&n1, &n3);

 }