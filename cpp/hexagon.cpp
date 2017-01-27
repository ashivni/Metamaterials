#include "hierarchicalMaterials.h"


hexagon::hexagon(int i, int j, int a){
  node * n1 = new node(i,j);
  this->nodes.push_back(n1);

  n1 = new node(i + 2*a, j);
  this->nodes.push_back(n1);

  n1 = new node(i + a, j + a);
  this->nodes.push_back(n1);

  n1 = new node(i - a, j + a);
  this->nodes.push_back(n1);

  n1 = new node(i - 2*a, j);
  this->nodes.push_back(n1);

  n1 = new node(i - a, j - a);
  this->nodes.push_back(n1);

  n1 = new node(i + a, j - a);
  this->nodes.push_back(n1);

  int k = 0;
  for (k=0; k<6; k++)
  {
    bond *b = new bond(this->nodes[0], this->nodes[k+1], "I");
    this->bonds.push_back(b);
	}

  for (k=0; k<6; k++)
  {
    bond *b = new bond(this->nodes[1 + (k+1)%6], this->nodes[1+(k+2)%6], "I");
    this->bonds.push_back(b);
	}
}
