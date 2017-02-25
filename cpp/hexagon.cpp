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

coordinates_2d hexagon::coordinates(){
  coordinates_2d coords;

  for(auto const& value: this->nodes) {
    coords.append(value->coordinates());
  }
  return coords;
}

std::ostream& operator<< (std::ostream &out, const hexagon &h){
  out << "Hexagon: i, j, a: " << h.i << ", " << h.i << ", " << h.a <<endl;
  out << h.nodes.size() << " Nodes:" << endl;
  for(auto const& n: h.nodes) {
    out << "\t" << (*n);
  }
  out << h.bonds.size() << " Bonds:" << endl;
  for(auto const& b: h.bonds) {
    out << "\t" << (*b);
  }
  return out;
}

void hexagon::destruct(){
  for(auto const & t: this->nodes){
    delete t;
  }

  for(auto const & t: this->bonds){
    delete t;
  }
}

hexagon::~hexagon(){
  this->destruct();
}
