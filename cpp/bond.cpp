#include "hierarchicalMaterials.h"

bond::bond(node *n1, node *n2, string bond_type){
  this->n1 = n1;
  this->n2 = n2;
  ostringstream stringStream;
  stringStream << this->n1->get_id() << "_" << this->n2->get_id();
  this->id = stringStream.str();
  this->bond_type = bond_type;
  this->n1->add_neighbor(this->n2);
};

bond::bond(const bond &obj){
  this->n1 = obj.n1;
  this->n2 = obj.n2;
  this->id = obj.id;
  this->bond_type = obj.bond_type;
};

string bond::get_id(){
  return this->id;
};


vector<coordinates_2d> bond::coordinates(int nx, int ny){
  vector<coordinates_2d>  coords;
  vector<double> c3 = {1.0, 1.0};
  vector<double> c1 = {this->n1->coordinates().x[0], this->n1->coordinates().y[0]};
  vector<double> c2 = {this->n2->coordinates().x[0], this->n2->coordinates().y[0]};
  double rt3 = pow(3.0,0.5);
  if (this->bond_type == "R" || this->bond_type == "U" || this->bond_type == "UR" || this->bond_type == "I"){
    vector<double> x0 = {c1[0], c2[0]};
    vector<double> y0 = {c1[1], c2[1]};
    coordinates_2d c(x0, y0);
    coords.push_back(c);
  }
  else if(this->bond_type == "PR"){
    vector<double> x00 = {c1[0], (c1[0] + nx + c2[0])*0.5};
    vector<double> x01 = {(c1[0] - nx + c2[0])*0.5, c2[0]};
    vector<double> y00 = {c1[1], (c1[1] + c2[1])*0.5};
    vector<double> y01 = {(c1[1] + c2[1])*0.5, c2[1]};
    coordinates_2d c1(x00, y00), c2(x01, y01);
    coords.push_back(c1);
    coords.push_back(c2);
  }
  else if(this->bond_type == "PU"){
    vector<double> x00 = {c1[0], (c1[0] + c2[0])*0.5};
    vector<double> x01 = {(c1[0] + c2[0])*0.5, c2[0]};
    vector<double> y00 = {c1[1], (c1[1] + c2[1] + ny*rt3)*0.5};
    vector<double> y01 = {(c1[1] + c2[1] -ny*rt3)*0.5, c2[1]};
    coordinates_2d c1(x00, y00), c2(x01, y01);
    coords.push_back(c1);
    coords.push_back(c2);
  }
  else if(this->bond_type == "C"){
    vector<double> x00 = {c1[0], (c1[0] + nx + c2[0])*0.5};
    vector<double> x01 = {(c1[0] - nx + c2[0])*0.5, c2[0]};
    vector<double> y00 = {c1[1], (c1[1] + c2[1] + ny*rt3)*0.5};
    vector<double> y01 = {(c1[1] + c2[1] -ny*rt3)*0.5, c2[1]};
    coordinates_2d c1(x00, y00), c2(x01, y01);
    coords.push_back(c1);
    coords.push_back(c2);
  }
  else if(this->bond_type == "PUR"){
    vector<double> x00 = {c1[0], (c1[0] + nx + c2[0])*0.5};
    vector<double> x01 = {(c1[0] - nx + c2[0])*0.5, c2[0]};
    vector<double> y00 = {c1[1], (c1[1] + c2[1])*0.5};
    vector<double> y01 = {(c1[1] + c2[1])*0.5, c2[1]};
    coordinates_2d c1(x00, y00), c2(x01, y01);
    coords.push_back(c1);
    coords.push_back(c2);
  }

	return coords;
}

nodes_and_bonds bond::refine(int magnification){
  array<double, 2> c1 = {this->n1->coordinates().x[0], this->n1->coordinates().y[0]};
  array<double, 2> c2 = {this->n2->coordinates().x[0], this->n2->coordinates().y[0]};
  array<int, 2> t_vec = {(this->n2->i - this->n1->i)/magnification, (this->n2->j - this->n1->j)/magnification};
  array<double, 2> t_vec_real = { (c2[0] - c1[0])/magnification, (c2[1] - c1[1])/magnification };
  int l = (int) round(pow(pow(t_vec_real[0], 2.0) + pow(t_vec_real[1], 2.0), 0.5 ));
  nodes_and_bonds nb;
  int k = 0;
  for(k=0; k<magnification; k++){
    hexagon * h = new hexagon(this->n1->i + k*t_vec[0], this->n1->j + k*t_vec[1], l);
    for(auto const& n : h->nodes){
      //if (nb.nodes.find(n->id) ==  nb.nodes.end()){
        nb.nodes[n->id] = n;
      //}
    }

    for(auto const& b : h->bonds){
      //if (nb.bonds.find(b->id) ==  nb.bonds.end()){
        nb.bonds[b->id] = b;
        //try{
          nb.nodes[b->n1->id]->add_neighbor(nb.nodes[b->n2->id]);
        //} catch (const invalid_argument & e){}
      //}
    }
  }

  return nb;
}

void bond::sever(){
  this->n1->debond_from(this->n2);
}

std::ostream& operator<< (std::ostream &out, const bond &b){
  out << "Bond ID, type: " << b.id << ", " << b.bond_type << endl;
  out << "\tNode 1: " << *(b.n1);
  out << "\tNode 2: " << *(b.n2);
  return out;
}
