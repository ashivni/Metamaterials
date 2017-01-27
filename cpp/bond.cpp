#include "hierarchicalMaterials.h"

bond::bond(node *n1, node *n2, string bond_type){
  this->n1 = n1;
  this->n2 = n2;
  ostringstream stringStream;
  stringStream << this->n1->get_id() << "_" << this->n2->get_id();
  this->id = stringStream.str();
  this->bond_type = bond_type;
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


array < vector < array<float, 2> >,2> bond::coordinates(int nx, int ny){
  array < vector < array<float, 2> >,2>  coords;
  array <float, 2> c1 = this->n1->coordinates(), c2 = this->n2->coordinates();
  float rt3 = pow(3.0,0.5);
  if (this->bond_type == "R" || this->bond_type == "U" || this->bond_type == "UR" || this->bond_type == "I"){
    array<float, 2> x0, y0;
    x0[0] = c1[0];  x0[1] = c2[0];
    y0[0] = c1[1];  y0[1] = c2[1];
    coords[0].push_back(x0);
    coords[1].push_back(y0);
  }
  else if(this->bond_type == "PR"){
    array<float, 2> x00, x01, y00, y01;
    x00[0] = c1[0];                       x00[1] = (c1[0] + nx + c2[0])*0.5;
    x01[0] = (c1[0] - nx + c2[0])*0.5;    x01[1] = c2[0];
    coords[0].push_back(x00);
    coords[0].push_back(x01);

    y00[0] = c1[1];                       y00[1] = (c1[1] + c2[1])*0.5;
    y01[0] = (c1[1] + c2[1])*0.5;         y01[1] = c2[1];
    coords[1].push_back(y00);
    coords[1].push_back(y01);

  }
  else if(this->bond_type == "PU"){
    array<float, 2> x00, x01, y00, y01;
    x00[0] = c1[0];                       x00[1] = (c1[0] + c2[0])*0.5;
    x01[0] = (c1[0] + c2[0])*0.5;         x01[1] = c2[0];
    coords[0].push_back(x00);
    coords[0].push_back(x01);

    y00[0] = c1[1];                               y00[1] = (c1[1] + c2[1] + ny*rt3)*0.5;
    y01[0] = (c1[1] + c2[1] - ny*rt3)*0.5;        y01[1] = c2[1];
    coords[1].push_back(y00);
    coords[1].push_back(y01);
  }
  else if(this->bond_type == "C"){
    array<float, 2> x00, x01, y00, y01;
    x00[0] = c1[0];                     x00[1] = (c1[0] + nx + c2[0])*0.5;
    x01[0] = (c1[0] -nx + c2[0])*0.5;       x01[1] = c2[0];
    coords[0].push_back(x00);
    coords[0].push_back(x01);

    y00[0] = c1[1];                           y00[1] = (c1[1] + c2[1] + ny*rt3)*0.5;
    y01[0] = (c1[1] + c2[1] - ny*rt3)*0.5;    y01[1] = c2[1];
    coords[1].push_back(y00);
    coords[1].push_back(y01);
  }
  else if(this->bond_type == "PUR"){
    array<float, 2> x00, x01, y00, y01;
    x00[0] = c1[0];                     x00[1] = (c1[0] + nx + c2[0])*0.5;
    x01[0] = (c1[0] -nx + c2[0])*0.5;       x01[1] = c2[0];
    coords[0].push_back(x00);
    coords[0].push_back(x01);

    y00[0] = c1[1];                           y00[1] = (c1[1] + c2[1])*0.5;
    y01[0] = (c1[1] + c2[1])*0.5;    y01[1] = c2[1];
    coords[1].push_back(y00);
    coords[1].push_back(y01);
  }

	return coords;
}

nodes_and_bonds bond::refine(int magnification){
  array<float, 2> c1 = this->n1->coordinates(), c2 = this->n2->coordinates();
  array<int, 2> t_vec = {(this->n2->i - this->n1->i)/magnification, (this->n2->j - this->n1->j)/magnification};
  array<float, 2> t_vec_real = { (c2[0] - c1[0])/magnification, (c2[1] - c1[1])/magnification };
  int l = (int) round(pow(pow(t_vec_real[0], 2.0) + pow(t_vec_real[1], 2.0), 0.5 ));
  nodes_and_bonds nb;
  int k = 0;
  for(k=0; k<magnification; k++){
    hexagon h(this->n1->i + k*t_vec[0], this->n1->j + k*t_vec[1], l);
    for(auto const& n : h.nodes){
      nb.nodes[n->id] = n;
    }

    for(auto const& b : h.bonds){
      nb.bonds[b->id] = b;
      nb.nodes[b->n1->id]->add_neighbor(nb.nodes[b->n2->id]);
    }
  }

  return nb;
}

void bond::sever(){
  this->n1->debond_from(this->n2);
}