#include "hierarchicalMaterials.h"

node::node(int i, int j){
  this->i = i;
  this->j = j;
  this->n_neighbors = 0;
  ostringstream stringStream;
  stringStream << i << "_" << j;
  this->id = stringStream.str();
};

node::node(const node &obj){
  this->i = obj.i;
  this->j = obj.j;
  this->id = obj.id;
};

array <float, 2> node::coordinates(){
  array <float, 2> coords;
  coords[0] = 0.5*this->i;
  coords[1] = 0.5*this->j*pow(3.0,0.5);
  return coords;
}

bool node::is_connected_to(node *target){
  // Check if there is a path from this node to the target node
  deque<node *> dq;
  dq.push_front(this);
  set<node *> visited, queued;

  while (dq.size() > 0){
    node *n = dq.front();
    dq.pop_front();
    if (n == target) return true;
    else  visited.insert(n);

    for(auto const& x : n->neighbors){
      if( (visited.find(x.second) == visited.end()) && (queued.find(x.second) == queued.end())){
        dq.push_back(x.second);
        visited.insert(x.second);
      }
    }
  }

  return false;
};

void node::debond_from(node *n1){
  if (this->bonded_to(n1)){
    this->neighbors.erase(n1->get_id());
    this->n_neighbors -= 1;
  }

  if (n1->bonded_to(this)){
    n1->debond_from(this);
  }
}

bool node::bonded_to(node *n1){
  return this->neighbors.find(n1->get_id()) != this->neighbors.end();
};

void node::add_neighbor(node *n1){
  if(! this->bonded_to(n1)){
    this->neighbors.insert(std::make_pair(n1->get_id(),n1));
    this->n_neighbors +=1;
  }

  if(! n1->bonded_to(this)){
    n1->add_neighbor(this);
  }
};

string node::get_id(){
  return this->id;
};
