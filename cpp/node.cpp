#include "hierarchicalMaterials.h"

node::node(int i, int j){
  this->i = i;
  this->j = j;
  ostringstream stringStream;
  stringStream << "_" << i << "_" << j;
  this->id = stringStream.str();
  this->neighbors = new map<string, node *>;
};

node::node(const node &obj){
  this->i = obj.i;
  this->j = obj.j;
  this->id = obj.id;
  this->neighbors = obj.neighbors;
};

coordinates_2d node::coordinates() const {
  coordinates_2d coords;
  coords.x.push_back(0.5*this->i);
  coords.y.push_back(0.5*this->j*pow(3.0,0.5));
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

    for(auto const& x : *(n->neighbors)){
      if( (visited.find(x.second) == visited.end()) && (queued.find(x.second) == queued.end())){
        dq.push_back(x.second);
        visited.insert(x.second);
      }
    }
  }

  return false;
};

void node::debond_from(node *n1){
  /*if(!this->bonded_to(n1) && !n1->bonded_to(this)){
      ostringstream stringStream;
      stringStream << "[debond_from] Nodes not bonded: " << this->id << " " << n1->id << endl;
      throw invalid_argument(stringStream.str());
  }

  if(!this->bonded_to(n1) || !n1->bonded_to(this)){
      ostringstream stringStream;
      stringStream << "[debond_from] Nodes partially bonded: " << this->id << " " << n1->id << endl;
      throw invalid_argument(stringStream.str());
  }*/
    this->neighbors->erase(n1->get_id());
    n1->neighbors->erase(this->get_id());
}


bool node::bonded_to(node *n1){
  return this->neighbors->find(n1->get_id()) != this->neighbors->end();
};

void node::add_neighbor(node *n1){
  /*if(this->bonded_to(n1) && n1->bonded_to(this)){
      ostringstream stringStream;
      stringStream << "[add_neighbor] Nodes already bonded: " << this->id << " " << n1->id << endl;
      throw invalid_argument(stringStream.str());
  }

  if(this->bonded_to(n1) || n1->bonded_to(this)){
      ostringstream stringStream;
      stringStream << "[add_neighbor] Nodes partially bonded: " << this->id << " " << n1->id << endl;
      throw invalid_argument(stringStream.str());
  }*/
  (*(this->neighbors))[n1->get_id()] = n1;
  (*(n1->neighbors))[this->get_id()] = this;
};

string node::get_id(){
  return this->id;
};

std::ostream& operator<< (std::ostream &out, const node &n){
  out << "Node ID, i, j: " << n.id << ",\t" << n.i << ",\t" << n.j << ", coordinates: " << n.coordinates();
  return out;
}
