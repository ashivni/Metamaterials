#include "hierarchicalMaterials.h"


int my_node_map::size(){
  return this->keys.size();
};

string my_node_map::key_at(int i){
  if(i < this->size()){
    return this->keys[i];
  }
  else{
    throw invalid_argument("Index out of range");
  }
};

node * my_node_map::value_at(int i){
  if(i < this->size()){
    return this->values[i];
  }
  else{
    throw invalid_argument("Index out of range");
  }
};

void my_node_map::erase(string key){
  if(this->has(key)){
    int pos = this->index[key];
    this->keys.erase(this->keys.begin()+pos);
    this->values.erase(this->values.begin()+pos);
    this->index.erase(key);
  }
  else{
    throw invalid_argument("Key does not exists");
  }
};

bool my_node_map::has(string key){
  return this->index.find(key) != this->index.end();
};

node * my_node_map::value(string key){
  if(this->has(key)){
    return this->values[this->index[key]];
  }
  throw invalid_argument("Key Error");
};


void my_node_map::insert(string key,node * val){
  if(this->has(key)){
    this->values[this->index[key]] = val;
  }
  else{
    this->keys.push_back(key);
    this->values.push_back(val);
    this->index[key] = this->keys.size()-1;
  }
};