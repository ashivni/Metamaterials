#include "hierarchicalMaterials.h"

coordinates_2d::coordinates_2d(vector <double> x, vector <double> y){
  this->append(x, y);
}

coordinates_2d::coordinates_2d(){

}

void coordinates_2d::append(coordinates_2d c){
  this->append(c.x, c.y);
}

void coordinates_2d::append(vector <double> x, vector <double> y){
  for(auto const& value: x) {
    this->x.push_back(value);
  }
  for(auto const& value: y) {
    this->y.push_back(value);
  }
}

std::ostream& operator<< (std::ostream &out, const coordinates_2d &c){
  for(int i = 0; i < c.x.size(); i++){
    out << c.x[i] << ",\t" << c.y[i] << endl;
  }

  return out;
}