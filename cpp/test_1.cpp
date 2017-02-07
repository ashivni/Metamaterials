#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <map>
#include <stdexcept>

using namespace std;

template <class K, class V>
class my_map {
  private:
    vector <K> keys;
    vector <V> values;
    map <K, int> index;
  public:
    int size();
    void insert(K, V);
    bool has(K);
    V value(K);
};

template <class K, class V>
int my_map<K,V>::size(){
  return this->keys.size();
};

template <class K, class V>
bool my_map<K,V>::has(K key){
  return this->index.find(key) != this->index.end();
};

template <class K, class V>
V my_map<K,V>::value(K key){
  if(this->has(key)){
    return this->values[this->index[key]];
  }
  throw invalid_argument("Key Error");
};


template <class K, class V>
void my_map<K,V>::insert(K key,V val){
  if(this->has(key)){
    this->values[this->index[key]] = val;
  }
  else{
    this->keys.push_back(key);
    this->values.push_back(val);
    this->index[key] = this->keys.size()-1;
  }
};

int main(){
  int Nmax = 100;
  int i = 0;

  cout << "My map initilization" << endl;
  i = 0;
  my_map <int, int> mm;
  while(i < Nmax){
    mm.insert(i,i+1);
    i+=1;
  }
  cout << mm.value(Nmax-2) << endl;
  cout << mm.size() << endl;

  cout << "Map initilization" << endl;
  i = 0;
  map <int, int> m;
  while(i < Nmax){
    m[i] = i+1;
    i+=1;
  }
  cout << m[Nmax-2] << endl;

  cout << "Vector initilization" << endl;
  i = 0;
  vector <int> v;
  while(i < Nmax){
    v.push_back(i+1);
    i+=1;
  }
  cout << v[Nmax-2] << endl;


  return 1;
}