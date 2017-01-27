#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <map>

using namespace std;

class node {
  public:
    int n;
    node(int n){
      this->n = n;
    }
};

class nodes {
  public:
   map<int, node *> n;
};

int main()
{

  nodes ns;
  node n(12);
  ns.n[0] = &n;

  for(auto const& x : ns.n){
    cout << x.first << ": " << x.second->n << endl;
  }

  return 1;
}
