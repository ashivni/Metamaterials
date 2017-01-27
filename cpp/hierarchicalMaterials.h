#ifndef _Includes
#define _Includes

#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <iterator>
#include <vector>
#include <deque>
#include <set>
#include <array>
#include <cmath>
using namespace std;
#endif

#ifndef _Node
#define _Node
class node {

  public:
    int i, j, n_neighbors;
    string id;
    map<string, node *> neighbors;

    node (int, int);
    node (const node &obj);
    string get_id();
    bool bonded_to(node *);
    void add_neighbor(node *);
    void debond_from(node *);
    bool is_connected_to(node *);
    array <float, 2> coordinates();
};
#endif

#ifndef _Hexagon
#define _Hexagon
class node;
class bond;
class hexagon {
  public:
    int i, j, a;
    vector <node *> nodes;
    vector <bond *> bonds;
    hexagon(int, int, int);
};
#endif

#ifndef _Nodes_And_Bonds
#define _Nodes_And_Bonds
class node;
class bond;
class nodes_and_bonds{
  public:
    map<string, node *> nodes;
    map<string, bond *> bonds;
};
#endif

#ifndef _Bond
#define _Bond
class bond {
  node *n1, *n2;
  string id;
  string bond_type;
  public:
    bond (node *, node *, string);
    bond (const bond &obj);
    string get_id();
    array< vector< array<float, 2> >,2> coordinates(int, int);
    void sever();
    nodes_and_bonds refine(int magnification);
};
#endif