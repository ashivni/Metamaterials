#ifndef _Includes
#define _Includes

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <iterator>
#include <vector>
#include <deque>
#include <set>
#include <array>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>
using namespace std;
#endif

#ifndef _Node
#define _Node
class coordinates_2d;
class node {

  public:
    int i, j;
    string id;
    map<string, node *> * neighbors;

    node (int, int);
    node (const node &obj);
    string get_id();
    bool bonded_to(node *);
    void add_neighbor(node *);
    void debond_from(node *);
    bool is_connected_to(node *);
    coordinates_2d coordinates() const;
    friend std::ostream& operator<< (std::ostream &, const node &);
};
#endif

#ifndef _Coordinates2D
#define _Coordinates2D
class coordinates_2d{
  public:
    vector <double> x, y;
    coordinates_2d(vector <double>, vector <double>);
    coordinates_2d();
    void append(vector<double>, vector<double>);
    void append(coordinates_2d );
    friend std::ostream& operator<< (std::ostream &, const coordinates_2d &);
};
#endif

#ifndef _Hexagon
#define _Hexagon
class node;
class bond;
class coordinates_2d;
class hexagon {
  public:
    int i, j, a;
    vector <node *> nodes;
    vector <bond *> bonds;
    hexagon(int, int, int);
    coordinates_2d coordinates();
    friend std::ostream& operator<< (std::ostream &, const hexagon &);
};
#endif

#ifndef _Hexagonal_grid
#define _Hexagonal_grid
class hexagonal_unit_cell;
class coordinates_2d;
class hexagonal_grid {
  public:
    int nx, ny, level, l0, magnification, a;
    map<int, hexagonal_unit_cell *> master_grid;
    hexagonal_grid(int, int, int, int, int);
    map<string, node *> * node_map();
    map<string, bond *> * bond_map();
    set<string> * all_cells;
    coordinates_2d node_coordinates();
};
#endif

#ifndef _Hexagonal_unit_cell
#define _Hexagon_unit_cell
class node;
class bond;
class coordinates_2d;
class hexagonal_unit_cell {
  public:
    int i, j, a;
    string id;
    vector <node *> nodes;
    vector <bond *> interior_bonds;
    vector <bond *> exterior_bonds;
    vector <hexagonal_unit_cell *> neighbors;
    hexagonal_unit_cell(int, int, int, set<string> *);
    void add_neighbor(hexagonal_unit_cell *, string);
    bool has_neighbor(hexagonal_unit_cell *);
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
class node;
class coordinates_2d;
class nodes_and_bonds;

class bond {
  public:
    node *n1, *n2;
    string id;
    string bond_type;

    bond (node *, node *, string);
    bond (const bond &obj);
    string get_id();
    vector< coordinates_2d > coordinates(int, int);
    void sever();
    nodes_and_bonds refine(int magnification);
    friend std::ostream& operator<< (std::ostream &, const bond &);

};
#endif

#ifndef _Hierarchical_grid
#define _Hierarchical_grid
class hierarchical_grid {
  public:
    int n_step, nx, ny, levels, l0, magnification, a, lx, ly, ly_ind, lx_ind;
    bool notch, damage;
    double damage_fraction, EY;
    hexagonal_grid * outline;
    map<int, map<string, node *> *> * level_nodes;
    map<int, map<string, bond *> *> * level_bonds;
    map<int, map<string, bond *> *> * level_broken_bonds;
    map<int, bool> * level_is_build;
    hierarchical_grid(int, int, int, int, int, bool, bool, double, double);
    void build_eqns();
};
#endif