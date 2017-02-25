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
#include <umfpack.h>
#include "st_to_cc.hpp"
#include <cs.h>

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
    ~node();
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
    void destruct();
    ~hexagon();
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
class CPPSparse;
class hierarchical_grid {
  public:
    int n_step, nx, ny, levels, l0, magnification, a, lx, ly, ly_ind, lx_ind;
    bool notch, damage;
    double damage_fraction, EY, notch_len;
    hexagonal_grid * outline;
    map<int, map<string, node *> *> * level_nodes;
    map<int, map<string, bond *> *> * level_bonds, * level_notched_bonds;
    map<int, map<string, bond *> *> * level_broken_bonds;
    map<int, bool> * level_is_build;
    map<int, bool> * level_is_solved;
    map<int, bool> * level_not_broken;
    map<int, double> * level_curr_dn, *level_notch_len;
    map<int, vector <double> *> * level_voltage, * level_current, *level_scaled_stress, *level_stress, *level_strain;
    map<int, CPPSparse *> * level_L, * level_B, * level_V, * level_C, * level_curr;
    map<int, map <string, int> *> *level_interior_node_indices, *level_exterior_up_node_indices, *level_exterior_dn_node_indices;
    map<int, vector <string> *>  *level_interior_node_ids, *level_exterior_up_node_ids, *level_exterior_dn_node_ids;
    map<int, map<string, int> *> *level_int_vars;

    hierarchical_grid(int, int, int, int, int, bool, bool, double, double, double notch_len = 0.25);
    void build_eqns();
    void solve();
    void step();
    void simulate_fracture(int);
    bool is_broken();
    void dump(string pref = "level");
    void save();
    ~hierarchical_grid();
    void destruct();
    double bond_current(bond *, int );

};
#endif

#ifndef _Linear_system
#define _Linear_system
class linear_system{
  public:
    int *i, *j, nnz, m, n, ncc, *icc, *ccc;
    double *A, *b, *x, *acc;
    map <string, int> index_map;
    bool is_solved;
    linear_system(int, int, int, vector<int> *, vector<int> *, vector<double> *, vector<double> *);
    vector<double> solve();
    void update_mat(int, int, double);
    void update_mat(int, double);
    void update_load(int, double);
};
#endif

#ifndef _CPPSparse
#define _CPPSparse
class CPPSparse{
  public:
    cs_di * cs_t, *cs_c;
    map<string, int> sp_index_t, sp_index_c;
    CPPSparse(int, int, int, vector<int> *, vector<int> *, vector<double> *);
    CPPSparse();
    CPPSparse(const CPPSparse &);
    void CPPSprint();
    int data_index(int,int,char);
    void make_index();
    CPPSparse multiply(const CPPSparse &) const;
    void copy(const CPPSparse &);
    CPPSparse operator * (const CPPSparse &) const;
    CPPSparse add(const CPPSparse &)const;
    CPPSparse operator + (const CPPSparse &) const;
    CPPSparse subtract(const CPPSparse &)const;
    CPPSparse operator - (const CPPSparse &) const;
    CPPSparse & operator = (const CPPSparse &);
    vector <double> solve(vector <double>);
    void destruct();
    void update(int, int, double);
    double get(int, int);
    double max_data();
    ~CPPSparse();
};
cs_di * csc_to_cst(cs_di *);
#endif

#ifndef __cplusplus
#define __cplusplus
#endif

#ifndef _UtilFuncs
#define _UtilFuncs
vector<vector<double>> stress(hierarchical_grid, int);
#endif