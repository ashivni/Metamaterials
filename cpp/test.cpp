#include <iostream>
#include <vector>
#include <fstream>
#include <cs.h>
using namespace std;

int main()
{
  cs_di *c;
  int nnz = 5, n = 3, m = 3;
  c->p = new int[nnz];
  c->i = new int[nnz];
  c->x = new double[nnz];
  c->nzmax = nnz;
  c->m = m;
  c->n = n;
  c->nz = nnz;

  return 1;
}
