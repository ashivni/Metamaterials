#include "hierarchicalMaterials.h"

linear_system::linear_system(int m, int n, int nnz, vector<int> *i, vector<int> *j, vector<double> *A, vector <double> *b){

  //Validate
  if(i->size() != nnz) throw invalid_argument("i.size() must be nnz.");
  if(j->size() != nnz) throw invalid_argument("j.size() must be nnz.");
  if(A->size() != nnz) throw invalid_argument("A.size() must be nnz.");
  if(b->size() != m) throw invalid_argument("b.size() must be m.");

  // initialize
  this->i = new int[nnz];
  this->j = new int[nnz];
  this->A = new double[nnz];
  this->b = new double[m];
  this->x = new double[n];
  this->nnz = nnz;
  this->m = m;
  this->n = n;

  // copy
  for(int k = 0; k < nnz; k++){
    this->i[k] = (*i)[k];
    this->j[k] = (*j)[k];
    this->A[k] = (*A)[k];
  }

  for(int k = 0; k < m; k++){
    this->b[k] = (*b)[k];
  }

  // create column-compressed representation
  this->ncc = st_to_cc_size ( nnz, this->i, this->j);
  this->icc = new int[this->ncc];
  this->ccc = new int[this->n+1];
  st_to_cc_index ( this->nnz, this->i, this->j, this->ncc, this->n, this->icc, this->ccc );
  this->acc = st_to_cc_values ( this->nnz, this->i, this->j, this->A, this->ncc, this->n, this->icc, this->ccc);
  //st_print ( m, n, nnz, this->i, this->j, this->A, "  The matrix in ST format:" );
  //cc_print ( m, n, this->ncc, this->icc, this->ccc, this->acc, "  CC Matrix:" );

  // create map of index to positions in A
  for(int k = 0; k < nnz; k++){
    ostringstream stringStream;
    stringStream <<  this->i[k] << "_" << this->j[k];
    this->index_map[stringStream.str()] = k;
  }

  this->is_solved = false;
}

void linear_system::update_mat(int i, int j, double v){
  ostringstream stringStream;
  stringStream <<  i << "_" << j;
  this->A[this->index_map[stringStream.str()]] = v;
  this->is_solved = false;
}

void linear_system::update_mat(int k, double v){
  this->A[k] = v;
  this->is_solved = false;
}

void linear_system::update_load(int k, double v){
  this->b[k] = v;
  this->is_solved = false;
}

vector <double> linear_system::solve(){
  if(! this->is_solved){
    double *null = ( double * ) NULL;
    int status;
    void *Symbolic, *Numeric;
    status = umfpack_di_symbolic ( this->m, this->n, this->ccc, this->icc, this->acc, &Symbolic, null, null );
    status = umfpack_di_numeric (this->ccc, this->icc, this->acc, Symbolic, &Numeric, null, null );
    umfpack_di_free_symbolic ( &Symbolic );
    status = umfpack_di_solve ( UMFPACK_A, this->ccc, this->icc, this->acc, this->x, this->b, Numeric, null, null );
    umfpack_di_free_numeric ( &Numeric );
    this->is_solved = true;
  }
  vector <double> x;
  for(int k=0; k < this->n; k++){
    x.push_back(this->x[k]);
  }
  return x;
}
