#include "hierarchicalMaterials.h"

// Constructors
CPPSparse::CPPSparse(){
  this->cs_t = NULL;
  this->cs_c = NULL;
  /*
  this->cs_t = (cs_di *) malloc(sizeof(cs_di));
  this->cs_c = (cs_di *) malloc(sizeof(cs_di));

  this->cs_t->i = NULL;
  this->cs_t->p = NULL;
  this->cs_t->x = NULL;
  this->cs_t->m = 0;
  this->cs_t->n = 0;
  this->cs_t->nzmax = 0;
  this->cs_t->nz = 0;

  this->cs_c->i = NULL;
  this->cs_c->p = NULL;
  this->cs_c->x = NULL;
  this->cs_c->m = 0;
  this->cs_c->n = 0;
  this->cs_c->nzmax = 0;
  this->cs_c->nz = 0;
  */
}

CPPSparse::CPPSparse(int n, int m, int nnz, vector<int> *i, vector<int> *j, vector<double> *A){
  // A sparse matrix. Triplet format.
  //Validate
  if(i->size() != nnz) throw invalid_argument("i.size() must be nnz.");
  if(j->size() != nnz) throw invalid_argument("j.size() must be nnz.");
  if(A->size() != nnz) throw invalid_argument("A.size() must be nnz.");

  this->cs_t = (cs_di *) malloc(sizeof(cs_di));

  // initialize
  this->cs_t->p = (int *) malloc(sizeof(int)*nnz);
  this->cs_t->i = (int *) malloc(sizeof(int)*nnz);
  this->cs_t->x = (double *) malloc(sizeof(double)*nnz);
  this->cs_t->nzmax = nnz;
  this->cs_t->m = m;
  this->cs_t->n = n;
  this->cs_t->nz = nnz;

  // copy
  for(int k = 0; k < nnz; k++){
    this->cs_t->i[k] = (*i)[k];
    this->cs_t->p[k] = (*j)[k];
    this->cs_t->x[k] = (*A)[k];
  }

  // make column compressed
  this->cs_c = (cs_di *) malloc(sizeof(cs_di));
  this->cs_c->p = (int *) malloc(sizeof(int)*(n+1));
  this->cs_c->i = (int *) malloc(sizeof(int)*nnz);
  this->cs_c->nzmax = nnz;
  this->cs_c->m = m;
  this->cs_c->n = n;
  this->cs_c->nz = -1;

  st_to_cc_index (this->cs_c->nzmax, this->cs_t->i, this->cs_t->p, this->cs_c->nzmax, this->cs_c->n, this->cs_c->i, this->cs_c->p);
  this->cs_c->x = st_to_cc_values ( this->cs_c->nzmax, this->cs_t->i, this->cs_t->p, this->cs_t->x, this->cs_c->nzmax, this->cs_c->n, this->cs_c->i, this->cs_c->p);


  // make index
  this->make_index();
}

cs_di * csc_to_cst(cs_di * csc){
  cs_di * cst = (cs_di *)malloc(sizeof(cs_di));
  cst->m = csc->m;
  cst->n = csc->n;
  cst->nzmax = csc->nzmax;
  cst->nz = csc->nzmax;
  cst->p = (int *) malloc(sizeof(int)*cst->nzmax);
  cst->i = (int *) malloc(sizeof(int)*cst->nzmax);
  cst->x = (double *) malloc(sizeof(double)*cst->nzmax);

  int count = 0;
  for(int col = 0; col < csc->n; col++){
    for(int k = csc->p[col]; k < csc->p[col+1]; k++){
      cst->i[k] = csc->i[k];
      cst->p[k] = col;
      cst->x[k] = csc->x[k];
    }
  }

  return cst;
}

void CPPSparse::make_index(){
  this->sp_index_c.erase(this->sp_index_c.begin(), this->sp_index_c.end());
  for(int col = 0; col < this->cs_c->m; col++){
    for(int k = this->cs_c->p[col]; k < this->cs_c->p[col+1]; k++){
      int row = this->cs_c->i[k];
      ostringstream stringStream;
      stringStream << row << "_" << col;
      this->sp_index_c[stringStream.str()] = k;
    }
  }

  this->sp_index_t.erase(this->sp_index_t.begin(), this->sp_index_t.end());
  for(int k = 0; k < this->cs_t->nzmax; k++){
      int row = this->cs_t->i[k];
      int col = this->cs_t->p[k];
      ostringstream stringStream;
      stringStream << row << "_" << col;
      this->sp_index_t[stringStream.str()] = k;
  }
}

void CPPSparse::copy(const CPPSparse &B){
  this->destruct();
  if (B.cs_t != NULL){
    this->cs_t = (cs_di *) malloc(sizeof(cs_di));
    this->cs_t->i = (int *) malloc(sizeof(int)*B.cs_t->nzmax);
    this->cs_t->p = (int *) malloc(sizeof(int)*B.cs_t->nzmax);
    this->cs_t->x = (double *) malloc(sizeof(double)*B.cs_t->nzmax);

    for(int k = 0; k < B.cs_t->nzmax; k++){
      this->cs_t->p[k] = B.cs_t->p[k];
      this->cs_t->i[k] = B.cs_t->i[k];
      this->cs_t->x[k] = B.cs_t->x[k];
    }
    this->cs_t-> nzmax = B.cs_t->nzmax;
    this->cs_t-> n = B.cs_t->n;
    this->cs_t-> m = B.cs_t->m;
    this->cs_t-> nz = B.cs_t->nz;
  }

  if(B.cs_c != NULL){
    this->cs_c = (cs_di *) malloc(sizeof(cs_di));
    this->cs_c->i = (int *) malloc(sizeof(int)*B.cs_c->nzmax);
    this->cs_c->p = (int *) malloc(sizeof(int)*B.cs_c->n+1);
    this->cs_c->x = (double *) malloc(sizeof(double)*B.cs_c->nzmax);

    for(int k = 0; k < B.cs_c->nzmax; k++){
      this->cs_c->i[k] = B.cs_c->i[k];
      this->cs_c->x[k] = B.cs_c->x[k];
    }
    for(int k = 0; k < B.cs_c->n+1; k++){
      this->cs_c->p[k] = B.cs_c->p[k];
    }
    this->cs_c-> nzmax = B.cs_c->nzmax;
    this->cs_c-> n = B.cs_c->n;
    this->cs_c-> m = B.cs_c->m;
    this->cs_c-> nz = B.cs_c->nz;
  }

  this->sp_index_t.erase(this->sp_index_t.begin(), this->sp_index_t.end());
  for(auto const & t: B.sp_index_t){
    this->sp_index_t[t.first] = t.second;
  }

  this->sp_index_c.erase(this->sp_index_c.begin(), this->sp_index_c.end());
  for(auto const & t: B.sp_index_c){
    this->sp_index_c[t.first] = t.second;
  }
}

CPPSparse & CPPSparse::operator = (const CPPSparse &B){
  if(&B != this){
    this->destruct();
    this->copy(B);
  }
  return *this;
}
CPPSparse::CPPSparse(const CPPSparse &B){
  if (B.cs_t != NULL){
    this->cs_t = (cs_di *) malloc(sizeof(cs_di));
    this->cs_t->i = (int *) malloc(sizeof(int)*B.cs_t->nzmax);
    this->cs_t->p = (int *) malloc(sizeof(int)*B.cs_t->nzmax);
    this->cs_t->x = (double *) malloc(sizeof(double)*B.cs_t->nzmax);

    for(int k = 0; k < B.cs_t->nzmax; k++){
      this->cs_t->p[k] = B.cs_t->p[k];
      this->cs_t->i[k] = B.cs_t->i[k];
      this->cs_t->x[k] = B.cs_t->x[k];
    }
    this->cs_t-> nzmax = B.cs_t->nzmax;
    this->cs_t-> n = B.cs_t->n;
    this->cs_t-> m = B.cs_t->m;
    this->cs_t-> nz = B.cs_t->nz;
  }

  if(B.cs_c != NULL){
    this->cs_c = (cs_di *) malloc(sizeof(cs_di));
    this->cs_c->i = (int *) malloc(sizeof(int)*B.cs_c->nzmax);
    this->cs_c->p = (int *) malloc(sizeof(int)*B.cs_c->n+1);
    this->cs_c->x = (double *) malloc(sizeof(double)*B.cs_c->nzmax);

    for(int k = 0; k < B.cs_c->nzmax; k++){
      this->cs_c->i[k] = B.cs_c->i[k];
      this->cs_c->x[k] = B.cs_c->x[k];
    }
    for(int k = 0; k < B.cs_c->n; k++){
      this->cs_c->p[k] = B.cs_c->p[k];
    }
    this->cs_c-> nzmax = B.cs_c->nzmax;
    this->cs_c-> n = B.cs_c->n;
    this->cs_c-> m = B.cs_c->m;
    this->cs_c-> nz = B.cs_c->nz;
  }

  this->sp_index_t.erase(this->sp_index_t.begin(), this->sp_index_t.end());
  for(auto const & t: B.sp_index_t){
    this->sp_index_t[t.first] = t.second;
  }

  this->sp_index_c.erase(this->sp_index_c.begin(), this->sp_index_c.end());
  for(auto const & t: B.sp_index_c){
    this->sp_index_c[t.first] = t.second;
  }
}

void CPPSparse::destruct(){

  if(this->cs_t != NULL){
    if(this->cs_t->i != NULL) free(this->cs_t->i);
    if(this->cs_t->p != NULL) free(this->cs_t->p);
    if(this->cs_t->x != NULL) free(this->cs_t->x);
    free(cs_t);
  }
  if(this->cs_c != NULL){
    if(this->cs_c->i != NULL) free(this->cs_c->i);
    if(this->cs_c->p != NULL) free(this->cs_c->p);
    if(this->cs_c->x != NULL) free(this->cs_c->x);
    free(cs_c);
  }
  this->sp_index_t.erase(this->sp_index_t.begin(), this->sp_index_t.end());
  this->sp_index_c.erase(this->sp_index_c.begin(), this->sp_index_c.end());
}

void CPPSparse::update(int row, int col, double x){
  if (this->data_index(row, col, 'C') != -1){
    this->cs_c->x[data_index(row, col, 'C')] = x;
  }

  if (this->data_index(row, col, 'T') != -1){
    this->cs_t->x[data_index(row, col, 'T')] = x;
  }
}

double CPPSparse::max_data(){
  double max = -numeric_limits<double>::infinity();
  if(this->cs_c != NULL && this->cs_c->nzmax > 0){
    for(int k = 0; k < this->cs_c->nzmax; k++){
      max = this->cs_c->x[k] > max ? this->cs_c->x[k] : max;
    }
    return max;
  }
  throw invalid_argument("Null structure has no max.");
}

double CPPSparse::get(int row, int col){
  if (row < 0 || row >= this->cs_c->n || col < 0 || col >= this->cs_c->m){
    cout << row << ", " << col << ", " << this->cs_c->m << ", " << this->cs_c->n << endl;
    throw invalid_argument("Index out of range in CPPSparse::get.");
  }

  if (this->data_index(row, col, 'C') != -1){
    return this->cs_c->x[data_index(row, col, 'C')];
  }
  return 0;
}


CPPSparse CPPSparse::multiply(const CPPSparse &B)const {
  if (this->cs_c->n != B.cs_c->m){
    throw invalid_argument("Inconsistent dimensions in multiply.");
  }
  CPPSparse m;
  m.cs_c = cs_multiply(this->cs_c, B.cs_c);
  m.cs_t = csc_to_cst(m.cs_c);

  m.make_index();
  return m;
}

CPPSparse CPPSparse::operator * (const CPPSparse &B) const {
  return this->multiply(B);
}

CPPSparse CPPSparse::add(const CPPSparse &B)const {
  if (this->cs_c->n != B.cs_c->n || this->cs_c->m != B.cs_c->m){
    throw invalid_argument("Inconsistent dimensions in add.");
  }  CPPSparse m;
  m.cs_c = cs_add(this->cs_c, B.cs_c, 1.0, 1.0);
  m.cs_t = csc_to_cst(m.cs_c);
  m.make_index();
  return m;
}

CPPSparse CPPSparse::operator + (const CPPSparse &B) const {
  return this->add(B);
}

CPPSparse CPPSparse::subtract(const CPPSparse &B)const {
  if (this->cs_c->n != B.cs_c->n || this->cs_c->m != B.cs_c->m){
    throw invalid_argument("Inconsistent dimensions in add.");
  }
  CPPSparse m;
  m.cs_c = cs_add(this->cs_c, B.cs_c, 1.0, -1.0);
  m.cs_t = csc_to_cst(m.cs_c);
  m.make_index();
  return m;
}

CPPSparse CPPSparse::operator - (const CPPSparse &B) const {
  return this->subtract(B);
}

vector <double> CPPSparse::solve(vector <double> b){
  // Solve Ax = b
  double *null = ( double * ) NULL;
  int status;
  void *Symbolic, *Numeric;
  double *sol = (double *) malloc(sizeof(double)*this->cs_c->n);

  status = umfpack_di_symbolic ( this->cs_c->m, this->cs_c->n, this->cs_c->p, this->cs_c->i, this->cs_c->x, &Symbolic, null, null );
  status = umfpack_di_numeric (this->cs_c->p, this->cs_c->i, this->cs_c->x, Symbolic, &Numeric, null, null );
  umfpack_di_free_symbolic ( &Symbolic );
  status = umfpack_di_solve ( UMFPACK_A, this->cs_c->p, this->cs_c->i, this->cs_c->x, sol, b.data(), Numeric, null, null );
  umfpack_di_free_numeric ( &Numeric );

  vector <double> x;
  for(int k=0; k < this->cs_c->n; k++){
    x.push_back(sol[k]);
  }
  free(sol);
  return x;
}


int CPPSparse::data_index(int row, int col, char m_type){
  // Return the sparse index k for the element at row, col
  // i.e return k such that A[row,col] = cs_
  ostringstream stringStream;
  stringStream << row << "_" << col;
  if(m_type == 'C'){
    if(this->sp_index_c.find(stringStream.str()) != this->sp_index_c.end()){
      return this->sp_index_c[stringStream.str()];
    }
    else{
      return -1;
    }
  }
  else if(m_type == 'T'){
    if(this->sp_index_t.find(stringStream.str()) != this->sp_index_t.end()){
      return this->sp_index_t[stringStream.str()];
    }
    else{
      return -1;
    }
  }
  else{
    throw invalid_argument("m_type should be C (column-compressed) or T (triplet)");
  }
}

CPPSparse::~CPPSparse(){
  this->destruct();
}

void CPPSparse::CPPSprint(){
  cout << "Triplet: " << endl;
  cout << "i:\t";
  for(int k = 0; k < this->cs_t->nzmax; k++){
    cout << this->cs_t->i[k] << "\t";
  }
  cout << endl;

  cout << "j:\t";
  for(int k = 0; k < this->cs_t->nzmax; k++){
    cout << this->cs_t->p[k] << "\t";
  }
  cout << endl;

  cout << "x:\t";
  for(int k = 0; k < this->cs_t->nzmax; k++){
    cout << this->cs_t->x[k] << "\t";
  }
  cout << endl << endl;

  cout << "CCS: " << endl;
  cout << "i:\t";
  for(int k = 0; k < this->cs_c->nzmax; k++){
    cout << this->cs_c->i[k] << "\t";
  }
  cout << endl;

  cout << "p:\t";
  for(int k = 0; k < this->cs_c->n+1; k++){
    cout << this->cs_c->p[k] << "\t";
  }
  cout << endl;

  cout << "x:\t";
  for(int k = 0; k < this->cs_c->nzmax; k++){
    cout << this->cs_c->x[k] << "\t";
  }
  cout << endl;
}