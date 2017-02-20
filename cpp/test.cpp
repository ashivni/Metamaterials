#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

int main()
{
  vector <int *> v;
  int * n = new int(2);
  v.push_back(n);
  n = new int(3);
  v.push_back(n);

  for(auto const & t: v){
    cout << *t << endl;
  }

  return 1;
}
