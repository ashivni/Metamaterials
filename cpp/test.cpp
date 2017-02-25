#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>

using namespace std;

int main()
{
  vector<int> v1 = {7,5,6};
  vector<int> v2 = {0,1,2};

  sort(v2.begin(), v2.end(), [&v1](int a, int b){return v1[a] < v1[b];});
  sort(v1.begin(), v1.end());

  cout << "v1: ";
  for(auto const & a: v1){
    cout << a << "\t";
  }
  cout << endl;

  cout << "v2: ";
  for(auto const & a: v2){
    cout << a << "\t";
  }
  cout << endl;

  return 1;
}
