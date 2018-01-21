#include <iostream>
#include <bitset>
#include <limits>
int main(){   
  int a = 124156;
  std::bitset<std::numeric_limits<int>::digits> bset;
  std::cout<<a<<std::endl;
  // convert a to bitset bset
  bset = std::bitset<std::numeric_limits<int>::digits>(a); 
  std::cout<<bset<<std::endl;
  // print 3rd bit of bitset 
  std::cout<<bset[3]<<std::endl;
  // flip 3rd bit of bitset
  bset.flip(3) ; 
  // print 3rd bit of bitset again
  std::cout<<bset[3]<<std::endl;
  std::cout<<bset<<std::endl;
  // convert bitset back to integer a
  a=static_cast<int>(bset.to_ulong());
  std::cout<<a<<std::endl;
}
