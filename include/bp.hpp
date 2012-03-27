#ifndef BP_HPP
#define BP_HPP
#include <map>
#include <vector>
#include "math.hpp"
#include <list>
#include "mesh.hpp"
class BP{
public:
  /**@brief msg[to][from][label]=value*/
  std::vector< std::map<int ,std::map<int,real_t> > >msg;
  void maxProd();
  std::list<int> schedule;
  BP(Mesh & m);
  std::vector<Plane> plane;
  float distw;
  real_t bpdistance(Trig & t , Plane & p);
  void update_distance(std::vector<real_t > & dist,
                       std::vector<real_t > & cdf,
                      std::vector<Plane> & plane,
                      Mesh & m, int ll);
};
#endif
