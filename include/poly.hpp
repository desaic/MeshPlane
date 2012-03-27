#ifndef POLY_HPP
#define POLY_HPP
#include "mesh.hpp"
#include <vector>
//mesh is labeled with plane assignment
class Poly{
public:
	Mesh * m;
	Poly(Mesh * _m );
	std::vector<std::vector<int> >p;
};
#endif
