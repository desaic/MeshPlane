#include "poly.hpp"

Poly::Poly(Mesh * _m):m(_m)
{
	std::vector<bool > visited;
	visited.resize(m->t.size());
	unsigned int idx = 0;
	m->adjlist();
	while(idx<visited.size()){
		if(visited[idx]){
			idx++;
			continue;
		}

	}
}
