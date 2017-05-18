#ifndef MINCUT_HPP
#define MINCUT_HPP
#include "GCoptimization.h"
#include "mesh.hpp"
extern float dataCostW;
extern float smoothW;
extern float distw;
extern float distCenterW;
extern real_t saliency_weight;
void runMincut(Mesh & m);
extern int MC_ITER;
#endif
