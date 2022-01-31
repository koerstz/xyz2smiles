#ifndef ADJM_H
#define ADJM_H

#include <Eigen/Dense>
#include "Mol.h"

double AtomDistance(Atom A, Atom B);

Eigen::MatrixXd DistanceMatrix(std::vector<Atom> atoms);

Eigen::MatrixXi AdjacencyMatrixDist(std::vector<Atom> atoms, double covalentFactor=1.3);

#endif