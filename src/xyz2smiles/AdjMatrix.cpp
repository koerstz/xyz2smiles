#include <vector>
#include <Eigen/Dense>
#include <GraphMol/PeriodicTable.h>

#include "AdjMatrix.h"
#include "Mol.h"

double AtomDistance(Atom A, Atom B) {
    double dist = std::sqrt(
        (A.x - B.x)*(A.x - B.x) +
        (A.y - B.y)*(A.y - B.y) +
        (A.z - B.z)*(A.z - B.z)
    );
    return dist;
}

Eigen::MatrixXd DistanceMatrix(std::vector<Atom> atoms) {

    Eigen::MatrixXd distM(atoms.size(), atoms.size());

    for (int i = 0; i < atoms.size(); i++) {
        for (int j = i; j < atoms.size(); j++) {
            distM(j, i) = AtomDistance(atoms[i], atoms[j]);
            distM(i, j) = distM(j, i);
        }
    }
    return distM;
}

Eigen::MatrixXi AdjacencyMatrixDist(const std::vector<Atom> atoms, double covalentFactor) {
    
    RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();   

    Eigen::MatrixXi AC = Eigen::MatrixXi::Zero(atoms.size(), atoms.size());
    Eigen::MatrixXd distM = DistanceMatrix(atoms);

    for (int i = 0; i < atoms.size(); i++) {
        for (int j = i + 1; j < atoms.size(); j++) {
            double atomIRcov = tbl->getRcovalent(atoms[i].symbol) * covalentFactor;
            double atomJRcov = tbl->getRcovalent(atoms[j].symbol) * covalentFactor;

            if (distM(i, j) <= atomIRcov + atomJRcov) {
                AC(i, j) = 1;
                AC(j, i) = 1;
            }
        }
    }
    return AC;
}                 