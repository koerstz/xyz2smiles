#include "xyz2smiles.h"
#include "Mol.h"

std::string Xyz2Smiles::from_xyzfile(std::string filename, int charge) {
    Mol m = Mol::FromXYZ(filename, charge, chargedFragments, ignoreChirality, useAtomMap);
    return m.GetSmiles();
}