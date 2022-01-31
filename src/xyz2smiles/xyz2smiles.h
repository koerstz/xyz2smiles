#ifndef XYZ2SMILES_H
#define XYZ2SMILES_H

#include <GraphMol/GraphMol.h>

class Xyz2Smiles {

    public:
        bool chargedFragments = true;
        bool ignoreChirality = false;
        bool useAtomMap = false;

        Xyz2Smiles() {};
        Xyz2Smiles(bool chargedFragments, bool ignoreChirality, bool useAtomMap)
            : chargedFragments(chargedFragments), ignoreChirality(ignoreChirality), useAtomMap(useAtomMap) {};
        
        std::string from_xyzfile(std::string filename, int charge);

};

#endif