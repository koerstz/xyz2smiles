#ifndef MOL_H
#define MOL_H

#include <vector>
#include <string>
#include <unordered_map>
#include <Eigen/Dense>
#include <GraphMol/GraphMol.h>


typedef std::unordered_map<std::string, std::vector<int> >  ValenceMap;
typedef std::vector<std::pair<int,int> > ListofPairs;

struct Atom {
    std::string symbol;
    double x, y, z;
};


class Mol 
{

public:
    // Class constructors
    Mol() {};
    static Mol FromXYZ(
        std::string fileName, 
        int charge,
        bool chargedFragments = true,
        bool ignoreChirality = true,
        bool useAtomMap = false
    );

    // member functions
    void SetAC(double covalentFactor = 1.3);
    void SetBO(double covalentFactor = 1.3);

    std::string GetSmiles();

    Eigen::MatrixXi GetBO() {return BO;};
    Eigen::MatrixXi GetAC() {return AC;};

private:
    std::vector<Atom> atoms;

    int charge;
    bool chargedFragments;
    bool ignoreChirality;
    bool useAtomMap;

    Eigen::MatrixXi AC;
    Eigen::MatrixXi BO;

    std::vector<int> UA;
    std::vector<int> DU;

    Mol(
        std::vector<Atom> atoms, 
        int charge, 
        bool chargedFragments = true,
        bool ignoreChirality = true, 
        bool useAtomMap = false
    );

    void SetUAaDU(std::vector<int> maxValenceList, Eigen::VectorXi valanceList);
    bool CheckBO(std::vector<int> valences);
    bool ChargeOK(std::vector<int> valences);
    bool BOValencesToLarge(std::vector<int> valences);
    void ComputeBO(std::vector<int> valences, std::vector<std::pair<int, int> > &UAPairs);
    int GetAtomicCharge(int atomID, int BOvalence);
    ListofPairs GetUAPairs();

    void SetRDKitMolAtomicCharges(RDKit::RWMol &mol);
    void SetRDKitMolRadicals(RDKit::RWMol &mol);
};


#endif