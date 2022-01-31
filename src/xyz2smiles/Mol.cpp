#include <iostream>
#include <set>
#include <map>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>

#include "Mol.h"
#include "Fileio.h"
#include "AdjMatrix.h"
#include "VecUtils.h"

using namespace boost;

typedef property <edge_weight_t, int, property<edge_index_t, int>> EdgeProperty;
typedef adjacency_list<vecS, vecS, undirectedS, no_property, EdgeProperty> Graph;

ValenceMap atomicValence = {
    {"H", {1}},
    {"B", {3,4}},
    {"C", {4}},
    {"N", {3,4}},
    {"O", {2,1,3}},
    {"F", {1}},
    {"Si", {4}},
    {"P", {5,3}},
    {"S", {6,3,2}},
    {"Cl", {1}},
    {"Ge", {4}},
    {"Br", {1}},
    {"I", {1}}
};

std::unordered_map<int, RDKit::Bond::BondType> bondMap = { 
        {1, RDKit::Bond::BondType::SINGLE},
        {2, RDKit::Bond::BondType::DOUBLE},
        {3, RDKit::Bond::BondType::TRIPLE}
    };

// Constructor functions
Mol::Mol(std::vector<Atom> atoms, int charge, bool chargedFragments, bool ignoreChirality, bool useAtomMap) 
    : atoms(atoms), charge(charge), chargedFragments(chargedFragments), ignoreChirality(ignoreChirality), useAtomMap(useAtomMap)
{
    SetAC();
    SetBO();
}

Mol Mol::FromXYZ(std::string FileName, int charge, bool chargedFragments, bool ignoreChirality, bool useAtomMap ) 
{   
    return Mol(ReadXYZ(FileName), charge, chargedFragments, ignoreChirality, useAtomMap);
}

// Member functions
void Mol::SetAC(double covalentFactor) {
        AC = AdjacencyMatrixDist(atoms, covalentFactor);
        BO = AC;
}

void Mol::SetBO(double covalentFactor) {
    // 1) Find all possible combinations of allowed valences.
    Eigen::VectorXi ACvalence = AC.colwise().sum();
    std::vector<std::vector<int> > allValences;
    
    for (int i = 0; i < atoms.size(); i++) {        
        std::vector<int> possibleValence;
            
        // valence can not be smaller than degree.
        for (int val : atomicValence[atoms[i].symbol]) {  
            if ( val >= ACvalence[i] ) {
                possibleValence.push_back(val);
            } 
        }
        allValences.push_back(possibleValence);
    }
    
    // 2: Validate the valence possibilities, and identify the best BO.
    std::vector<std::vector<int> > allValenceProduct = IterProduct(allValences); 
    for (std::vector<int> valences: allValenceProduct) {
        SetUAaDU(valences, ACvalence); 

        // All atoms fully saturated.
        if ( UA.empty() && CheckBO(valences) ) { 
            BO = AC;
            break;
        }

        ListofPairs UAPairs = GetUAPairs();
        ComputeBO(valences, UAPairs);
        if ( CheckBO(valences) )  {
            break;
        }
    }
}

bool Mol::CheckBO(std::vector<int> valences) {
    if ( BOValencesToLarge(valences) ) {
        return false;
    }
    bool checkSum = ( (BO - AC).sum() == SumVector(DU) );
    bool checkCharge = ChargeOK(valences);
    if (checkSum && checkCharge) {
        return true;
    }
    return false;
}


void Mol::SetUAaDU(std::vector<int> maxValenceList, Eigen::VectorXi valanceList) {
    // Reset UA and DU.
    UA.clear();
    DU.clear();
    for (int i = 0; i < maxValenceList.size(); i++) {
        if ( valanceList[i] < maxValenceList[i]) {
            UA.push_back(i);
            DU.push_back(maxValenceList[i] - valanceList[i]);
        }
    }
}

bool Mol::BOValencesToLarge(std::vector<int> valences) {
    Eigen::VectorXi numberOfBonds = BO.colwise().sum();
    for ( int i = 0; i < numberOfBonds.size(); i++ ) {
        if (numberOfBonds[i] > valences[i]) {
            return true;
        }
    }
    return false;
}

bool Mol::ChargeOK(std::vector<int> valences) {
    int Q = 0; // total charge

    if ( ! chargedFragments ) {
        return (charge == Q);
    }

    Eigen::VectorXi BOvalences = BO.colwise().sum();
    for (int i = 0; i < atoms.size(); i++) {
        int q = GetAtomicCharge(i,  BOvalences[i]);
        Q += q;

        if ( atoms[i].symbol == "C" ) {
            // count single bonds to C.
            int cSingleBonds = 0;
            for (int bondBO: BO(Eigen::all, i) ) {if (bondBO == 1) {
                cSingleBonds += 1;}
            } 
            if ( cSingleBonds == 2 && BOvalences[i] == 2 ) {Q += 1;}
            if ( cSingleBonds == 3 && Q + 1 < charge) {Q += 2;}
        }
        
    }   
    return (charge == Q);
}

ListofPairs Mol::GetUAPairs() {
    ListofPairs UAPairs;
    ListofPairs bonds = {};
    for (int i = 0; i < UA.size(); i++) {
        for (int j =  i + 1; j < UA.size(); j++) {
            if (AC(UA[i], UA[j]) == 1 ) {
                if (UA[i] < UA[j]) {
                    bonds.push_back(std::make_pair(UA[i], UA[j]));
                } else {
                    bonds.push_back(std::make_pair(UA[j], UA[i]));
                }
            }
        }
    }
    if ( bonds.empty() ) {
        return bonds;
    }
        
    // build graph and run maximum_weighted_matching, and extract pairs.
    std::set<int> s;
    for (std::pair<int,int> bond: bonds) {
        s.insert(bond.first);
        s.insert(bond.second);
    }

    std::map<int,int> molToGraphID;
    std::map<int,int> graphTomolID;
    int graphID = 0;
    for (int molAtomId: s) {
        molToGraphID[molAtomId] = graphID;
        graphTomolID[graphID] = molAtomId;
        graphID++;
    }

    Graph g(s.size());
    for (std::pair<int,int> bond: bonds) {            
        boost::add_edge(
            molToGraphID[bond.first], molToGraphID[bond.second], EdgeProperty(1), g);
    }

    // Maximum weighted matching
    std::vector< graph_traits<Graph>::vertex_descriptor > mate1(s.size());
    maximum_weighted_matching(g, &mate1[0]);

    graph_traits< Graph >::vertex_iterator vi, vi_end;
    for ( boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi ) {
        if ( mate1[*vi] != graph_traits< Graph >::null_vertex() && *vi < mate1[*vi] ) {
            UAPairs.push_back(std::make_pair(graphTomolID[*vi], graphTomolID[mate1[*vi]]));
        }
    }
    return UAPairs;
}

void Mol::ComputeBO(std::vector<int> valences, std::vector<std::pair<int, int> > &UAPairs) {
    
    BO = AC; // Reset BO.
    std::vector<int> DUsave = {};

    while (DUsave != DU) {
        for (std::pair<int, int> p: UAPairs) { 
            BO(p.first, p.second) += 1;
            BO(p.second, p.first) += 1;
        }

        Eigen::VectorXi BOvalence = BO.colwise().sum();
        
        // copy DU vector.
        DUsave = {};
        for (int i: DU) {
            DUsave.push_back(i);
        }
        UAPairs = GetUAPairs();
    }
}

int Mol::GetAtomicCharge(int atomID, int BOvalence) {

    RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();

    if (atoms[atomID].symbol == "H") {
        return 1 - BOvalence;
    } else if (atoms[atomID].symbol == "B") {
        return 3 - BOvalence;
    } else if (atoms[atomID].symbol == "P" && BOvalence == 5) {
        return 0;
    } else if (atoms[atomID].symbol == "S" && BOvalence == 6) {
        return 0;
    } else {
        return tbl->getNouterElecs(atoms[atomID].symbol) - 8 + BOvalence;
    }
}

std::string Mol::GetSmiles() {
    
    if (atoms[0].symbol == "H") {
        std::cout << "Don't start with H" << std::endl;
        exit(1);
    }

    RDKit::RWMol m = RDKit::RWMol( *RDKit::SmilesToMol(atoms[0].symbol));
    // Add Atoms
    for (int i = 1; i < atoms.size(); i++) {
        RDKit::Atom a = RDKit::Atom(atoms[i].symbol);
        int idx = m.addAtom(&a);
    }
    
    // Add Bonds
    for (int i = 0; i < atoms.size(); i++) {
        for (int j = i + 1; j < atoms.size(); j++) {
            if ( BO(i,j) != 0) {
                m.addBond(i, j, bondMap[BO(i,j)]);
            }
        }
    }
    
    // TODO: Set atomic charge/radical if needed.
    if (chargedFragments) {
        SetRDKitMolAtomicCharges(m);   
    } else {
        SetRDKitMolRadicals(m);
    }

    if (! ignoreChirality ) {
        RDKit::MolOps::sanitizeMol(m);
    }


    // Add Conformer - with xyz coords
    RDKit::ROMol *mol = new RDKit::ROMol(m);
    RDKit::Conformer conf = RDKit::Conformer(atoms.size());
    for (int i = 0; i < atoms.size(); i++) {
        conf.setAtomPos(i, RDGeom::Point3D(atoms[i].x, atoms[i].y, atoms[i].z));
    }
    mol->addConformer(&conf);

    if (! ignoreChirality  ) {
        RDKit::MolOps::assignStereochemistryFrom3D(*mol, -1, true);
    }

    if (useAtomMap) {
        for (int i = 0; i < atoms.size(); i++ ) {
            mol->getAtomWithIdx(i)->setAtomMapNum(i + 1); 
        }
    }

    // Make Mol Canonical
    RDKit::ROMol *outm = RDKit::SmilesToMol(RDKit::MolToSmiles(*mol));

    // check charge.
    return RDKit::MolToSmiles(*outm);
}

void Mol::SetRDKitMolAtomicCharges(RDKit::RWMol &mol) {

    Eigen::VectorXi  BOvalence = BO.colwise().sum();
    int q = 0;

    for (int i = 0; i < atoms.size(); i++) {
        int charge = GetAtomicCharge(i, BOvalence[i]);
        q += charge;
        if (atoms[i].symbol == "C") {
            int cSingleBonds = 0;
            for (int bondBO: BO(Eigen::all, i) ) {
                if (bondBO == 1) {
                    cSingleBonds += 1;
                }
            } 

            if (cSingleBonds == 2 && BOvalence[i] == 2) {
                q += 1;
                charge = 0;
            }

            if ( cSingleBonds == 3 && q + 1 < charge ) {
                q += 2;
                charge = 1;
            }
        }
                
        if ( std::abs(charge) > 0 ) {
            mol.getAtomWithIdx(i)->setFormalCharge(charge);
        }
    }
}

void Mol::SetRDKitMolRadicals(RDKit::RWMol &mol) {

    Eigen::VectorXi  BOvalence = BO.colwise().sum();

    for (int i = 0; i < atoms.size(); i++) {
        int charge = GetAtomicCharge(i, BOvalence[i]);

        if ( std::abs(charge) > 0 ) {
            mol.getAtomWithIdx(i)->setNumRadicalElectrons(charge);
        }
    }
}   