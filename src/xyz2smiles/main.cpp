#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <iostream>

#include "boost/program_options.hpp" 
#include "xyz2smiles.h"

namespace po = boost::program_options;

int main(int argc, char** argv) {

    bool ignoreChirality = false;
    bool noChargedFragments = false;
    bool atomMap = false;

    po::options_description desc("Allowed options"); 
    desc.add_options() 
        ("help,h", ": Print help message") 
        ("filename", po::value<std::string>()->required(), ": Input file")
        ("charge,c", po::value<int>()->default_value(0), ": Total charge of the system")
        ("ignore-chiral", po::bool_switch(&ignoreChirality), ": Ignore chiral centers")
        ("no-charged-fragments", po::bool_switch(&noChargedFragments), ": Allow radicals to be made")
        ("atom-map", po::bool_switch(&atomMap), ": Adds atom map");

    po::positional_options_description positionalOptions; 
    positionalOptions.add("filename", 1); 

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv)
                    .options(desc)
                    .positional(positionalOptions)
                    .run(), 
                  vm);
       
        if ( vm.count("help") ) {
            std::cout << "usage: ./xyz2mol [options] molecule.xyz \n" << std::endl;
            std::cout << desc << std::endl; 
            return 1; 
        }

        po::notify(vm);

    } catch (const std::exception& e) {
        std::cerr << "Unhandled Exception reached the top of main: " 
                  << e.what() << ", application will now exit" << std::endl; 
        return 1;
    }


    // Run program":"bool_switch(&chargedFragments), "
    std::string file = vm["filename"].as<std::string>();
    int charge = vm["charge"].as<int>();
    bool chargedFragments = ! noChargedFragments; // switch : hack becaouse of po::bool_switch
 
    // std::cout << "Charge: " << charge << std::endl;
    // std::cout << "Charge Fragments: " << true << " " <<chargedFragments << std::endl;
    // std::cout << "ignore ignoreChiral: " << true << " " << ignoreChirality   << std::endl;
    // std::cout << "atom Map: " << true << " " << atomMap << std::endl;


    Xyz2Smiles x2s(chargedFragments, ignoreChirality, atomMap);
    std::cout << file << ": " << x2s.from_xyzfile(file, charge) << std::endl;
    
   return 0;
}