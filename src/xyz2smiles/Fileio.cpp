#include <fstream>
#include <iostream>
#include <sstream>

#include "Fileio.h"
#include "Mol.h"

std::vector<Atom> ReadXYZ(std::string FileName) {

    std::vector<Atom> atoms;

    std::ifstream infile(FileName);
    if (!infile) {
        std::cout << "!!!! ERROR opening xyz file: ";
        std::cout << FileName << " !!!!" << std::endl;
        exit(-1);
    }

    std::string line;
    std::string cell;

    // Read number of atoms
    getline( infile, line );
    int natoms = std::stoi(line);
    
    // skip header line
    getline( infile, line );

    // read coordinates
    for (int i = 0; i < natoms; i++) {
        std::vector<std::string> data;

        std::getline(infile, line);

        std::stringstream lineStream(line);
        while (std::getline(lineStream, cell, ' ')) {
            if (cell.empty()) {
                continue;
            }
            data.push_back(cell);
        }

        atoms.push_back(
            Atom {data[0], std::stod(data[1]), std::stod(data[2]), std::stod(data[3])}
        );
    }
    infile.close();
    return atoms;
}