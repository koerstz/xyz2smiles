#ifndef FILEIO_H
#define FILEIO_H

#include <string>
#include <vector>

#include "Mol.h"

std::vector<Atom> ReadXYZ(std::string FileName);

#endif