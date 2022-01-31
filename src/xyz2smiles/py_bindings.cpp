#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Mol.h"
#include "xyz2smiles.h"

namespace py = pybind11;

PYBIND11_MODULE(xyz2smiles, m)
{
    py::class_<Xyz2Smiles>(m, "xyz2smiles")
        .def(py::init< >())
        .def(py::init<bool, bool, bool>())
        .def_readwrite("chargedFragments", &Xyz2Smiles::chargedFragments)
        .def_readwrite("ignoreChirality", &Xyz2Smiles::ignoreChirality)
        .def_readwrite("useAtomMap", &Xyz2Smiles::useAtomMap)
        .def("from_xyzfile", &Xyz2Smiles::from_xyzfile, "a function", py::arg("file"), py::arg("charge"), py::return_value_policy::automatic);
}
