#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py=pybind11;

#include "mmh.h"

PYBIND11_MODULE(mmh, m)
{
  py::class_<McMillanHe>(m, "McMillanHe")
    .def(py::init<>())
    .def("diffuse", &McMillanHe::diffuse)
    .def("get_acc", &McMillanHe::get_acc)
    .def("get_lbox", &McMillanHe::get_lbox)
    .def("wfval", &McMillanHe::wfval)
    ;
}
