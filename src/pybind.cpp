/*
 * pybind.cpp
 *
 *  Created on: Feb 11, 2012
 *      Author: Moises Osorio [WCoder]
 */

#include <boost/python.hpp>

#include "momet.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(momet)
{
    // Create the Python type object for our extension class and define __init__ function.
    class_<Momet>("Momet")
        .def("errorRatio", &Momet::errorRatio)
        .def("generationalDistance", &Momet::genDistance)
        .def("spacing", &Momet::spacing)
        .def("coverage", &Momet::coverage)
        .def("additiveEpsilon", &Momet::addEpsilonIndicator)
        .def("multiplicativeEpsilon", &Momet::multEpsilonIndicator)
    ;
}
