/*
 * pybind.cpp
 *
 *  Created on: Feb 11, 2012
 *      Author: Moises Osorio [WCoder]
 */

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/indexing_suite.hpp>

#include "momet.h"

using namespace boost::python;
using namespace std;

BOOST_PYTHON_MODULE(momet)
{
    class_<Momet>("Momet")
        .def("errorRatio", &Momet::errorRatio)
        .def("spacing", &Momet::spacing)
        .def("coverage", &Momet::coverage)
        .def("additiveEpsilon", &Momet::addEpsilonIndicator)
        .def("multiplicativeEpsilon", &Momet::multEpsilonIndicator)
        .def("hypervolume", &Momet::hypervolume)
        .def("generationalDistance", &Momet::genDistance)
        .def("generationalDistanceP", &Momet::genDistanceP)
        .def("invertedGenerationalDistance", &Momet::invertedGenDistance)
        .def("invertedGenerationalDistanceP", &Momet::invertedGenDistanceP)
        .def("deltaP", &Momet::deltaP)
    ;
    
    class_<vector<double> >("dList")
            .def(vector_indexing_suite<vector<double> >());
    class_<vector<vector<double> > >("ddList")
                .def(vector_indexing_suite<vector<vector<double> > >());
}
