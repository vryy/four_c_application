/*
LICENSE: see four_c_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: 25 Sep 2024 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes


// External includes
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>


// Project includes
#include "mpi_python/mpi_python.h"
#include "custom_python/add_four_c_model_to_python.h"
#include "custom_interfaces/four_c_model.h"
#include "custom_interfaces/four_c_problem.h"

namespace Kratos
{

namespace Python
{

static FourCModel::Pointer FourCModel_init(PythonMPIComm& rPythonMPIComm, const int dim)
{
    return FourCModel::Pointer(new FourCModel(rPythonMPIComm.GetMPIComm(), dim));
}

void FourCModel_CreateElement(FourCModel& rDummy, const std::string& dis_name, const std::string& element_type,
        const int id, const boost::python::list& list_nodes)
{
    typedef FourCModel::IndexType IndexType;
    std::vector<IndexType> nodes;
    typedef boost::python::stl_input_iterator<int> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& i,
        std::make_pair(iterator_value_type(list_nodes), iterator_value_type() ) )
    {
        nodes.push_back(i);
    }
    rDummy.CreateElement(dis_name, element_type, id, nodes);
}

static FourCProblem::Pointer FourCProblem_init(const boost::python::list& arguments)
{
    int argc = 0;
    std::vector<char*> argv;

    typedef boost::python::stl_input_iterator<std::string> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& str,
        std::make_pair(iterator_value_type(arguments), iterator_value_type() ) )
    {
        ++argc;
        argv.push_back(const_cast<char*>(str.c_str()));
    }

    return FourCProblem::Pointer(new FourCProblem(argc, argv.data()));
}

void FourCApplication_AddFourCModelToPython()
{
    namespace bp = boost::python;

    bp::class_<FourCModel, FourCModel::Pointer, boost::noncopyable>
    ("FourCModel", bp::no_init)
    .def("__init__", bp::make_constructor(&FourCModel_init))
    .def("CreateDiscretization", &FourCModel::CreateDiscretization)
    .def("CreateNode", &FourCModel::CreateNode)
    .def("CreateElement", &FourCModel_CreateElement)
    .def("FillCompleteDiscretization", &FourCModel::FillCompleteDiscretization)
    // .def("EvaluateSystem", &FourCModel::EvaluateSystem)
    .def(bp::self_ns::str(bp::self))
    ;

    FourCModel::Pointer(FourCProblem::*pointer_to_pGetModel)() = &FourCProblem::pGetModel;

    bp::class_<FourCProblem, FourCProblem::Pointer, boost::noncopyable>
    ("FourCProblem", bp::no_init)
    .def("__init__", bp::make_constructor(&FourCProblem_init))
    .def("ReadInputFile", &FourCProblem::ReadInputFile)
    .def("GetModel", pointer_to_pGetModel)
    .def(bp::self_ns::str(bp::self))
    ;
}

} // namespace Python.

} // namespace Kratos.
