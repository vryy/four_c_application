/*
LICENSE: see four_c_application/LICENSE.txt
*/

//
//   Project Name:        KratosFourCApplication
//   Last modified by:    $Author: hbui $
//   Date:                $Date: 25 Sep 2024 $
//
//

// System includes


// External includes
#include <boost/python.hpp>


// Project includes
#include "mpi_python/mpi_python.h"
#include "python/python_utils.h"
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
    PythonUtils::Unpack<IndexType>(list_nodes, nodes);
    rDummy.CreateElement(dis_name, element_type, id, nodes);
}

////////////////////////////////////////////////////////////////////////////////

static FourCProblem::Pointer FourCProblem_init(const boost::python::list& arguments)
{
    std::vector<std::string> args;
    args.push_back(""); // this is for the first argument ("4C") as CLI requires it
    PythonUtils::UnpackAndAppend<std::string>(arguments, args);

    int argc = args.size();

    std::vector<char*> argv;
    for (int i = 0; i < argc; ++i)
        argv.push_back(const_cast<char*>(args[i].c_str()));

    return FourCProblem::Pointer(new FourCProblem(argc, argv.data()));
}

boost::python::list FourCProblem_GetDiscretizationNames(FourCProblem& rDummy)
{
    boost::python::list a;
    auto names = rDummy.GetDiscretizationNames();
    for (auto it = names.begin(); it != names.end(); ++it)
        a.append(*it);
    return std::move(a);
}

////////////////////////////////////////////////////////////////////////////////

void FourCApplication_AddFourCModelToPython()
{
    namespace bp = boost::python;

    void(FourCModel::*pointer_to_Evaluate)(Teuchos::ParameterList&, const std::string&) = &FourCModel::Evaluate;
    FourCModel::disc_data_t(FourCModel::*pointer_to_GetDiscretizationData)(const std::string&) const = &FourCModel::GetDiscretizationData;

    bp::class_<FourCModel::disc_data_t>
    ("FourCDiscData", bp::no_init)
    .def_readonly("disc", &FourCModel::disc_data_t::disc)
    .def_readonly("mat1", &FourCModel::disc_data_t::mat1)
    .def_readonly("mat2", &FourCModel::disc_data_t::mat2)
    .def_readonly("vec1", &FourCModel::disc_data_t::vec1)
    .def_readonly("vec2", &FourCModel::disc_data_t::vec2)
    .def_readonly("vec3", &FourCModel::disc_data_t::vec3)
    ;

    bp::class_<FourCModel, FourCModel::Pointer, boost::noncopyable>
    ("FourCModel", bp::no_init)
    .def("__init__", bp::make_constructor(&FourCModel_init))
    .def("CreateDiscretization", &FourCModel::CreateDiscretization)
    .def("CreateNode", &FourCModel::CreateNode)
    .def("CreateElement", &FourCModel_CreateElement)
    .def("FillComplete", &FourCModel::FillComplete)
    .def("SetZeroState", &FourCModel::SetZeroState)
    .def("Evaluate", pointer_to_Evaluate)
    .def("GetDiscretizationData", pointer_to_GetDiscretizationData)
    .def(bp::self_ns::str(bp::self))
    ;

    FourCModel::Pointer(FourCProblem::*pointer_to_pGetModel)() = &FourCProblem::pGetModel;

    bp::class_<FourCProblem, FourCProblem::Pointer, boost::noncopyable>
    ("FourCProblem", bp::no_init)
    .def("__init__", bp::make_constructor(&FourCProblem_init))
    .def("GetModel", pointer_to_pGetModel)
    .def("GetDiscretizationNames", &FourCProblem_GetDiscretizationNames)
    .def("Run", &FourCProblem::Run)
    .def(bp::self_ns::str(bp::self))
    ;
}

} // namespace Python.

} // namespace Kratos.
