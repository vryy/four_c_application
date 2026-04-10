/*
LICENSE: see four_c_application/LICENSE.txt
*/

//
//   Project Name:        KratosFourCApplication
//   Last modified by:    $Author: hbui $
//   Date:                $Date: 10 Apr 2026 $
//
//

// System includes


// External includes
#include <boost/python.hpp>

#include "4C_fem_discretization.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_linalg_fevector.hpp"

// Project includes
#include "includes/define_python.h"
#include "python/python_utils.h"
#include "custom_python/add_four_c_classes_to_python.h"


namespace Kratos
{

namespace Python
{

template< class T>
std::string PrintFourCObject(const T& rObject)
{
    std::stringstream ss;
    rObject.print(ss);
    return ss.str();
}

void FourCApplication_AddFourCClassesToPython()
{
    namespace bp = boost::python;

    typedef FourC::Core::FE::Discretization         _4C_Discretization;
    typedef FourC::Core::LinAlg::SparseOperator     _4C_SparseOperator;
    typedef FourC::Core::LinAlg::SparseMatrix       _4C_SparseMatrix;
    typedef FourC::Core::LinAlg::Vector<double>     _4C_Vector;
    typedef FourC::Core::LinAlg::FEVector<double>   _4C_FEVector;

    bp::class_<_4C_Discretization, std::shared_ptr<_4C_Discretization>, boost::noncopyable>
    ("_4C_Discretization", bp::no_init)
    .def("__str__", PrintFourCObject<_4C_Discretization>)
    ;

    bp::class_<_4C_SparseOperator, std::shared_ptr<_4C_SparseOperator>, boost::noncopyable>
    ("_4C_SparseOperator", bp::no_init)
    ;

    bp::class_<_4C_SparseMatrix, std::shared_ptr<_4C_SparseMatrix>, bp::bases<_4C_SparseOperator>, boost::noncopyable>
    ("_4C_SparseMatrix", bp::no_init)
    .def("num_global_rows", &_4C_SparseMatrix::num_global_rows)
    .def("num_global_cols", &_4C_SparseMatrix::num_global_cols)
    .def("num_global_entries", &_4C_SparseMatrix::num_global_entries)
    .def("num_global_nonzeros", &_4C_SparseMatrix::num_global_nonzeros)
    .def("__str__", PrintFourCObject<_4C_SparseMatrix>)
    ;

    bp::class_<_4C_Vector, std::shared_ptr<_4C_Vector>, boost::noncopyable>
    ("_4C_Vector", bp::no_init)
    .def("local_length", &_4C_Vector::local_length)
    .def("global_length", &_4C_Vector::global_length)
    .def("__str__", PrintFourCObject<_4C_Vector>)
    ;

    bp::class_<_4C_FEVector, std::shared_ptr<_4C_FEVector>, boost::noncopyable>
    ("_4C_FEVector", bp::no_init)
    .def("__str__", PrintFourCObject<_4C_FEVector>)
    ;

}

} // namespace Python.

} // namespace Kratos.
