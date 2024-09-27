/*
LICENSE: see four_c_application/LICENSE.txt
*/
//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: 25/09/2024 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes

// External includes
#if defined(KRATOS_PYTHON)

// Project includes
#include "includes/define_python.h"
#include "four_c_application.h"
#include "four_c_application_variables.h"

namespace Kratos
{

namespace Python
{

PYBIND11_MODULE(KratosFourCApplication, m)
{
    namespace py = pybind11;

    py::class_<KratosFourCApplication, KratosFourCApplication::Pointer, KratosApplication >
    (m, "KratosFourCApplication")
    .def(py::init<>())
    ;

    /* KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, VARIABLE_NAME ) */ // Example registering scalar varibles to python
    /* KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( m, VARIABLE_NAME ) */ // Example registering 3D varibles to python

}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON
