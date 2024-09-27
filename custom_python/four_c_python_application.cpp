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
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "four_c_application.h"
#include "four_c_application_variables.h"
#include "custom_python/add_four_c_model_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;
BOOST_PYTHON_MODULE(KratosFourCApplication)
{

    class_<KratosFourCApplication, KratosFourCApplication::Pointer, bases<KratosApplication>, boost::noncopyable>
    ("KratosFourCApplication");

    // KRATOS_REGISTER_IN_PYTHON_VARIABLE( VARIABLE_NAME ) // Example registering scalar varibles to python
    // KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( VARIABLE_NAME ) // Example registering 3D varibles to python

    FourCApplication_AddFourCModelToPython();
}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON
