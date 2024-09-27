/*
LICENSE: see four_c_application/LICENSE.txt
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 25/09/2024 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_FOUR_C_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_FOUR_C_APPLICATION_VARIABLES_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#ifndef SD_APP_FORWARD_COMPATIBILITY
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"
#endif

namespace Kratos
{

#ifdef SD_APP_FORWARD_COMPATIBILITY

// Variables definition
/* KRATOS_DEFINE_APPLICATION_VARIABLE(FOUR_C_APPLICATION, double, VARIABLE_NAME) */ // Example declaring scalar varibles
/* KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(FOUR_C_APPLICATION, VARIABLE_NAME) */ // Example declaring 3D varibles

#else // SD_APP_FORWARD_COMPATIBILITY

// Variables definition
/* KRATOS_DEFINE_VARIABLE(double, VARIABLE_NAME) */ // Example declaring scalar varibles
/* KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(VARIABLE_NAME) */ // Example declaring 3D varibles

#endif

} // namespace Kratos

#endif // KRATOS_FOUR_C_APPLICATION_H_INCLUDED defined
