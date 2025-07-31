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
//Change log:
//  +   25/09/2024: create four_c_application.cpp

// System includes

// External includes
#include "4C_utils_singleton_owner.hpp"

// Project includes
#include "version.h"
#include "four_c_application.h"
#include "four_c_application_variables.h"
/* Including as needed the geometries */
// #include "geometries/triangle_2d_3.h"
// #include "geometries/triangle_3d_3.h"
// #include "geometries/triangle_2d_6.h"
// #include "geometries/triangle_3d_6.h"
// #include "geometries/quadrilateral_2d_4.h"
// #include "geometries/quadrilateral_3d_4.h"
// #include "geometries/quadrilateral_2d_8.h"
// #include "geometries/quadrilateral_3d_8.h"
// #include "geometries/quadrilateral_2d_9.h"
// #include "geometries/quadrilateral_3d_9.h"
// #include "geometries/tetrahedra_3d_4.h"
// #include "geometries/tetrahedra_3d_10.h"
// #include "geometries/prism_3d_6.h"
// #include "geometries/prism_3d_15.h"
// #include "geometries/hexahedra_3d_8.h"
// #include "geometries/hexahedra_3d_20.h"
// #include "geometries/hexahedra_3d_27.h"

#ifdef SD_APP_FORWARD_COMPATIBILITY
#define FOUR_C_APP_CREATE_ELEMENT(element_type, geometry_type, number_of_nodes) \
    element_type( 0, Element::GeometryType::Pointer( new geometry_type <Node>( Element::GeometryType::PointsArrayType( number_of_nodes ) ) ) )
#define FOUR_C_APP_CREATE_CONDITION(condition_type, geometry_type, number_of_nodes) \
    condition_type( 0, Condition::GeometryType::Pointer( new geometry_type <Node>( Condition::GeometryType::PointsArrayType( number_of_nodes ) ) ) )
#else
#define FOUR_C_APP_CREATE_ELEMENT(element_type, geometry_type, number_of_nodes) \
    element_type( 0, Element::GeometryType::Pointer( new geometry_type <Node<3> >( Element::GeometryType::PointsArrayType( number_of_nodes, Node<3>() ) ) ) )
#define FOUR_C_APP_CREATE_CONDITION(condition_type, geometry_type, number_of_nodes) \
    condition_type( 0, Condition::GeometryType::Pointer( new geometry_type <Node<3> >( Condition::GeometryType::PointsArrayType( number_of_nodes, Node<3>() ) ) ) )
#endif

#define FOUR_C_APP_CREATE_ELEMENT_ALL_GEOMETRIES(element_type) \
    FOUR_C_APP_CREATE_ELEMENT(m##element_type##2D3N, Triangle2D3, 3), \
    FOUR_C_APP_CREATE_ELEMENT(m##element_type##2D6N, Triangle2D6, 6), \
    FOUR_C_APP_CREATE_ELEMENT(m##element_type##2D4N, Quadrilateral2D4, 4), \
    FOUR_C_APP_CREATE_ELEMENT(m##element_type##2D8N, Quadrilateral2D8, 8), \
    FOUR_C_APP_CREATE_ELEMENT(m##element_type##2D9N, Quadrilateral2D9, 9), \
    FOUR_C_APP_CREATE_ELEMENT(m##element_type##3D4N, Tetrahedra3D4, 4), \
    FOUR_C_APP_CREATE_ELEMENT(m##element_type##3D10N, Tetrahedra3D10, 10), \
    FOUR_C_APP_CREATE_ELEMENT(m##element_type##3D8N, Hexahedra3D8, 8), \
    FOUR_C_APP_CREATE_ELEMENT(m##element_type##3D20N, Hexahedra3D20, 20), \
    FOUR_C_APP_CREATE_ELEMENT(m##element_type##3D27N, Hexahedra3D27, 27), \
    FOUR_C_APP_CREATE_ELEMENT(m##element_type##3D6N, Prism3D6, 6), \
    FOUR_C_APP_CREATE_ELEMENT(m##element_type##3D15N, Prism3D15, 15) \

#define FOUR_C_APP_CREATE_CONDITION_ALL_GEOMETRIES(condition_type) \
    FOUR_C_APP_CREATE_CONDITION(m##condition_type##2D2N, Line2D2, 2), \
    FOUR_C_APP_CREATE_CONDITION(m##condition_type##2D3N, Line2D3, 3), \
    FOUR_C_APP_CREATE_CONDITION(m##condition_type##3D3N, Triangle3D3, 3), \
    FOUR_C_APP_CREATE_CONDITION(m##condition_type##3D6N, Triangle3D6, 6), \
    FOUR_C_APP_CREATE_CONDITION(m##condition_type##3D4N, Quadrilateral3D4, 4), \
    FOUR_C_APP_CREATE_CONDITION(m##condition_type##3D8N, Quadrilateral3D8, 8), \
    FOUR_C_APP_CREATE_CONDITION(m##condition_type##3D9N, Quadrilateral3D9, 9) \

#define FOUR_C_APP_REGISTER_ELEMENT_ALL_GEOMETRIES(element_type) \
    KRATOS_REGISTER_ELEMENT( #element_type"2D3N", m##element_type##2D3N ) \
    KRATOS_REGISTER_ELEMENT( #element_type"2D6N", m##element_type##2D6N ) \
    KRATOS_REGISTER_ELEMENT( #element_type"2D4N", m##element_type##2D4N ) \
    KRATOS_REGISTER_ELEMENT( #element_type"2D8N", m##element_type##2D8N ) \
    KRATOS_REGISTER_ELEMENT( #element_type"2D9N", m##element_type##2D9N ) \
    KRATOS_REGISTER_ELEMENT( #element_type"3D4N", m##element_type##3D4N ) \
    KRATOS_REGISTER_ELEMENT( #element_type"3D10N", m##element_type##3D10N ) \
    KRATOS_REGISTER_ELEMENT( #element_type"3D8N", m##element_type##3D8N ) \
    KRATOS_REGISTER_ELEMENT( #element_type"3D20N", m##element_type##3D20N ) \
    KRATOS_REGISTER_ELEMENT( #element_type"3D27N", m##element_type##3D27N ) \
    KRATOS_REGISTER_ELEMENT( #element_type"3D6N", m##element_type##3D6N ) \
    KRATOS_REGISTER_ELEMENT( #element_type"3D15N", m##element_type##3D15N ) \

#define FOUR_C_APP_REGISTER_CONDITION_ALL_GEOMETRIES(condition_type) \
    KRATOS_REGISTER_CONDITION( #condition_type"2D2N", m##condition_type##2D2N ) \
    KRATOS_REGISTER_CONDITION( #condition_type"2D3N", m##condition_type##2D3N ) \
    KRATOS_REGISTER_CONDITION( #condition_type"3D3N", m##condition_type##3D3N ) \
    KRATOS_REGISTER_CONDITION( #condition_type"3D6N", m##condition_type##3D6N ) \
    KRATOS_REGISTER_CONDITION( #condition_type"3D4N", m##condition_type##3D4N ) \
    KRATOS_REGISTER_CONDITION( #condition_type"3D8N", m##condition_type##3D8N ) \
    KRATOS_REGISTER_CONDITION( #condition_type"3D9N", m##condition_type##3D9N ) \

namespace Kratos
{

KratosFourCApplication::KratosFourCApplication()
#ifdef SD_APP_FORWARD_COMPATIBILITY
    : KratosApplication("FourCApplication")
#else
    : KratosApplication()
#endif
      // , FOUR_C_APP_CREATE_ELEMENT_ALL_GEOMETRIES( ElementType ) // Example creating elements for all geometries
      // , FOUR_C_APP_CREATE_CONDITION_ALL_GEOMETRIES( ConditionType ) // Example creating conditions for all geometries
{
    FourC::Core::Utils::SingletonOwnerRegistry::initialize();
    // here we manually initialize and finalize the SingletonOwnerRegistry
    // instead of using SingletonOwnerRegistry::ScopeGuard as in 4C.
    // Using that's only applicable if everything is called within a main() function
}

KratosFourCApplication::~KratosFourCApplication()
{
    FourC::Core::Utils::SingletonOwnerRegistry::finalize();
    std::cout << "KratosFourCApplication is unloaded" << std::endl;
}

void KratosFourCApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosFourCApplication..."
              << " Revision: " << FOUR_C_APPLICATION_GIT_REV
              << " Tag: " << FOUR_C_APPLICATION_GIT_TAG
              << " Branch: " << FOUR_C_APPLICATION_GIT_BRANCH
              << std::endl;

    // KRATOS_REGISTER_VARIABLE( VARIABLE_NAME ) // Example registering scalar varibles
    // KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( VARIABLE_NAME ) // Example registering 3D varibles

    // FOUR_C_APP_REGISTER_ELEMENT_ALL_GEOMETRIES( ElementType ) // Example registering elements for all geometries
    // FOUR_C_APP_REGISTER_CONDITION_ALL_GEOMETRIES( ConditionType ) // Example registering conditions for all geometries
}

} // namespace Kratos
