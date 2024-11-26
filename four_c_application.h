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
//  +   25/09/2024: create four_c_application.h

#if !defined(KRATOS_FOUR_C_APPLICATION_H_INCLUDED)
#define KRATOS_FOUR_C_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#define FOUR_C_APP_DEFINE_ELEMENT_ALL_GEOMETRIES(element_type) \
    const element_type m##element_type##2D3N; \
    const element_type m##element_type##2D6N; \
    const element_type m##element_type##2D4N; \
    const element_type m##element_type##2D8N; \
    const element_type m##element_type##2D9N; \
    const element_type m##element_type##3D4N; \
    const element_type m##element_type##3D10N; \
    const element_type m##element_type##3D8N; \
    const element_type m##element_type##3D20N; \
    const element_type m##element_type##3D27N; \
    const element_type m##element_type##3D6N; \
    const element_type m##element_type##3D15N;

#define FOUR_C_APP_DEFINE_CONDITION_ALL_GEOMETRIES(condition_type) \
    const condition_type m##condition_type##2D2N; \
    const condition_type m##condition_type##2D3N; \
    const condition_type m##condition_type##3D3N; \
    const condition_type m##condition_type##3D6N; \
    const condition_type m##condition_type##3D4N; \
    const condition_type m##condition_type##3D8N; \
    const condition_type m##condition_type##3D9N;

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Enum's
///@{

///@}
///@name Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * Application provides access to 4C interface
 */
class KRATOS_API(FOUR_C_APPLICATION) KratosFourCApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosMultiphaseApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosFourCApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosFourCApplication();

    /// Destructor.
    ~KratosFourCApplication() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Description of your application";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "KratosFourCApplication contains the following components:" << std::endl;
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosFourCApplication& operator=(KratosFourCApplication const& rOther);

    /// Copy constructor.
    KratosFourCApplication(KratosFourCApplication const& rOther);

    ///@}

}; // Class KratosFourCApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos

#endif // KRATOS_FOUR_C_APPLICATION_H_INCLUDED defined
