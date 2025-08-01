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

#if !defined(KRATOS_FOUR_C_MODEL_H_INCLUDED)
#define KRATOS_FOUR_C_MODEL_H_INCLUDED

// System includes
#include <unordered_map>

// External includes
#include "4C_fem_discretization.hpp"

// Project includes
#include "includes/define.h"

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
 * Interface model to provide access to 4C discretizations
 * Model keeps a list of Discretization and perform operations on them
 */
class FourCModel
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FourCModel);

    typedef int IndexType;

    struct disc_data_t
    {
        std::shared_ptr<FourC::Core::FE::Discretization> disc;

        std::shared_ptr<FourC::Core::LinAlg::SparseOperator> mat1;
        std::shared_ptr<FourC::Core::LinAlg::SparseOperator> mat2;
        std::shared_ptr<FourC::Core::LinAlg::Vector<double> > vec1;
        std::shared_ptr<FourC::Core::LinAlg::Vector<double> > vec2;
        std::shared_ptr<FourC::Core::LinAlg::Vector<double> > vec3;
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FourCModel(const MPI_Comm& comm, const int dim);

    /// Destructor.
    virtual ~FourCModel();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /*********************************/
    /****** SETUP OPERATIONS *********/
    /*********************************/

    /// Create a discretization with name
    void CreateDiscretization(const std::string& name);

    /// Create node for a specific discretization
    void CreateNode(const std::string& dis_name, const IndexType id, const double x, const double y, const double z);

    /// Create element for a specific discretization
    void CreateElement(const std::string& dis_name, const std::string& element_type, const IndexType id, const std::vector<IndexType>& nodes);

    /// Add the discretization to the model
    void AddDiscretization(std::shared_ptr<FourC::Core::FE::Discretization> pdisc);

    /// Fill complete the discretization
    void FillComplete();

    /*********************************/
    /****** COMPUTE OPERATIONS *******/
    /*********************************/

    /// Set the (zero) state for discretization
    void SetZeroState(const std::string& dis_name, const unsigned nds,
            const std::string& state_name);

    /// Set the state for discretization
    void SetState(const std::string& dis_name, const unsigned nds,
            const std::string& state_name, std::shared_ptr<const FourC::Core::LinAlg::Vector<double> > state);

    /// Evaluate the discretization providing the parameter list
    void Evaluate(Teuchos::ParameterList& params, const std::string& dis_name);

    ///@}
    ///@name Access
    ///@{

    /// Access the discretization (public)
    std::shared_ptr<const FourC::Core::FE::Discretization> pGetDiscretization(const std::string& dis_name) const
    {
        return GetDiscretizationData(dis_name).disc;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        return "FourCModel";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
        for (auto it = mDiscretizationData.begin(); it != mDiscretizationData.end(); ++it)
        {
            it->second.disc->print(rOStream);
            rOStream << std::endl;
        }
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

    /// Access the discretization (protected)
    std::shared_ptr<FourC::Core::FE::Discretization> pGetDiscretization(const std::string& dis_name)
    {
        return GetDiscretizationData(dis_name).disc;
    }

    /// Access the discretization
    const disc_data_t GetDiscretizationData(const std::string& dis_name) const
    {
        if (mDiscretizationData.find(dis_name) == mDiscretizationData.end())
            KRATOS_ERROR << "4C Discretization " << dis_name << " does not exist";
        return mDiscretizationData.at(dis_name);
    }

    /// Access the discretization (protected)
    disc_data_t GetDiscretizationData(const std::string& dis_name)
    {
        if (mDiscretizationData.find(dis_name) == mDiscretizationData.end())
            KRATOS_ERROR << "4C Discretization " << dis_name << " does not exist";
        return mDiscretizationData.at(dis_name);
    }

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

    int mDimension;
    MPI_Comm mComm;
    std::unordered_map<std::string, disc_data_t> mDiscretizationData;

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
    FourCModel& operator=(FourCModel const& rOther);

    /// Copy constructor.
    FourCModel(FourCModel const& rOther);

    ///@}

}; // Class FourCModel

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const FourCModel& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

} // namespace Kratos

#endif // KRATOS_FOUR_C_MODEL_H_INCLUDED defined
