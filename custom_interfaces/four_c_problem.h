/*
LICENSE: see four_c_application/LICENSE.txt
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 26/09/2024 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_FOUR_C_PROBLEM_H_INCLUDED)
#define KRATOS_FOUR_C_PROBLEM_H_INCLUDED

// System includes
#include <filesystem>

// External includes
#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"

// Project includes
#include "includes/define.h"
#include "custom_interfaces/four_c_model.h"

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
 * Wrapper around 4C GlobalProblem to invoke 4C analysis
 */
class FourCProblem
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FourCProblem);

    typedef int IndexType;

    /**
     * Structure to hold command line arguments.
     */
    struct CommandlineArguments
    {
      bool help = false;
      int n_groups = 1;
      bool parameters = false;
      std::vector<int> group_layout = {};
      FourC::Core::Communication::NestedParallelismType nptype =
          FourC::Core::Communication::NestedParallelismType::no_nested_parallelism;
      int diffgroup = -1;
      int restart = 0;
      std::string restart_file_identifier = "";
      std::vector<int> restart_per_group = {};
      std::vector<std::string> restart_identifier_per_group = {};
      bool interactive = false;
      std::vector<std::pair<std::filesystem::path, std::string>> io_pairs;
      std::filesystem::path input_file_name = "";
      std::string output_file_identifier = "";
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FourCProblem(int argc, char** argv);

    /// Destructor.
    virtual ~FourCProblem();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Run the analysis defined in input file
    void Run();

    ///@}
    ///@name Access
    ///@{

    std::vector<std::string> GetDiscretizationNames() const;

    FourCModel::ConstPointer pGetModel() const {return mpModel;}

    FourCModel::Pointer pGetModel() {return mpModel;}

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        return "FourCProblem";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const;

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

    FourC::Global::Problem* mpProblem;
    FourCModel::Pointer mpModel;
    CommandlineArguments mArguments;
    std::shared_ptr<FourC::Core::Communication::Communicators> mpCommunicators;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void setup_parallel_output(const CommandlineArguments& arguments,
            const FourC::Core::Communication::Communicators& communicators);

    void setup_global_problem(FourC::Core::IO::InputFile& input_file, const CommandlineArguments& arguments,
            const FourC::Core::Communication::Communicators& communicators);

    void entrypoint_switch();

    CommandlineArguments parse_command_line(int argc, char** argv) const;

    double walltime_in_seconds() const;

    void print_help_message() const;

    void run(CommandlineArguments& cli_args, const FourC::Core::Communication::Communicators& communicators);

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
    FourCProblem& operator=(FourCProblem const& rOther);

    /// Copy constructor.
    FourCProblem(FourCProblem const& rOther);

    ///@}

}; // Class FourCProblem

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const FourCProblem& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

} // namespace Kratos

#endif // KRATOS_FOUR_C_PROBLEM_H_INCLUDED defined
