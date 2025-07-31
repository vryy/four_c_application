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
     * Gather some arguments from the command line and store them in a struct. This should be
     * changed to a proper parser
     */
    struct CommandlineArguments
    {
      std::string input_file_name;
      std::string output_file_identifier;
      std::string restart_file_identifier;
      int restart_step = 0;

      std::shared_ptr<FourC::Core::Communication::Communicators> comms;
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void setup_parallel_output(const CommandlineArguments& arguments);

    void setup_global_problem(FourC::Core::IO::InputFile& input_file, const CommandlineArguments& arguments);

    void entrypoint_switch();

    std::vector<std::string> parse_input_output_files(int argc, char** argv, int my_rank) const;

    void parse_commandline_arguments(int argc, char** argv, CommandlineArguments& arguments) const;

    void parse_restart_definition(const std::vector<std::string>& inout, const int in_out_args,
            std::string& restart_file_identifier, const std::string& outfile_identifier,
            const int restart_group, CommandlineArguments& arguments) const;

    double walltime_in_seconds() const;

    void print_help_message() const;

    void run(const CommandlineArguments& arguments);

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
