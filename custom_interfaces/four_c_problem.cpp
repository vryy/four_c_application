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

// System includes

// External includes
#include "4C_comm_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_global_data_read.hpp"
#include "4C_global_legacy_module.hpp"
#include "4C_io_inputreader.hpp"
#include "4C_io_pstream.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Kokkos_Core.hpp>

// Project includes
#include "four_c_problem.h"


namespace Kratos
{

FourCProblem::FourCProblem(int argc, char** argv)
{
    // MPI_Init(&argc, &argv);
    Kokkos::ScopeGuard kokkos_guard(argc, argv); // DO NOT KNOW WHAT IS THIS

    using namespace FourC;

    Teuchos::RCP<Core::Communication::Communicators> communicators =
        Core::Communication::create_comm(std::vector<std::string>(argv, argv + argc));
    Global::Problem::instance()->set_communicators(communicators);
    // Teuchos::RCP<Epetra_Comm> lcomm = communicators->local_comm();
    // Teuchos::RCP<Epetra_Comm> gcomm = communicators->global_comm();
    // int ngroups = communicators->num_groups();

    ///
    mpProblem = FourC::Global::Problem::instance();
}

FourCProblem::~FourCProblem()
{
    // MPI_Finalize();
}

void FourCProblem::PrintData(std::ostream& rOStream) const
{
    rOStream << "Problem name: " << mpProblem->problem_name() << std::endl;
    rOStream << "Number of fields: " << mpProblem->num_fields() << std::endl;

    const std::vector<std::string> dis_names = mpProblem->get_dis_names();
    for (auto it = dis_names.begin(); it != dis_names.end(); ++it)
    {
        auto pdisc = mpProblem->get_dis(*it);
        pdisc->print(rOStream);
    }
}

void FourCProblem::ReadInputFile(const std::string& inputfile_name,
        const std::string& outputfile_kenner, const std::string& restartfile_kenner)
{
    using namespace FourC;

    Teuchos::RCP<Epetra_Comm> lcomm = mpProblem->get_communicators()->local_comm();
    Teuchos::RCP<Epetra_Comm> gcomm = mpProblem->get_communicators()->global_comm();
    int group = mpProblem->get_communicators()->group_id();
    Core::Communication::NestedParallelismType npType = mpProblem->get_communicators()->np_type();

    // and now the actual reading
    Core::IO::DatFileReader reader(inputfile_name, lcomm);

    Global::read_parameter(*mpProblem, reader);

    setup_parallel_output(outputfile_kenner, lcomm, group);

    // create control file for output and read restart data if required
    mpProblem->open_control_file(*lcomm, inputfile_name, outputfile_kenner, restartfile_kenner);

    // input of materials
    Global::read_materials(*mpProblem, reader);

    // input of materials
    Global::read_contact_constitutive_laws(*mpProblem, reader);

    // input of materials of cloned fields (if needed)
    Global::read_cloning_material_map(*mpProblem, reader);

    {
        Core::UTILS::FunctionManager function_manager;
        global_legacy_module_callbacks().AttachFunctionDefinitions(function_manager);
        function_manager.read_input(reader);
        mpProblem->set_function_manager(std::move(function_manager));
    }

    // input of particles
    Global::read_particles(*mpProblem, reader);

    switch (npType)
    {
        case Core::Communication::NestedParallelismType::no_nested_parallelism:
        case Core::Communication::NestedParallelismType::every_group_read_dat_file:
        case Core::Communication::NestedParallelismType::separate_dat_files:
            // input of fields
            Global::read_fields(*mpProblem, reader);

            // read result tests
            Global::read_result(*mpProblem, reader);

            // read all types of geometry related conditions (e.g. boundary conditions)
            // Also read time and space functions and local coord systems
            Global::read_conditions(*mpProblem, reader);

            // read all knot information for isogeometric analysis
            // and add it to the (derived) nurbs discretization
            Global::read_knots(*mpProblem, reader);
        break;
        default:
            KRATOS_ERROR << "nptype (nested parallelity type) not recognized";
        break;
    }

    // all reading is done at this point!
    if (lcomm->MyPID() == 0) mpProblem->write_input_parameters();

    // before we destroy the reader we want to know about unused sections
    reader.print_unknown_sections();

    /***************************************************/
    /**************** create the model *****************/
    /***************************************************/
    Teuchos::RCP<Epetra_MpiComm> mpicomm = Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(gcomm);
    if (mpicomm == Teuchos::null)
        KRATOS_ERROR << "gcomm is not Epetra_MpiComm";
    mpModel = FourCModel::Pointer(new FourCModel(mpicomm->Comm(), mpProblem->n_dim()));

    // add the discretization
    const std::vector<std::string> dis_names = mpProblem->get_dis_names();
    for (auto it = dis_names.begin(); it != dis_names.end(); ++it)
    {
        auto pdisc = mpProblem->get_dis(*it);
        mpModel->AddDiscretization(pdisc);
    }
}  // end of ntainp_ccadiscret()

void FourCProblem::setup_parallel_output(const std::string& outputfile_kenner, Teuchos::RCP<Epetra_Comm> lcomm, int group)
{
    using namespace FourC;

    // configure the parallel output environment
    const Teuchos::ParameterList& io = Global::Problem::instance()->io_params();
    bool screen = io.get<bool>("WRITE_TO_SCREEN");
    bool file = io.get<bool>("WRITE_TO_FILE");
    bool preGrpID = io.get<bool>("PREFIX_GROUP_ID");
    int oproc = io.get<int>("LIMIT_OUTP_TO_PROC");
    auto level = Teuchos::getIntegralValue<Core::IO::Verbositylevel>(io, "VERBOSITY");

    Core::IO::cout.setup(
        screen, file, preGrpID, level, std::move(lcomm), oproc, group, outputfile_kenner);
}

} // namespace Kratos
