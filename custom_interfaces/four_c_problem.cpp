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

#include "4C_ale_dyn.hpp"
#include "4C_art_net_dyn_drt.hpp"
#include "4C_ehl_dyn.hpp"
#include "4C_elch_dyn.hpp"
#include "4C_elemag_dyn.hpp"
#include "4C_fluid_dyn_nln_drt.hpp"
#include "4C_fpsi_dyn.hpp"
#include "4C_fs3i_dyn.hpp"
#include "4C_fsi_dyn.hpp"
#include "4C_global_data.hpp"
#include "4C_immersed_problem_dyn.hpp"
#include "4C_levelset_dyn.hpp"
#include "4C_loma_dyn.hpp"
#include "4C_lubrication_dyn.hpp"
#include "4C_particle_algorithm_sim.hpp"
#include "4C_pasi_dyn.hpp"
#include "4C_poroelast_dyn.hpp"
#include "4C_poroelast_scatra_dyn.hpp"
#include "4C_porofluidmultiphase_dyn.hpp"
#include "4C_poromultiphase_dyn.hpp"
#include "4C_poromultiphase_scatra_dyn.hpp"
#include "4C_red_airways_dyn_drt.hpp"
#include "4C_scatra_cardiac_monodomain_dyn.hpp"
#include "4C_scatra_dyn.hpp"
#include "4C_ssi_dyn.hpp"
#include "4C_ssti_dyn.hpp"
#include "4C_sti_dyn.hpp"
#include "4C_stru_multi_microstatic_npsupport.hpp"
#include "4C_structure_dyn_nln_drt.hpp"
#include "4C_thermo_dyn.hpp"
#include "4C_tsi_dyn.hpp"

#include <Kokkos_Core.hpp>

// Project includes
#include "four_c_problem.h"


namespace Kratos
{

FourCProblem::FourCProblem(int argc, char** argv)
{
    // MPI_Init(&argc, &argv); // this is called outside by Kratos
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
    // MPI_Finalize(); // this is called outside by Kratos
}

std::vector<std::string> FourCProblem::GetDiscretizationNames() const
{
    return mpProblem->get_dis_names();
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
    /**************** 4C setup part *****************/
    /** This part is copied from                   **/
    /**   4C/apps/global_full/4C_global_full_main  **/
    /** Please keep track of possible changes      **/
    /************************************************/
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

    /*******************************************************************/
    /**************** interface part: create the model *****************/
    /*******************************************************************/
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
    const Teuchos::ParameterList& io = mpProblem->io_params();
    bool screen = io.get<bool>("WRITE_TO_SCREEN");
    bool file = io.get<bool>("WRITE_TO_FILE");
    bool preGrpID = io.get<bool>("PREFIX_GROUP_ID");
    int oproc = io.get<int>("LIMIT_OUTP_TO_PROC");
    auto level = Teuchos::getIntegralValue<Core::IO::Verbositylevel>(io, "VERBOSITY");

    Core::IO::cout.setup(
        screen, file, preGrpID, level, std::move(lcomm), oproc, group, outputfile_kenner);
}

void FourCProblem::Run()
{
    /********************* 4C setup part ********************/
    /** This part is copied from                           **/
    /**   apps/global_full/4C_global_full_cal_control.cpp  **/
    /** Please keep track of possible updates              **/
    /********************************************************/
  using namespace FourC;

  int restart = mpProblem->restart();

  // choose the entry-routine depending on the problem type
  switch (mpProblem->get_problem_type())
  {
    case Core::ProblemType::structure:
    case Core::ProblemType::polymernetwork:
      caldyn_drt();
      break;
    case Core::ProblemType::fluid:
    case Core::ProblemType::fluid_redmodels:
      dyn_fluid_drt(restart);
      break;
    case Core::ProblemType::lubrication:
      lubrication_dyn(restart);
      break;
    case Core::ProblemType::ehl:
      ehl_dyn();
      break;
    case Core::ProblemType::scatra:
      scatra_dyn(restart);
      break;
    case Core::ProblemType::cardiac_monodomain:
      scatra_cardiac_monodomain_dyn(restart);
      break;
    case Core::ProblemType::sti:
      sti_dyn(restart);
      break;
    case Core::ProblemType::fluid_xfem:
      fluid_xfem_drt();
      break;
      break;
    case Core::ProblemType::fluid_ale:
      fluid_ale_drt();
      break;

    case Core::ProblemType::fsi:
    case Core::ProblemType::fsi_redmodels:
      fsi_ale_drt();
      break;
    case Core::ProblemType::fsi_xfem:
      xfsi_drt();
      break;
    case Core::ProblemType::fpsi_xfem:
      xfpsi_drt();
      break;
    case Core::ProblemType::gas_fsi:
    case Core::ProblemType::biofilm_fsi:
    case Core::ProblemType::thermo_fsi:
    case Core::ProblemType::fps3i:
      fs3i_dyn();
      break;
    case Core::ProblemType::fbi:
      fsi_immersed_drt();
      break;

    case Core::ProblemType::ale:
      dyn_ale_drt();
      break;

    case Core::ProblemType::thermo:
      thr_dyn_drt();
      break;

    case Core::ProblemType::tsi:
      tsi_dyn_drt();
      break;

    case Core::ProblemType::loma:
      loma_dyn(restart);
      break;

    case Core::ProblemType::elch:
      elch_dyn(restart);
      break;

    case Core::ProblemType::art_net:
      dyn_art_net_drt();
      break;

    case Core::ProblemType::red_airways:
      dyn_red_airways_drt();
      break;

    case Core::ProblemType::immersed_fsi:
      immersed_problem_drt();
      break;

    case Core::ProblemType::poroelast:
      poroelast_drt();
      break;
    case Core::ProblemType::poroscatra:
      poro_scatra_drt();
      break;
    case Core::ProblemType::porofluidmultiphase:
      porofluidmultiphase_dyn(restart);
      break;
    case Core::ProblemType::poromultiphase:
      poromultiphase_dyn(restart);
      break;
    case Core::ProblemType::poromultiphasescatra:
      poromultiphasescatra_dyn(restart);
      break;
    case Core::ProblemType::fpsi:
      fpsi_drt();
      break;
    case Core::ProblemType::ssi:
      ssi_drt();
      break;
    case Core::ProblemType::ssti:
      ssti_drt();
      break;
    case Core::ProblemType::redairways_tissue:
      redairway_tissue_dyn();
      break;

    case Core::ProblemType::particle:
      particle_drt();
      break;

    case Core::ProblemType::pasi:
      pasi_dyn();
      break;

    case Core::ProblemType::level_set:
      levelset_dyn(restart);
      break;

    case Core::ProblemType::np_support:
      MultiScale::np_support_drt();
      break;

    case Core::ProblemType::elemag:
      electromagnetics_drt();
      break;

    default:
      FOUR_C_THROW("solution of unknown problemtyp %d requested",
          mpProblem->get_problem_type());
      break;
  }
}

} // namespace Kratos
