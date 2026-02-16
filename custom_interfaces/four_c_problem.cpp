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
#include "4C_config_revision.hpp"
#include "4C_comm_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_global_data_read.hpp"
#include "4C_global_legacy_module.hpp"
#include "4C_io_command_line_helpers.hpp"
#include "4C_io_input_file.hpp"
#include "4C_io_pstream.hpp"
#include "4C_utils_parameter_list.hpp"
#include "4C_utils_singleton_owner.hpp"

#include "4C_ale_dyn.hpp"
#include "4C_art_net_dyn_drt.hpp"
#include "4C_ehl_dyn.hpp"
#include "4C_elch_dyn.hpp"
#include "4C_fluid_dyn_nln_drt.hpp"
#include "4C_fpsi_dyn.hpp"
#include "4C_fs3i_dyn.hpp"
#include "4C_fsi_dyn.hpp"
#include "4C_global_data.hpp"
#include "4C_levelset_dyn.hpp"
#include "4C_loma_dyn.hpp"
#include "4C_lubrication_dyn.hpp"
#include "4C_particle_algorithm_sim.hpp"
#include "4C_pasi_dyn.hpp"
#include "4C_poroelast_dyn.hpp"
#include "4C_poroelast_scatra_dyn.hpp"
#include "4C_porofluid_pressure_based_dyn.hpp"
#include "4C_porofluid_pressure_based_elast_dyn.hpp"
#include "4C_porofluid_pressure_based_elast_scatra_dyn.hpp"
#include "4C_red_airways_dyn_drt.hpp"
#include "4C_reduced_lung_main.hpp"
#include "4C_scatra_cardiac_monodomain_dyn.hpp"
#include "4C_scatra_dyn.hpp"
#include "4C_ssi_dyn.hpp"
#include "4C_ssti_dyn.hpp"
#include "4C_sti_dyn.hpp"
#include "4C_stru_multi_microstatic_npsupport.hpp"
#include "4C_structure_dyn_nln_drt.hpp"
#include "4C_thermo_dyn.hpp"
#include "4C_tsi_dyn.hpp"
#include "4C_utils_enum.hpp"

// #include <Kokkos_Core.hpp>
#include <CLI/CLI.hpp>

// Project includes
#include "four_c_problem.h"


namespace Kratos
{

// Custom CLI11 formatter to add extra spacing between options
class SpacedFormatter : public CLI::Formatter
{
 public:
  SpacedFormatter() : Formatter() {}

  std::string make_option(const CLI::Option* opt, bool in_sub) const override
  {
    std::string s = Formatter::make_option(opt, in_sub);
    if (!s.empty() && s.back() == '\n')
      s += '\n';
    else
      s += "\n\n";
    return s;
  }
};

/**
 * \brief Build canonical input/output pairs from positional command line arguments.
 * Due to the legacy argument structure, we separate between the primary input and output (first two
 * positional arguments) and the rest (io_pairs). The latter are only required when using nested
 * parallelism with separate input files.
 * \param io_pairs Vector of strings from the command line representing input/output pairs.
 * \param primary_input The primary input file name (first positional argument).
 * \param primary_output The primary output file identifier (second positional argument).
 * \return A vector of pairs of input file paths and output file identifiers.
 */
std::vector<std::pair<std::filesystem::path, std::string>> build_io_pairs(
    const std::vector<std::string>& io_pairs, const std::filesystem::path& primary_input,
    const std::string& primary_output)
{
  using namespace FourC;

  std::vector<std::pair<std::filesystem::path, std::string>> io_pairs_new;

  io_pairs_new.emplace_back(primary_input, primary_output);

  if (!io_pairs.empty())
  {
    if (io_pairs.size() % 2 != 0)
    {
      FOUR_C_THROW("Positional arguments must be provided as pairs: <input> <output>.\n");
    }
    for (size_t i = 0; i < io_pairs.size(); i += 2)
      io_pairs_new.emplace_back(std::filesystem::path(io_pairs[i]), io_pairs[i + 1]);
  }
  return io_pairs_new;
}

/**
 * \brief Validates cross-compatibility of command line options.
 * \param arguments The parsed command line arguments.
 */
void validate_argument_cross_compatibility(const FourCProblem::CommandlineArguments& arguments)
{
  using namespace FourC;
  using NPT = Core::Communication::NestedParallelismType;

  if (!arguments.group_layout.empty())
  {
    const int layout_len = static_cast<int>(arguments.group_layout.size());
    if (arguments.n_groups != layout_len)
    {
      FOUR_C_THROW(
          "When --glayout is provided its number of entries must equal --ngroup.\n "
          "Example mpirun -np 4 ./4C --ngroup=2 --glayout=1,3 \n");
    }
  }

  if (arguments.n_groups > 1 && arguments.nptype == NPT::no_nested_parallelism)
  {
    FOUR_C_THROW("when --ngroup > 1, a nested parallelism type must be specified via --nptype.\n");
  }

  if (!arguments.parameters)
  {
    const size_t num_pairs = arguments.io_pairs.size();
    if (arguments.nptype == NPT::no_nested_parallelism ||
        arguments.nptype == NPT::every_group_read_input_file)
    {
      if (num_pairs != 1)
      {
        FOUR_C_THROW(
            "when using 'no_nested_parallelism' or 'everyGroupReadInputFile' the "
            "number of <input> <output> pairs must be exactly 1.\n");
      }
    }
    else if (arguments.nptype == NPT::separate_input_files ||
             arguments.nptype == NPT::nested_multiscale)
    {
      if (static_cast<int>(num_pairs) != arguments.n_groups)
      {
        FOUR_C_THROW(
            "when using 'separateInputFiles' or 'nestedMultiscale' the number of "
            "<input> <output> pairs must equal --ngroup {}.\n",
            arguments.n_groups);
      }
    }
  }

  if (arguments.nptype != NPT::separate_input_files &&
      (arguments.restart_per_group.size() > 1 || arguments.restart_identifier_per_group.size() > 1))
  {
    FOUR_C_THROW(
        "When using --nptype other than 'separateInputFiles', only one restart step and one "
        "restartfrom identifier must be given.");
  }

  for (size_t i = 0; i < arguments.restart_identifier_per_group.size(); ++i)
  {
    if (i >= arguments.restart_per_group.size())
    {
      FOUR_C_THROW("You need to specify a restart step when using restartfrom.");
    }
  }
}

/**
 * \brief Updates input/output identifiers based on group id and nested parallelism type.
 * \param arguments The command line arguments to update.
 * \param group The group id of the current process.
 */
void update_io_identifiers(FourCProblem::CommandlineArguments& arguments, int group)
{
  using namespace FourC;
  using NPT = Core::Communication::NestedParallelismType;

  std::filesystem::path input_filename;
  std::string output_file_identifier;

  int restart_input_index = (arguments.nptype == NPT::separate_input_files) ? group : 0;

  arguments.restart =
      arguments.restart_per_group.empty() ? 0 : arguments.restart_per_group[restart_input_index];
  std::string restart_file_identifier =
      arguments.restart_identifier_per_group.empty()
          ? ""
          : arguments.restart_identifier_per_group[restart_input_index];

  switch (arguments.nptype)
  {
    case NPT::no_nested_parallelism:
      input_filename = arguments.io_pairs[0].first;
      output_file_identifier = arguments.io_pairs[0].second;
      if (restart_file_identifier == "")
      {
        restart_file_identifier = output_file_identifier;
      }
      break;
    case NPT::every_group_read_input_file:
    {
      input_filename = arguments.io_pairs[0].first;
      std::string output_file_identifier_temp = arguments.io_pairs[0].second;
      // check whether output_file_identifier includes a dash and in case separate the number at the
      // end
      size_t pos = output_file_identifier_temp.rfind('-');
      auto extract_number_and_identifier = [](const std::string& str, size_t pos)
      {
        std::string number_str = str.substr(pos + 1);
        std::string identifier = str.substr(0, pos);
        int number = 0;
        try
        {
          size_t idx = 0;
          number = std::stoi(number_str, &idx);
          if (idx != number_str.size())
          {
            FOUR_C_THROW("Invalid numeric value in output identifier: '{}'", number_str);
          }
        }
        catch (const std::exception& e)
        {
          FOUR_C_THROW(
              "Failed to parse number in output identifier '{}': {}", number_str, e.what());
        }
        return std::make_pair(identifier, number);
      };
      if (pos != std::string::npos)
      {
        auto [identifier, number] = extract_number_and_identifier(output_file_identifier_temp, pos);
        output_file_identifier = std::format("{}_group_{}_{}", identifier, group, number);
      }
      else
      {
        output_file_identifier = std::format("{}_group_{}", output_file_identifier_temp, group);
      }
      size_t pos_r = restart_file_identifier.rfind('-');
      if (restart_file_identifier == "")
      {
        restart_file_identifier = output_file_identifier;
      }
      else if (pos_r != std::string::npos)
      {
        auto [identifier, number] = extract_number_and_identifier(restart_file_identifier, pos_r);
        restart_file_identifier = std::format("{}_group_{}-{}", identifier, group, number);
      }
      else
      {
        restart_file_identifier = std::format("{}_group_{}", restart_file_identifier, group);
      }
      break;
    }
    case NPT::separate_input_files:
    case NPT::nested_multiscale:
      input_filename = arguments.io_pairs[group].first;
      output_file_identifier = arguments.io_pairs[group].second;
      if (restart_file_identifier == "")
      {
        restart_file_identifier = output_file_identifier;
      }
      break;
    default:
      FOUR_C_THROW("-nptype value {} is not valid.", static_cast<int>(arguments.nptype));
      break;
  }
  arguments.input_file_name = input_filename;
  arguments.output_file_identifier = output_file_identifier;
  arguments.restart_file_identifier = restart_file_identifier;
}

/**
 * \brief Ensures a valid group layout exists.
 *
 * If group_layout is empty and n_groups > 1, computes an equal-size layout based
 * on the size of MPI_COMM_WORLD. Aborts if the processor count is not divisible by the number of
 * groups.
 * \param n_groups Requested number of groups.
 * \param group_layout In/out container for processors per group; populated when initially empty.
 */
void assign_group_layout(const int& n_groups, std::vector<int>& group_layout)
{
  int myrank = -1;
  int num_procs = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (n_groups > 1)
  {
    if (group_layout.size() == 0)
    {
      if (myrank == (num_procs - 1))  // myrank == 0 is eventually not within 4C (i.e. coupling to
                                      // external codes)
      {
        printf(
            "\n\n\nINFO: Group layout is not specified. Default is equal size of the "
            "groups.\n");
      }
      if ((num_procs % n_groups) != 0)
      {
        if (myrank == (num_procs - 1))
        {
          printf("\n\nNumber of processors (%d) cannot be divided by the number of groups (%d)!\n",
              num_procs, n_groups);
          printf("Try again!\n");
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }

      // equal size of the groups
      for (int k = 0; k < n_groups; k++)
      {
        group_layout.push_back(num_procs / n_groups);
      }
    }
  }
}

/**
 * \brief Writes the Teuchos::TimeMonitor information to std::cout
 */
void write_timemonitor(MPI_Comm comm)
{
  using namespace FourC;
  std::shared_ptr<const Teuchos::Comm<int>> TeuchosComm =
      Core::Communication::to_teuchos_comm<int>(comm);
  Teuchos::TimeMonitor::summarize(Teuchos::Ptr(TeuchosComm.get()), std::cout, false, true, false);
}

///////////////////////////////////////////////////////////////

FourCProblem::FourCProblem(int argc, char** argv)
{
    /************************* 4C setup part *********************************/
    /** This constructor is adapted from                                    **/
    /**   4C/apps/global_full/4C_global_full_main.cpp                       **/
    /** The maintainer is responsible for keeping track of possible changes **/
    /*************************************************************************/
    // MPI_Init(&argc, &argv); // this is called outside by Kratos
    // Kokkos::ScopeGuard kokkos_guard{}; // DO NOT KNOW WHAT IS THIS

    using namespace FourC;

    // // Initialize our own singleton registry to ensure we clean up all singletons properly.
    // Core::Utils::SingletonOwnerRegistry::ScopeGuard singleton_owner_guard{};

    mArguments = parse_command_line(argc, argv);

    Core::Communication::CommConfig config{
        .group_layout = mArguments.group_layout,
        .np_type = mArguments.nptype,
        .diffgroup = mArguments.diffgroup,
    };
    mpCommunicators = std::make_shared<Core::Communication::Communicators>(Core::Communication::create_comm(config));

    ///
    mpProblem = Global::Problem::instance();
    ///

    if (mArguments.interactive)
    {
      char hostname[256];
      gethostname(hostname, sizeof(hostname));
      printf("Global rank %d with PID %d on %s is ready for attach\n",
          Core::Communication::my_mpi_rank(mpCommunicators->global_comm()), getpid(), hostname);
      if (Core::Communication::my_mpi_rank(mpCommunicators->global_comm()) == 0)
      {
        printf("\n** Enter a character to continue > \n");
        fflush(stdout);
        char go = ' ';
        if (scanf("%c", &go) == EOF)
        {
          FOUR_C_THROW("Error while reading input.\n");
        }
      }
    }

    Core::Communication::barrier(mpCommunicators->global_comm());

    if (mArguments.parameters)
    {
      if (Core::Communication::my_mpi_rank(mpCommunicators->local_comm()) == 0)
      {
        ryml::Tree tree = Core::IO::init_yaml_tree_with_exceptions();
        ryml::NodeRef root = tree.rootref();
        root |= ryml::MAP;
        Core::IO::YamlNodeRef root_ref(root, "");

        // Write the non-user input metadata that is defined globally for 4C.
        Global::emit_general_metadata(root_ref);

        // Write the user input defined for various physics module.
        Core::IO::InputFile input_file = Global::set_up_input_file(mpCommunicators->local_comm());
        input_file.emit_metadata(root_ref);

        // Finally, dump everything.
        std::cout << tree;
      }
    }
    else
    {
      if (Core::Communication::my_mpi_rank(mpCommunicators->global_comm()) == 0)
      {
        constexpr int box_width = 54;

        const auto print_centered = [&](const std::string& str)
        {
          // Subtract 2 for the asterisks on either side
          constexpr int width = box_width - 2;
          FOUR_C_ASSERT(str.size() < width, "String is too long to be centered.");
          std::cout << '*' << std::format("{:^{}}", str, width) << "*\n";
        };

        std::cout << '\n';
        std::cout << std::string(box_width, '*') << '\n';
        print_centered("");
        print_centered("4C");
        print_centered("");
        print_centered("version " FOUR_C_VERSION_FULL);
        print_centered("");
        print_centered("git SHA1");
        print_centered(VersionControl::git_hash);
        print_centered("");
        std::cout << std::string(box_width, '*') << '\n';
        std::cout << '\n';

        std::cout << "Trilinos Version: " << FOUR_C_TRILINOS_HASH << " (git SHA1)\n";
        std::cout << "Total number of MPI ranks: "
                  << Core::Communication::num_mpi_ranks(mpCommunicators->global_comm()) << '\n';
      }
    }
}

FourCProblem::~FourCProblem()
{
  if (mpCommunicators != nullptr)
    mpCommunicators->finalize();
    // FourC::Core::Utils::SingletonOwnerRegistry::finalize();
    // FourC::Global::Problem::done();
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

void FourCProblem::Run()
{
    try
    {
        run(mArguments, *mpCommunicators);
    }
    catch (FourC::Core::Exception& err)
    {
        char line[] = "=========================================================================\n";
        std::cout << "\n\n"
                  << line << err.what_with_stacktrace() << "\n"
                  << line << "\n"
                  << std::endl;

        if (mpCommunicators->num_groups() > 1)
        {
          printf("Global processor %d has thrown an error and is waiting for the remaining procs\n\n",
              FourC::Core::Communication::my_mpi_rank(mpCommunicators->global_comm()));
          FourC::Core::Communication::barrier(mpCommunicators->global_comm());
        }

        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void FourCProblem::setup_global_problem(FourC::Core::IO::InputFile& input_file, const CommandlineArguments& arguments,
    const FourC::Core::Communication::Communicators& communicators)
{
    /************************* 4C setup part *********************************/
    /** This function is adapted from                                       **/
    /**   4C/apps/global_full/4C_global_full_io.cpp                         **/
    /** The maintainer is responsible for keeping track of possible changes **/
    /*************************************************************************/
    using namespace FourC;
    mpProblem->set_restart_step(arguments.restart);
    mpProblem->set_communicators(communicators);
    Global::read_parameter(*mpProblem, input_file);

    setup_parallel_output(arguments, communicators);

    // create control file for output and read restart data if required
    mpProblem->open_control_file(communicators.local_comm(), arguments.input_file_name,
        arguments.output_file_identifier, arguments.restart_file_identifier);

    // input of materials
    Global::read_materials(*mpProblem, input_file);

    // input for multi-scale rough-surface contact
    Global::read_contact_constitutive_laws(*mpProblem, input_file);

    // input of materials of cloned fields (if needed)
    Global::read_cloning_material_map(*mpProblem, input_file);

    {
        Core::Utils::FunctionManager function_manager;
        global_legacy_module_callbacks().AttachFunctionDefinitions(function_manager);
        function_manager.read_input(input_file);
        mpProblem->set_function_manager(std::move(function_manager));
    }

    // input of particles
    Global::read_particles(*mpProblem, input_file);

    // input of fields
    auto mesh_reader = Global::read_discretization(*mpProblem, input_file);
    FOUR_C_ASSERT(mesh_reader, "Internal error: nullptr.");

    // read result tests
    Global::read_result(*mpProblem, input_file);

    // read all types of geometry related conditions (e.g. boundary conditions)
    // Also read time and space functions and local coord systems
    Global::read_conditions(*mpProblem, input_file, *mesh_reader);

    // read all knot information for isogeometric analysis
    // and add it to the (derived) nurbs discretization
    Global::read_knots(*mpProblem, input_file);

    Global::read_fields(*mpProblem, input_file, *mesh_reader);

    /*******************************************************************/
    /**************** interface part: create the model *****************/
    /*******************************************************************/
    MPI_Comm gcomm = mpProblem->get_communicators().global_comm();
    mpModel = FourCModel::Pointer(new FourCModel(gcomm, mpProblem->n_dim()));

    // add the discretization
    const std::vector<std::string> dis_names = mpProblem->get_dis_names();
    for (auto it = dis_names.begin(); it != dis_names.end(); ++it)
    {
        auto pdisc = mpProblem->get_dis(*it);
        mpModel->AddDiscretization(pdisc);
    }
}

void FourCProblem::setup_parallel_output(const CommandlineArguments& arguments,
      const FourC::Core::Communication::Communicators& communicators)
{
    using namespace FourC;

    // configure the parallel output environment
    const Teuchos::ParameterList& io = mpProblem->io_params();
    bool screen = io.get<bool>("WRITE_TO_SCREEN");
    bool file = io.get<bool>("WRITE_TO_FILE");
    bool preGrpID = io.get<bool>("PREFIX_GROUP_ID");
    int oproc = io.get<int>("LIMIT_OUTP_TO_PROC");
    auto level = Teuchos::getIntegralValue<Core::IO::Verbositylevel>(io, "VERBOSITY");

    Core::IO::cout.setup(screen, file, preGrpID, level, std::move(communicators.local_comm()),
        oproc, communicators.group_id(), arguments.output_file_identifier);
}

void FourCProblem::run(CommandlineArguments& cli_args, const FourC::Core::Communication::Communicators& communicators)
{
    /************************* 4C setup part *********************************/
    /** This function is adapted from                                       **/
    /**   4C/apps/global_full/4C_global_full_main.cpp                       **/
    /** The maintainer is responsible for keeping track of possible changes **/
    /*************************************************************************/
    using namespace FourC;

    update_io_identifiers(cli_args, communicators.group_id());

  /* input phase, input of all information */
    global_legacy_module_callbacks().RegisterParObjectTypes();
    double t0 = walltime_in_seconds();

    // and now the actual reading
    Core::IO::InputFile input_file = Global::set_up_input_file(communicators.local_comm());
    input_file.read(cli_args.input_file_name);
    setup_global_problem(input_file, cli_args, communicators);

    // we wait till all procs are here. Otherwise a hang up might occur where
    // one proc ended with FOUR_C_THROW but other procs were not finished and waited...
    // we also want to have the printing above being finished.
    Core::Communication::barrier(communicators.local_comm());

    const double ti = walltime_in_seconds() - t0;
    if (Core::Communication::my_mpi_rank(communicators.global_comm()) == 0)
    {
      Core::IO::cout << "\nTotal wall time for INPUT:       " << std::setw(10) << std::setprecision(3)
                     << std::scientific << ti << " sec \n\n";
    }

    /*--------------------------------------------------calculation phase */
    t0 = walltime_in_seconds();

    entrypoint_switch();

    write_timemonitor(communicators.local_comm());

    const double tc = walltime_in_seconds() - t0;
    if (Core::Communication::my_mpi_rank(communicators.global_comm()) == 0)
    {
      Core::IO::cout << "\nTotal wall time for CALCULATION: " << std::setw(10) << std::setprecision(3)
                     << std::scientific << tc << " sec \n\n";
    }
}

void FourCProblem::entrypoint_switch()
{
    /********************* 4C analysis switch *********************/
    /** This part is copied from                                 **/
    /**   apps/global_full/4C_global_full_entrypoint_switch.cpp  **/
    /** Please keep track of possible updates                    **/
    /**************************************************************/

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
          thermo_dyn_drt();
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

        case Core::ProblemType::reduced_lung:
          ReducedLung::reduced_lung_main();
          break;

        case Core::ProblemType::poroelast:
          poroelast_drt();
          break;
        case Core::ProblemType::poroscatra:
          poro_scatra_drt();
          break;
        case Core::ProblemType::porofluid_pressure_based:
          porofluid_pressure_based_dyn(restart);
          break;
        case Core::ProblemType::porofluid_pressure_based_elast:
          porofluid_elast_dyn(restart);
          break;
        case Core::ProblemType::porofluid_pressure_based_elast_scatra:
          porofluid_pressure_based_elast_scatra_dyn(restart);
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

        default:
          FOUR_C_THROW("solution of unknown problemtyp %d requested",
              mpProblem->get_problem_type());
          break;
    }
}

FourCProblem::CommandlineArguments FourCProblem::parse_command_line(int argc, char** argv) const
{
  using namespace FourC;

  CLI::App cli_app{"4C - Multiphysics \nComprehensive Computational Community Code"};
  cli_app.formatter(std::make_shared<SpacedFormatter>());
  CommandlineArguments arguments;
  cli_app.add_flag("-p,--parameters", arguments.parameters,
      "Dumps information about the parameters for consumption by additional tools.");
  cli_app.add_flag("-i,--interactive", arguments.interactive,
      "4C waits at the beginning for keyboard input. "
      "Helpful for parallel debugging when attaching to a single job.");
  cli_app
      .add_option("--ngroup", arguments.n_groups,
          "Specify the number of groups for nested parallelism. (default: 1)")
      ->check(CLI::PositiveNumber);
  cli_app
      .add_option(
          "--glayout",
          [&arguments](const std::vector<std::string>& tokens) -> bool
          {
            if (tokens.empty()) return false;
            arguments.group_layout.clear();

            for (const auto& tok : tokens)
            {
              std::stringstream ss(tok);
              std::string part;
              while (std::getline(ss, part, ','))
              {
                if (part.empty()) continue;
                try
                {
                  size_t pos = 0;
                  int val = std::stoi(part, &pos);
                  if (pos != part.size() || val <= 0)
                  {
                    throw CLI::ValidationError("glayout", "Entries must be positive integers.");
                  }
                  arguments.group_layout.push_back(val);
                }
                catch (const std::invalid_argument&)
                {
                  throw CLI::ValidationError("glayout", "Entries must be positive integers.");
                }
                catch (const std::out_of_range&)
                {
                  throw CLI::ValidationError("glayout", "glayout entry out of range.");
                }
              }
            }
            return true;
          },
          "Specify the number of processors per group. Comma-separated list without spaces, e.g. "
          "--glayout=<a>,<b>.\n"
          "Argument --ngroup is mandatory if a glayout is provided. (default: equal distribution)")
      ->allow_extra_args(false);
  cli_app.add_option(
      "--nptype",
      [&arguments](const std::vector<std::string>& tokens) -> bool
      {
        if (tokens.empty()) return false;
        const std::string& input = tokens.front();
        if (input == "separateInputFiles")
        {
          arguments.nptype = Core::Communication::NestedParallelismType::separate_input_files;
          return true;
        }
        else if (input == "everyGroupReadInputFile")
        {
          arguments.nptype =
              Core::Communication::NestedParallelismType::every_group_read_input_file;
          return true;
        }
        else if (input == "nestedMultiscale")
        {
          arguments.nptype = Core::Communication::NestedParallelismType::nested_multiscale;
          return true;
        }
        else if (input.rfind("diffgroup", 0) == 0)
        {
          // Expect exactly an integer suffix after "diffgroup", and only allow 0 or 1
          const std::string suffix = input.substr(9);
          if (suffix.empty())
          {
            throw CLI::ValidationError(
                "nptype", "Missing suffix for 'diffgroup'; expected 'diffgroup0' or 'diffgroup1'.");
          }
          for (char c : suffix)
          {
            if (!std::isdigit(static_cast<unsigned char>(c)))
            {
              throw CLI::ValidationError(
                  "nptype", "Invalid diffgroup suffix; expected integer after 'diffgroup'.");
            }
          }
          int val = std::stoi(suffix);
          if (val != 0 && val != 1)
          {
            throw CLI::ValidationError("nptype", "Only diffgroup0 and diffgroup1 are allowed.");
          }
          arguments.nptype = Core::Communication::NestedParallelismType::no_nested_parallelism;
          arguments.diffgroup = val;
          return true;
        }
        else
        {
          throw CLI::ValidationError("nptype",
              "Only 'everyGroupReadInputFile', 'separateInputFiles', 'nestedMultiscale', and "
              "'diffgroupx' are available for nptype.");
        }
      },
      "Specify nested parallelism type: "
      "separateInputFiles|everyGroupReadInputFile|\nnestedMultiscale|diffgroup<N> \n"
      "Must be set if --ngroup > 1. \n"
      "'diffgroupx' can be used to compare vectors/matrices/results between two separate "
      "(serial/parallel) 4C runs; x must be 0 and 1 for the respective run");
  cli_app.add_option(
      "--restart",
      [&arguments](const std::vector<std::string>& tokens) -> bool
      {
        if (tokens.empty()) return false;
        for (const auto& tok : tokens)
        {
          // allow comma-separated lists in a single token (e.g. "--restart=3,4")
          std::stringstream ss(tok);
          std::string part;
          while (std::getline(ss, part, ','))
          {
            if (part.empty()) continue;
            if (part == "last_possible")
            {
              arguments.restart_per_group.push_back(-1);  // legacy sentinel for last_possible
              continue;
            }
            try
            {
              size_t pos = 0;
              int val = std::stoi(part, &pos);
              if (pos != part.size() || val < 0)
              {
                throw CLI::ValidationError(
                    "restart", "Restart step must be a non-negative integer or 'last_possible'.");
              }
              arguments.restart_per_group.push_back(val);
            }
            catch (const std::invalid_argument&)
            {
              throw CLI::ValidationError(
                  "restart", "Restart step must be a non-negative integer or 'last_possible'.");
            }
            catch (const std::out_of_range&)
            {
              throw CLI::ValidationError("restart", "Restart step out of range.");
            }
          }
        }
        return true;
      },
      "Restart the simulation from step <y>. Accepts a non-negative integer or 'last_possible'.\n"
      "If nested parallelism with separate input files is used, each group can have a different "
      "restart step defined as a comma-separated list, e.g., --restart=<a>,<b>.");
  cli_app.add_option(
      "--restartfrom",
      [&arguments](const std::vector<std::string>& tokens) -> bool
      {
        if (tokens.empty()) return false;
        for (const auto& tok : tokens)
        {
          std::stringstream ss(tok);
          std::string part;
          while (std::getline(ss, part, ','))
          {
            if (part.empty()) continue;
            arguments.restart_identifier_per_group.push_back(part);
          }
        }
        return true;
      },
      "Restart the simulation from the files prefixed with <restart_file_name>.\n If nested "
      "parallelism with separate input files is used, each group can have a different file prefix "
      "defined as a comma-separated list.");
  std::string primary_input;
  std::string primary_output;
  cli_app.add_option("input", primary_input, "Name of the input file, including the suffix");
  cli_app.add_option("output", primary_output, "Prefix of your output files.");

  std::vector<std::string> io_pairs;
  cli_app
      .add_option("io_pairs", io_pairs,
          "More pairs of simulation <input> and <output> names. Only necessary when using nested "
          "parallelism with multiple groups and separate input files.")
      ->expected(-1);


  std::vector<std::string> raw_args;
  raw_args.reserve(argc);
  for (int i = 1; i < argc; ++i) raw_args.emplace_back(argv[i]);

  LegacyCliOptions legacy_options = {.single_dash_legacy_names = {"ngroup", "glayout", "nptype"},
      .nodash_legacy_names = {"restart", "restartfrom"}};
  std::vector<std::string> sanitized_args = adapt_legacy_cli_arguments(raw_args, legacy_options);

  // Reversed order required when parsing std::vector<string> with CLI11
  std::reverse(sanitized_args.begin(), sanitized_args.end());
  try
  {
    cli_app.parse(sanitized_args);
  }
  catch (const CLI::ParseError& e)
  {
    std::exit(cli_app.exit(e));
  }

  if (!arguments.parameters)
  {
    if (primary_input.empty() || primary_output.empty())
    {
      FOUR_C_THROW("Please provide both <input> and <output> arguments.");
    }
  }

  arguments.io_pairs = build_io_pairs(io_pairs, primary_input, primary_output);
  validate_argument_cross_compatibility(arguments);
  assign_group_layout(arguments.n_groups, arguments.group_layout);
  return arguments;
}

double FourCProblem::walltime_in_seconds() const
{
  return std::chrono::duration_cast<std::chrono::milliseconds>(
             std::chrono::high_resolution_clock::now().time_since_epoch())
             .count() *
         1.0e-3;
}

void FourCProblem::print_help_message() const
{
  std::cout
      << "NAME\n"
      << "\t"
      << "4C - simulate just about anything\n"
      << "\n"
      << "SYNOPSIS\n"
      << "\t"
      << "4C [-h | --help] [-p | --parameters] [-d | --datfile] [-ngroup=<x>] \\ "
         "\n"
         "\t\t[-glayout=a,b,c,...] [-nptype=<parallelism_type>] \\ \n"
      << "\t\t<input_name> <output_name> [restart=<y>] [restartfrom=restart_file_name] \\ \n"
         "\t\t[ <input_name0> <output_name0> [restart=<y>] [restartfrom=restart_file_name] ... "
         "] \\ \n"
         "\t\t[--interactive]\n"
      << "\n"
      << "DESCRIPTION\n"
      << "\tThe am besten simulation tool in the world.\n"
      << "\n"
      << "OPTIONS\n"
      << "\t--help or -h\n"
      << "\t\tPrint this message.\n"
      << "\n"
      << "\t--parameters or -p\n"
      << "\t\tDumps information about the parameters for consumption by additional tools.\n"
      << "\n"
      << "\t-ngroup=<x>\n"
      << "\t\tSpecify the number of groups for nested parallelism. (default: 1)\n"
      << "\n"
      << "\t-glayout=<a>,<b>,<c>,...\n"
      << "\t\tSpecify the number of processors per group. \n"
         "\t\tArgument \"-ngroup\" is mandatory and must be preceding. \n"
         "\t\t(default: equal distribution)\n"
      << "\n"
      << "\t-nptype=<parallelism_type>\n"
      << "\t\tAvailable options: \"separateInputFiles\" and \"everyGroupReadInputFile\"; \n"
         "\t\tMust be set if \"-ngroup\" > 1.\n"
      << "\t\t\"diffgroupx\" can be used to compare results from separate but parallel 4C "
         "runs; \n"
         "\t\tx must be 0 and 1 for the respective run\n"
      << "\n"
      << "\t<input_name>\n"
      << "\t\tName of the input file, including the suffix\n"
      << "\n"
      << "\t<output_name>\n"
      << "\t\tPrefix of your output files.\n"
      << "\n"
      << "\trestart=<y>\n"
      << "\t\tRestart the simulation from step <y>. \n"
         "\t\tIt always refers to the previously defined <input_name> and <output_name>. \n"
         "\t\t(default: 0 or from <input_name>)\n"
         "\t\tIf y=last_possible, it will restart from the last restart step defined in the "
         "control file.\n"
      << "\n"
      << "\trestartfrom=<restart_file_name>\n"
      << "\t\tRestart the simulation from the files prefixed with <restart_file_name>. \n"
         "\t\t(default: <output_name>)\n"
      << "\n"
      << "\t--interactive\n"
      << "\t\t4C waits at the beginning for keyboard input. \n"
         "\t\tHelpful for parallel debugging when attaching to a single job. \n"
         "\t\tMust be specified at the end in the command line.\n"
      << "\n";
}

} // namespace Kratos
