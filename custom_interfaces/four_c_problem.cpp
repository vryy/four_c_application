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

#include <Kokkos_Core.hpp>

// Project includes
#include "four_c_problem.h"


namespace Kratos
{

FourCProblem::FourCProblem(int argc, char** argv)
{
    /************************* 4C setup part *********************************/
    /** This constructor is adapted from                                    **/
    /**   4C/apps/global_full/4C_global_full_main.cpp                       **/
    /** The maintainer is responsible for keeping track of possible changes **/
    /*************************************************************************/
    // MPI_Init(&argc, &argv); // this is called outside by Kratos
    Kokkos::ScopeGuard kokkos_guard(argc, argv); // DO NOT KNOW WHAT IS THIS

    std::shared_ptr<FourC::Core::Communication::Communicators> communicators =
        FourC::Core::Communication::create_comm(std::vector<std::string>(argv, argv + argc));

    mArguments.input_file_name = "";
    mArguments.output_file_identifier = "";
    mArguments.restart_file_identifier = "";
    mArguments.restart_step = 0;
    mArguments.comms = communicators;

    ///
    mpProblem = FourC::Global::Problem::instance();
    ///

    if (strcmp(argv[argc - 1], "--interactive") == 0)
    {
      char hostname[256];
      gethostname(hostname, sizeof(hostname));
      printf("Global rank %d with PID %d on %s is ready for attach\n",
          FourC::Core::Communication::my_mpi_rank(mArguments.comms->global_comm()), getpid(), hostname);
      if (FourC::Core::Communication::my_mpi_rank(mArguments.comms->global_comm()) == 0)
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

    FourC::Core::Communication::barrier(mArguments.comms->global_comm());

    if ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0)))
    {
      if (FourC::Core::Communication::my_mpi_rank(mArguments.comms->local_comm()) == 0)
      {
        printf("\n\n");
        print_help_message();
        printf("\n\n");
      }
    }
    else if ((argc == 2) && ((strcmp(argv[1], "-p") == 0) || (strcmp(argv[1], "--parameters") == 0)))
    {
      if (FourC::Core::Communication::my_mpi_rank(mArguments.comms->local_comm()) == 0)
      {
        ryml::Tree tree = FourC::Core::IO::init_yaml_tree_with_exceptions();
        ryml::NodeRef root = tree.rootref();
        root |= ryml::MAP;
        FourC::Core::IO::YamlNodeRef root_ref(root, "");

        // Write the non-user input metadata that is defined globally for 4C.
        FourC::Global::emit_general_metadata(root_ref);

        // Write the user input defined for various physics module.
        FourC::Core::IO::InputFile input_file = FourC::Global::set_up_input_file(mArguments.comms->local_comm());
        input_file.emit_metadata(root_ref);

        // Finally, dump everything.
        std::cout << tree;
      }
    }
    else
    {
      if (FourC::Core::Communication::my_mpi_rank(mArguments.comms->global_comm()) == 0)
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
        print_centered(FourC::VersionControl::git_hash);
        print_centered("");
        std::cout << std::string(box_width, '*') << '\n';
        std::cout << '\n';

        std::cout << "Trilinos Version: " << FOUR_C_TRILINOS_HASH << " (git SHA1)\n";
        std::cout << "Total number of MPI ranks: "
                  << FourC::Core::Communication::num_mpi_ranks(mArguments.comms->global_comm()) << '\n';
      }

      parse_commandline_arguments(argc, argv, mArguments);
      // KRATOS_WATCH(mArguments.input_file_name)
      // KRATOS_WATCH(mArguments.output_file_identifier)
      // KRATOS_WATCH(mArguments.restart_file_identifier)
      // KRATOS_WATCH(mArguments.restart_step)

      /* input phase, input of all information */
      FourC::global_legacy_module_callbacks().RegisterParObjectTypes();

      // and now the actual reading
      FourC::Core::IO::InputFile input_file = FourC::Global::set_up_input_file(mArguments.comms->local_comm());
      input_file.read(mArguments.input_file_name);
      setup_global_problem(input_file, mArguments);

      // we wait till all procs are here. Otherwise a hang up might occur where
      // one proc ended with FOUR_C_THROW but other procs were not finished and waited...
      // we also want to have the printing above being finished.
      FourC::Core::Communication::barrier(mArguments.comms->local_comm());
    }
}

FourCProblem::~FourCProblem()
{
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
        run(mArguments);
    }
    catch (FourC::Core::Exception& err)
    {
        char line[] = "=========================================================================\n";
        std::cout << "\n\n"
                  << line << err.what_with_stacktrace() << "\n"
                  << line << "\n"
                  << std::endl;

        if (mArguments.comms->num_groups() > 1)
        {
          printf("Global processor %d has thrown an error and is waiting for the remaining procs\n\n",
              FourC::Core::Communication::my_mpi_rank(mArguments.comms->global_comm()));
          FourC::Core::Communication::barrier(mArguments.comms->global_comm());
        }

        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void FourCProblem::setup_global_problem(FourC::Core::IO::InputFile& input_file, const CommandlineArguments& arguments)
{
    /************************* 4C setup part *********************************/
    /** This function is adapted from                                       **/
    /**   4C/apps/global_full/4C_global_full_io.cpp                         **/
    /** The maintainer is responsible for keeping track of possible changes **/
    /*************************************************************************/
    using namespace FourC;
    mpProblem->set_restart_step(arguments.restart_step);
    mpProblem->set_communicators(arguments.comms);
    Global::read_parameter(*mpProblem, input_file);

    setup_parallel_output(arguments);

    // create control file for output and read restart data if required
    mpProblem->open_control_file(arguments.comms->local_comm(), arguments.input_file_name,
        arguments.output_file_identifier, arguments.restart_file_identifier);

    // input of materials
    Global::read_materials(*mpProblem, input_file);

    // input of materials
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

    Global::read_fields(*mpProblem, input_file);

    /*******************************************************************/
    /**************** interface part: create the model *****************/
    /*******************************************************************/
    MPI_Comm gcomm = mpProblem->get_communicators()->global_comm();
    mpModel = FourCModel::Pointer(new FourCModel(gcomm, mpProblem->n_dim()));

    // add the discretization
    const std::vector<std::string> dis_names = mpProblem->get_dis_names();
    for (auto it = dis_names.begin(); it != dis_names.end(); ++it)
    {
        auto pdisc = mpProblem->get_dis(*it);
        mpModel->AddDiscretization(pdisc);
    }
}

void FourCProblem::setup_parallel_output(const CommandlineArguments& arguments)
{
    using namespace FourC;

    // configure the parallel output environment
    const Teuchos::ParameterList& io = mpProblem->io_params();
    bool screen = io.get<bool>("WRITE_TO_SCREEN");
    bool file = io.get<bool>("WRITE_TO_FILE");
    bool preGrpID = io.get<bool>("PREFIX_GROUP_ID");
    int oproc = io.get<int>("LIMIT_OUTP_TO_PROC");
    auto level = Teuchos::getIntegralValue<Core::IO::Verbositylevel>(io, "VERBOSITY");

    Core::IO::cout.setup(screen, file, preGrpID, level, std::move(arguments.comms->local_comm()),
        oproc, arguments.comms->group_id(), arguments.output_file_identifier);
}

void FourCProblem::run(const CommandlineArguments& arguments)
{
    /************************* 4C setup part *********************************/
    /** This function is adapted from                                       **/
    /**   4C/apps/global_full/4C_global_full_main.cpp                       **/
    /** The maintainer is responsible for keeping track of possible changes **/
    /*************************************************************************/
    using namespace FourC;

    double t0 = walltime_in_seconds();

    const double ti = walltime_in_seconds() - t0;
    if (Core::Communication::my_mpi_rank(arguments.comms->global_comm()) == 0)
    {
      Core::IO::cout << "\nTotal wall time for INPUT:       " << std::setw(10) << std::setprecision(3)
                     << std::scientific << ti << " sec \n\n";
    }

    /*--------------------------------------------------calculation phase */
    t0 = walltime_in_seconds();

    entrypoint_switch();

    const double tc = walltime_in_seconds() - t0;
    if (Core::Communication::my_mpi_rank(arguments.comms->global_comm()) == 0)
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

std::vector<std::string> FourCProblem::parse_input_output_files(const int argc, char** argv, const int my_rank) const
{
    if (argc < 1)
    {
      if (my_rank == 0)
      {
        printf("You forgot to give the input and output file names!\n");
        printf("Try again!\n");
      }
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    else if (argc < 2)
    {
      if (my_rank == 0)
      {
        printf("You forgot to give the output file name!\n");
        printf("Try again!\n");
      }
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    // parse command line and separate input/output arguments
    std::vector<std::string> inout;
    for (int i = 0; i < argc; i++)
    {
      std::string temp = argv[i];
      if (temp.substr(0, 1) != "-") inout.push_back(temp);
    }
    return inout;
}

void FourCProblem::parse_commandline_arguments(const int argc, char** argv, CommandlineArguments& arguments) const
{
    using namespace FourC;

    int group = arguments.comms->group_id();

    int restart_group = 0;
    int my_rank = FourC::Core::Communication::my_mpi_rank(arguments.comms->local_comm());

    std::vector<std::string> inout =
        parse_input_output_files(argc, argv, my_rank);
    KRATOS_WATCH_STD_CON(inout)

    // number of input/output arguments specified by the user
    auto inout_args = int(inout.size());

    std::string input_filename;
    std::string output_file_identifier;
    std::string restart_file_identifier;
    // set input file name in each group
    switch (arguments.comms->np_type())
    {
      case Core::Communication::NestedParallelismType::no_nested_parallelism:
        input_filename = inout[0];
        output_file_identifier = inout[1];
        restart_group = 0;
        break;
      case Core::Communication::NestedParallelismType::every_group_read_input_file:
      {
        if (inout_args > 4)
          FOUR_C_THROW(
              "You specified too many arguments ({}). A maximum of four args is allowed", inout_args);

        input_filename = inout[0];
        // check whether output_file_identifier includes a dash and in case separate the number at the
        // end
        size_t pos = inout[1].rfind('-');
        if (pos != std::string::npos)
        {
          int number = atoi(inout[1].substr(pos + 1).c_str());
          inout[1] = inout[1].substr(0, pos);
          output_file_identifier = std::format("{}_group_{}_{}", inout[1], group, number);
        }
        else
        {
          output_file_identifier = std::format("{}_group_{}", inout[1], group);
        }
        restart_group = 0;
      }
      break;
      case Core::Communication::NestedParallelismType::separate_input_files:
        if (inout_args % arguments.comms->num_groups() != 0)
          FOUR_C_THROW("Each group needs the same number of arguments for input/output.");
        inout_args /= arguments.comms->num_groups();
        input_filename = inout[group * inout_args];
        output_file_identifier = inout[group * inout_args + 1];
        restart_group = group;
        break;
      default:
        FOUR_C_THROW(
            "-nptype is not correct. Only everyGroupReadInputFile and separateInputFiles "
            "are available");
        break;
    }

    if (my_rank == 0)
    {
      std::cout << "input is read from     " << input_filename << std::endl;
    }
    parse_restart_definition(
        inout, inout_args, restart_file_identifier, output_file_identifier, restart_group, arguments);

    /// set IO file names and identifiers
    arguments.input_file_name = input_filename;
    arguments.output_file_identifier = output_file_identifier;
    arguments.restart_file_identifier = restart_file_identifier;
}

void FourCProblem::parse_restart_definition(const std::vector<std::string>& inout, const int in_out_args,
    std::string& restart_file_identifier, const std::string& outfile_identifier,
    const int restart_group, CommandlineArguments& arguments) const
{
    using namespace FourC;

    //  bool parameter defining if input argument is given
    bool restartIsGiven = false;
    bool restartfromIsGiven = false;

    // default case is an identical restartfile_identifier and outputfile_identifier
    restart_file_identifier = outfile_identifier;
    for (int i = 2; i < in_out_args; i++)
    {
      std::string restart = inout[restart_group * in_out_args + i];

      if (restart.substr(0, 8) == "restart=")
      {
        const std::string option = restart.substr(8, std::string::npos);
        int r;
        if (option.compare("last_possible") == 0)
        {
          r = -1;  // here we use a negative value to trigger the search in the control file in
                   // the later step. It does not mean a restart from a negative number is allowed
                   // from the user point of view.
        }
        else
        {
          r = atoi(option.c_str());
          if (r < 0) FOUR_C_THROW("Restart number must be a positive value");
        }
        // tell the global problem about the restart step given in the command line
        arguments.restart_step = r;
        restartIsGiven = true;
      }
      else if (restart.substr(0, 12) == "restartfrom=")
      {
        restart_file_identifier = (restart.substr(12, std::string::npos).c_str());

        switch (arguments.comms->np_type())
        {
          case Core::Communication::NestedParallelismType::no_nested_parallelism:
          case Core::Communication::NestedParallelismType::separate_input_files:
            // nothing to add to restartfileidentifier
            break;
          case Core::Communication::NestedParallelismType::every_group_read_input_file:
          {
            // check whether restartfileidentifier includes a dash and in case separate the number
            // at the end
            size_t pos = restart_file_identifier.rfind('-');
            if (pos != std::string::npos)
            {
              int number = atoi(restart_file_identifier.substr(pos + 1).c_str());
              std::string identifier = restart_file_identifier.substr(0, pos);
              restart_file_identifier =
                  std::format("{}_group_{}_-{}", identifier, arguments.comms->group_id(), number);
            }
            else
            {
              restart_file_identifier =
                  std::format("{}_group_{}", restart_file_identifier, arguments.comms->group_id());
            }
          }
          break;
          default:
            FOUR_C_THROW(
                "-nptype is not correct. Only everyGroupReadInputFile and "
                "separateInputFiles are available");
            break;
        }

        restartfromIsGiven = true;
      }
    }

    // throw error in case restartfrom is given but no restart step is specified
    if (restartfromIsGiven && !restartIsGiven)
    {
      FOUR_C_THROW("You need to specify a restart step when using restartfrom.");
    }
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
