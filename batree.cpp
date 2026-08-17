// Batree - An MFEM-based SPM, SPMe and P2D solver

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include "operators/EChemOperator.hpp"

using namespace mfem;

int
main(int argc, char * argv[])
{
  // Initialize MPI and HYPRE.
  Mpi::Init(argc, argv);
  Hypre::Init();

  // Parse command-line options.
  std::string model = "SPM";
  std::string cell = "LGM50";
  int order = 1;
  int ode_solver_type = 21;
  real_t t_final = 3600.0;
  real_t dt = 1.0;
  int output_steps = 5;

  int precision = 8;
  std::cout.precision(precision);

  OptionsParser args(argc, argv);
  args.AddOption(&model, "-m", "--model", "Electrochemical model: SPM, SPMe, or P2D.");
  args.AddOption(&cell, "-c", "--cell", "Cell model: LGM50 or Enertech.");
  args.AddOption(&order, "-o", "--order", "Order (degree) of the finite elements.");
  args.AddOption(&ode_solver_type, "-s", "--ode-solver", ODESolver::Types.c_str());
  args.AddOption(&t_final, "-tf", "--t-final", "Final time; start time is 0.");
  args.AddOption(&dt, "-dt", "--time-step", "Time step.");
  args.AddOption(&output_steps, "-os", "--output-steps", "Output every n-th timestep.");
  args.Parse();
  if (!args.Good())
  {
    args.PrintUsage(std::cout);
    return 1;
  }

  if (Mpi::Root())
    args.PrintOptions(std::cout);

  // Define the ODE solver used for time integration. Several implicit
  // singly diagonal implicit Runge-Kutta (SDIRK) methods, as well as
  // explicit Runge-Kutta methods are available in MFEM. For now, we only
  // support implicit methods and we have only tested Backward Euler.
  std::unique_ptr<mfem::ODESolver> ode_solver = ODESolver::Select(ode_solver_type);

  // Initialise grid and layout properties dependent on the electrochemical model and FE order
  init_settings(model, cell, order);

  // Build the 1d mesh for the macro problem and tag its elements according to their region.
  // Define the parallel mesh by a partitioning of the serial mesh.
  // Once the parallel mesh is defined, the serial mesh can be deleted.
  Mesh x_smesh = Mesh::MakeCartesian1D(NX);
  for (unsigned i = 0; i < NX; i++)
    x_smesh.SetAttribute(i, i < NNE ? NE : i < NNE + NSEP ? SEP : PE);
  ParMesh * x_pmesh = new ParMesh(MPI_COMM_WORLD, x_smesh);
  x_smesh.Clear(); // the serial mesh is no longer needed

  // Build the 1d mesh for the micro problems, i.e. the particles.
  // Define the parallel mesh by a partitioning of the serial mesh.
  // Once the parallel mesh is defined, the serial mesh can be deleted.
  Mesh r_smesh = Mesh::MakeCartesian1D(NR);
  ParMesh * r_pmesh = new ParMesh(MPI_COMM_WORLD, r_smesh);
  r_smesh.Clear(); // the serial mesh is no longer needed

  // Define the H1 finite element spaces representing concentrations/potentials
  H1_FECollection fe_coll(order, /*dim*/ 1);
  ParFiniteElementSpace * x_h1space = new ParFiniteElementSpace(x_pmesh, &fe_coll);
  ParFiniteElementSpace * r_h1space = new ParFiniteElementSpace(r_pmesh, &fe_coll);

  // Get the total number of dofs in the system (including boundaries), for
  // both the macro and micro problems, across all processors. This is for
  // reporting purposes only.
  HYPRE_BigInt fe_size_global = NMACRO * x_h1space->GlobalTrueVSize() +
                                NPAR * r_h1space->GlobalTrueVSize();

  // Get the number of dofs in the system (including boundaries), for
  // both the macro and micro problems, _owned_ by this processor.
  HYPRE_BigInt fe_size_owned = NMACRO * x_h1space->GetTrueVSize() +
                               NPAR * r_h1space->GetTrueVSize();

  if (Mpi::Root())
  {
    std::cout << std::endl;
    std::cout << "# vars: " << (SPM ? NPAR : SPMe ? NMACROC + NPAR : P2D ? NEQS : 0) << std::endl;
    std::cout << "# dofs (total): "
              << fe_size_global - (SPMe ? NMACROP * x_h1space->GlobalTrueVSize() : 0) << std::endl;
    std::cout << "# dofs (rank 0): "
              << fe_size_owned - (SPMe ? NMACROP * x_h1space->GetTrueVSize() : 0) << std::endl;
    std::cout << std::endl;
  }

  // Initialize the ElectroChemistry operator.
  real_t t = 0.0;
  BlockVector x;
  EChemOperator oper(x_h1space, r_h1space, fe_size_owned, x);

  // Perform time-integration (looping over the time iterations, ti, with a
  // time-step dt).
  oper.SetImplicitVariableType(TimeDependentOperator::STATE);
  ode_solver->Init(oper);

  bool last_step = false;
  for (int ti = 1; !last_step; ti++)
  {
    last_step = t + dt >= t_final - dt / 2;

    ode_solver->Step(x, t, dt);

    if (output_steps && ti == 1 && Mpi::Root())
      std::cout << "step\ttime[s]\tvoltage[V]\tSoC[%]" << std::endl;

    if (output_steps && (last_step || (ti % output_steps) == 0))
    {
      real_t V = oper.GetVoltage();
      real_t SoC = oper.GetSoC();
      if (Mpi::Root())
      {
        std::cout << std::left << ti << "\t" << t << "\t";
        std::cout << std::setw(8) << V << "\t" << SoC << std::endl;
      }
    }
  }

  // Free the used memory.
  delete x_pmesh;
  delete r_pmesh;

  return 0;
}
