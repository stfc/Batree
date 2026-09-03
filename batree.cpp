// Batree - An MFEM-based SPM, SPMe and P2D solver

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include "operators/EChemOperator.hpp"
#include "operators/CurrentCollectorOperator.hpp"

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
  real_t c_rate = 1;
  int order = 1;
  int ode_solver_type = 21;
  real_t t_final = -1.0;
  real_t dt = 1.0;
  int output_steps = 5;

  OptionsParser args(argc, argv);
  args.AddOption(&model, "-m", "--model", "Electrochemical model: SPM, SPMe, or P2D.");
  args.AddOption(&cell, "-c", "--cell", "Cell model: LGM50 or Enertech.");
  args.AddOption(&c_rate, "-cr", "--c-rate", "The C-rate to run a constant current (dis)charge");
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

  // Initialise properties dependent on the electrochemical model, cell type, current and FE order
  init_settings(model, cell, c_rate, order);

  // Build the 1d mesh for the macro problem and tag its elements according to their region.
  // Define the parallel mesh by a partitioning of the serial mesh.
  // Once the parallel mesh is defined, the serial mesh can be deleted.
  Mesh x_smesh = Mesh::MakeCartesian1D(NX);
  for (unsigned i = 0; i < NX; i++)
    x_smesh.SetAttribute(i, i < NNE ? NE : i < NNE + NSEP ? SEP : PE);
  ParMesh x_pmesh(MPI_COMM_WORLD, x_smesh);
  x_smesh.Clear(); // the serial mesh is no longer needed

  // Build the 1d mesh for the micro problems, i.e. the particles.
  // Define the parallel mesh by a partitioning of the serial mesh.
  // Once the parallel mesh is defined, the serial mesh can be deleted.
  Mesh r_smesh = Mesh::MakeCartesian1D(NR);
  ParMesh r_pmesh(MPI_COMM_WORLD, r_smesh);
  r_smesh.Clear(); // the serial mesh is no longer needed

  // Define the H1 finite element spaces representing concentrations/potentials
  H1_FECollection fe_coll(order, /*dim*/ 1);
  ParFiniteElementSpace x_h1space(&x_pmesh, &fe_coll);
  ParFiniteElementSpace r_h1space(&r_pmesh, &fe_coll);

  // Get the total number of dofs in the system (including boundaries), for
  // both the macro and micro problems, across all processors. This is for
  // reporting purposes only.
  HYPRE_BigInt fe_size_global =
      NMACRO * x_h1space.GlobalTrueVSize() + NPAR * r_h1space.GlobalTrueVSize();

  // Get the number of dofs in the system (including boundaries), for
  // both the macro and micro problems, _owned_ by this processor.
  HYPRE_BigInt fe_size_owned = NMACRO * x_h1space.GetTrueVSize() + NPAR * r_h1space.GetTrueVSize();

  if (Mpi::Root())
  {
    std::cout << std::endl;
    std::cout << "# vars: " << (SPM ? NPAR : SPMe ? NMACROC + NPAR : P2D ? NEQS : 0) << std::endl;
    std::cout << "# dofs (total): "
              << fe_size_global - (SPMe ? NMACROP * x_h1space.GlobalTrueVSize() : 0) << std::endl;
    std::cout << "# dofs (rank 0): "
              << fe_size_owned - (SPMe ? NMACROP * x_h1space.GetTrueVSize() : 0) << std::endl;
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

  //Collector positive potentials (pp)
  Mesh current_collector_positive_mesh("../mesh/rectangle_current_collector_pos.msh", 1, 0);
  ParMesh * pp_pmesh = new ParMesh(MPI_COMM_WORLD, current_collector_positive_mesh);

  // Define the H1 finite element spaces representing potentials in the positive current collector
  H1_FECollection fe_current_collector_positive_coll(order, pp_pmesh->Dimension());

  // Create FiniteElementSpace for the positive current collector
  ParFiniteElementSpace * current_collector_positive_h1space;
  current_collector_positive_h1space = new ParFiniteElementSpace(pp_pmesh, &fe_current_collector_positive_coll);

  //Collector negative potentials (np)
  Mesh current_collector_negative_mesh("../mesh/rectangle_current_collector_neg.msh", 1, 0);
  ParMesh * np_pmesh = new ParMesh(MPI_COMM_WORLD, current_collector_negative_mesh);

  // Define the H1 finite element spaces representing potentials in the negative current collector
  H1_FECollection fe_current_collector_negative_coll(order, np_pmesh->Dimension());

  // Create FiniteElementSpace for the negative current collector
  ParFiniteElementSpace * current_collector_negative_h1space;
  current_collector_negative_h1space = new ParFiniteElementSpace(np_pmesh, &fe_current_collector_negative_coll);

  //This is the domain attribute (without tab) for the positive/negative current collector
  mfem::Array<int> domain_attributes;
  domain_attributes.Append(1); 

  // domain currents mesh
  ParSubMesh dc_pmesh_currents(
    ParSubMesh::CreateFromDomain(*pp_pmesh, domain_attributes)
  );

  // Define the L2 finite element spaces representing currents in the positive/negative current collector
  L2_FECollection l2_fe_coll(0, dc_pmesh_currents.Dimension());

  // Create FiniteElementSpace for currents in the positive/negative current collector
  ParFiniteElementSpace * current_l2space;
  current_l2space = new ParFiniteElementSpace(&dc_pmesh_currents, &l2_fe_coll);
  
  HYPRE_BigInt fe_msmd_size_owned = 
  current_collector_positive_h1space->GetTrueVSize() + 
  current_collector_negative_h1space->GetTrueVSize() + 
  current_l2space->GetTrueVSize() + 
  1;
  // The last one is for the scaler dof for the reference potential 

  mfem::BlockVector y;
  CurrentCollectorOperator current_oper(current_collector_positive_h1space, current_collector_negative_h1space, current_l2space, fe_msmd_size_owned, y);

  bool last_step = false;
  for (int ti = 1; !last_step; ti++)
  {
    ode_solver->Step(x, t, dt);
    real_t V = oper.GetVoltage();

    last_step = (t_final >= 0 && t + dt / 2 >= t_final) || V <= CELL->lvoff() || V >= CELL->uvoff();

    // Print the time, voltage and SoC to the screen
    if (output_steps && ti == 1 && Mpi::Root())
      std::cout << "step\ttime[s]\tvoltage[V]\tSoC[%]" << std::endl;

    if (output_steps && (last_step || (ti % output_steps) == 0))
    {
      real_t SoC = oper.GetSoC();
      if (Mpi::Root())
      {
        std::cout << std::left << ti << "\t" << t << "\t";
        std::cout << std::setprecision(8) << std::setw(8) << V << "\t" << SoC << std::endl;
      }
    }
  }

  return 0;
}
