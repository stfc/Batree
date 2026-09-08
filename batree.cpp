// Batree - An MFEM-based SPM, SPMe and P2D solver

#include "operators/EChemOperator.hpp"

int
main(int argc, char * argv[])
{
  // Initialize MPI and HYPRE.
  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();

  // Parse command-line options.
  std::string model = "SPM";
  std::string cell = "LGM50";
  mfem::real_t c_rate = 1;
  int order = 1;
  int ode_solver_type = 21;
  mfem::real_t t_final = -1.0;
  mfem::real_t dt = 1.0;
  int output_steps = 5;

  mfem::OptionsParser args(argc, argv);
  args.AddOption(&model, "-m", "--model", "Electrochemical model: SPM, SPMe, or P2D.");
  args.AddOption(&cell, "-c", "--cell", "Cell model: LGM50 or Enertech.");
  args.AddOption(&c_rate, "-cr", "--c-rate", "The C-rate to run a constant current (dis)charge");
  args.AddOption(&order, "-o", "--order", "Order (degree) of the finite elements.");
  args.AddOption(&ode_solver_type, "-s", "--ode-solver", mfem::ODESolver::Types.c_str());
  args.AddOption(&t_final, "-tf", "--t-final", "Final time; start time is 0.");
  args.AddOption(&dt, "-dt", "--time-step", "Time step.");
  args.AddOption(&output_steps, "-os", "--output-steps", "Output every n-th timestep.");
  args.Parse();
  if (!args.Good())
  {
    args.PrintUsage(std::cout);
    return 1;
  }

  if (mfem::Mpi::Root())
    args.PrintOptions(std::cout);

  // Define the ODE solver used for time integration. Several implicit
  // singly diagonal implicit Runge-Kutta (SDIRK) methods, as well as
  // explicit Runge-Kutta methods are available in MFEM. For now, we only
  // support implicit methods and we have only tested Backward Euler.
  std::unique_ptr<mfem::ODESolver> ode_solver = mfem::ODESolver::Select(ode_solver_type);

  // Initialise properties dependent on the electrochemical model, cell type, current and FE order
  init_settings(model, cell, c_rate, order);

  // Initialize the ElectroChemistry operator.
  mfem::real_t t = 0.0;
  mfem::BlockVector x;
  EChemOperator oper(order, x);

  // Perform time-integration (looping over the time iterations, ti, with a
  // time-step dt).
  ode_solver->Init(oper);

  bool last_step = false;
  for (int ti = 1; !last_step; ti++)
  {
    ode_solver->Step(x, t, dt);
    mfem::real_t V = oper.GetVoltage();

    last_step = (t_final >= 0 && t + dt / 2 >= t_final) || V <= CELL->lvoff() || V >= CELL->uvoff();

    // Print the time, voltage and SoC to the screen
    if (output_steps && ti == 1 && mfem::Mpi::Root())
      std::cout << "step\ttime[s]\tvoltage[V]\tSoC[%]" << std::endl;

    if (output_steps && (last_step || (ti % output_steps) == 0))
    {
      mfem::real_t SoC = oper.GetSoC();
      if (mfem::Mpi::Root())
      {
        std::cout << std::left << ti << "\t" << t << "\t";
        std::cout << std::setprecision(8) << std::setw(8) << V << "\t" << SoC << std::endl;
      }
    }
  }

  return 0;
}
