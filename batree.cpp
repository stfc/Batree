// Batree - An MFEM-based SPM, SPMe and P2D solver

#include "operators/EChemOperator.hpp"
#include "operators/CurrentCollectorOperator.hpp"

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

  // Current-collector coupling test: a local model stands in for the per-element P2D
  std::string cc_model = "nonlinear"; // local model: linear or nonlinear
  std::string cc_slope = "exact";     // exact (Newton), frozen (chord) or none (report's Option 1)
  mfem::real_t cc_R0 = 1.0;           // ohmic resistance of the local model (scaled)
  mfem::real_t cc_j0 = 0.5;           // exchange current of the nonlinear local model (scaled)
  mfem::real_t cc_omega = 1.0;        // relaxation factor for the element currents
  int cc_refine = 0;                  // uniform refinements of the collector meshes
  bool cc_paraview = false;           // write the collector solution for ParaView

  mfem::OptionsParser args(argc, argv);
  args.AddOption(&model, "-m", "--model", "Electrochemical model: SPM, SPMe, or P2D.");
  args.AddOption(&cell, "-c", "--cell", "Cell model: LGM50 or Enertech.");
  args.AddOption(&c_rate, "-cr", "--c-rate", "C-rate to run a constant current (dis)charge.");
  args.AddOption(&order, "-o", "--order", "Order (degree) of the finite elements.");
  args.AddOption(&ode_solver_type, "-s", "--ode-solver", mfem::ODESolver::Types.c_str());
  args.AddOption(&t_final, "-tf", "--t-final", "Final time; start time is 0.");
  args.AddOption(&dt, "-dt", "--time-step", "Time step.");
  args.AddOption(&output_steps, "-os", "--output-steps", "Output every n-th timestep.");
  args.AddOption(
      &cc_model, "-ccm", "--cc-model", "Collector test local model: linear or nonlinear.");
  args.AddOption(&cc_slope,
                 "-ccs",
                 "--cc-slope",
                 "Collector coupling slope: exact (Newton), frozen (chord) or none (Picard).");
  args.AddOption(&cc_R0, "-ccr", "--cc-resistance", "Local model ohmic resistance (scaled).");
  args.AddOption(&cc_j0, "-ccj", "--cc-j0", "Nonlinear local model exchange current (scaled).");
  args.AddOption(&cc_omega, "-ccw", "--cc-omega", "Relaxation factor for the element currents.");
  args.AddOption(
      &cc_refine, "-ccref", "--cc-refine", "Uniform refinements of the collector meshes.");
  args.AddOption(&cc_paraview,
                 "-ccpv",
                 "--cc-paraview",
                 "-no-ccpv",
                 "--no-cc-paraview",
                 "Write the collector solution (phi+, phi-, i) to ParaView/.");
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
  init_settings(model, cell, order);

  // Initialize the ElectroChemistry operator.
  mfem::real_t t = 0.0;
  mfem::BlockVector x;
  EChemOperator oper(order, x);

  // Perform time-integration (looping over the time iterations, ti, with a
  // time-step dt).
  ode_solver->Init(oper);

  // Current collectors (MSMD): the operator builds its own meshes and spaces. The collector tab
  // flux uses the external current I, which init_settings does not set, so set it first.
  oper.SetCurrent(c_rate);
  mfem::BlockVector y;
  CurrentCollectorOperator current_oper(order, cc_refine, y);

  // Couple the collectors to the local model at one instant (the P2D is not connected yet).
  // OCV = 4.2 V, the same value the constant voltage coefficient used before.
  LocalModel local_model(cc_model, 4.2 / phi_scale, cc_R0, cc_j0);
  current_oper.SolveCoupled(local_model, cc_slope, cc_omega);
  if (cc_paraview)
    current_oper.SaveParaView(0, 0.0);

  bool last_step = false;
  for (int ti = 1; !last_step; ti++)
  {
    oper.SetCurrent(c_rate);
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
        std::cout << std::setw(8) << std::setprecision(7) << V << "\t";
        std::cout << std::setprecision(4) << SoC << std::endl;
      }
    }
  }

  return 0;
}
