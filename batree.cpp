// Batree - An MFEM-based SPM, SPMe and P2D solver

#include "operators/EChemOperator.hpp"
#include "operators/CurrentCollectorOperator.hpp"
#include "operators/EChemCellModel.hpp"

#include <chrono>

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
  std::string cc_model = "nonlinear"; // local model: linear, nonlinear or echem (the -m model)
  std::string cc_slope = "exact";     // exact (Newton), frozen (chord) or none (report's Option 1)
  mfem::real_t cc_R0 = 1.0;           // ohmic resistance of the local model (scaled)
  mfem::real_t cc_j0 = 0.5;           // exchange current of the nonlinear local model (scaled)
  mfem::real_t cc_omega = 1.0;        // relaxation factor for the element currents
  int cc_refine = 0;                  // uniform refinements of the collector meshes
  bool cc_paraview = false;           // write the collector solution for ParaView
  mfem::real_t cc_tol = 1e-10;        // coupling tolerance on the element currents
  mfem::real_t cc_scl = 1e-3;         // P2D self-consistency threshold
  mfem::real_t cc_h = 1e-5;           // finite-difference step for R̄ (echem coupling)

  mfem::OptionsParser args(argc, argv);
  args.AddOption(&model, "-m", "--model", "Electrochemical model: SPM, SPMe, or P2D.");
  args.AddOption(&cell, "-c", "--cell", "Cell model: LGM50 or Enertech.");
  args.AddOption(&c_rate, "-cr", "--c-rate", "C-rate to run a constant current (dis)charge.");
  args.AddOption(&order, "-o", "--order", "Order (degree) of the finite elements.");
  args.AddOption(&ode_solver_type, "-s", "--ode-solver", mfem::ODESolver::Types.c_str());
  args.AddOption(&t_final, "-tf", "--t-final", "Final time; start time is 0.");
  args.AddOption(&dt, "-dt", "--time-step", "Time step.");
  args.AddOption(&output_steps, "-os", "--output-steps", "Output every n-th timestep.");
  args.AddOption(&cc_model,
                 "-ccm",
                 "--cc-model",
                 "Collector cell model: linear, nonlinear (stand-ins) or echem (the -m model).");
  args.AddOption(&cc_tol, "-cctol", "--cc-tol", "Coupling tolerance on the element currents.");
  args.AddOption(&cc_scl, "-ccscl", "--cc-scl", "P2D self-consistency threshold.");
  args.AddOption(&cc_h, "-cch", "--cc-h", "Finite-difference step for the slope R (echem).");
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
  oper.SetSelfConsistencyThreshold(cc_scl);

  // Perform time-integration (looping over the time iterations, ti, with a
  // time-step dt).
  ode_solver->Init(oper);

  // Current collectors (MSMD): the operator builds its own meshes and spaces. c_rate is the
  // total current through the tabs.
  mfem::BlockVector y;
  CurrentCollectorOperator current_oper(order, cc_refine, c_rate, y);

  if (cc_model == "echem")
  {
    // Every collector element runs its own cell (the -m model), and the collectors couple them
    // at every time step.
    MFEM_VERIFY(mfem::Mpi::WorldSize() == 1, "-ccm echem runs in serial only for now.");
    mfem::Vector area;
    current_oper.ElementAreas(area);
    EChemCellModel cells(oper, *ode_solver, x, t, dt, current_oper.NumElements(), cc_h);

    const auto start = std::chrono::steady_clock::now();
    int total_passes = 0, unconverged = 0;
    bool last_step = false;
    for (int ti = 1; !last_step; ti++)
    {
      const CouplingInfo info =
          current_oper.SolveCoupled(cells, cc_slope, cc_omega, cc_tol, 40, ti == 1);
      cells.Commit();
      t += dt;
      total_passes += info.passes;
      unconverged += !info.converged;

      const mfem::real_t V = phi_scale * current_oper.TerminalVoltage();
      last_step =
          (t_final >= 0 && t + dt / 2 >= t_final) || V <= CELL->lvoff() || V >= CELL->uvoff();

      if (output_steps && ti == 1)
        std::cout << "\nstep\ttime[s]\tvoltage[V]\tSoC[%]\tSoC_min\tSoC_max\tpasses\tconv"
                     "\tresidual\ti_min\ti_max\tR_min\tR_max\n";
      if (output_steps && (last_step || (ti % output_steps) == 0 || !info.converged))
      {
        mfem::real_t soc, soc_min, soc_max;
        cells.SoC(area, soc, soc_min, soc_max);
        std::printf("%d\t%g\t%.7f\t%.4f\t%.4f\t%.4f\t%d\t%d\t%.2e\t%.6f\t%.6f\t%.4f\t%.4f\n",
                    ti,
                    t,
                    V,
                    soc,
                    soc_min,
                    soc_max,
                    info.passes,
                    int(info.converged),
                    info.residual,
                    info.imin,
                    info.imax,
                    info.rmin,
                    info.rmax);
      }
    }
    const std::chrono::duration<double> wall = std::chrono::steady_clock::now() - start;
    std::printf("\nCoupled run: %g s simulated, %d coupling passes, %d unconverged steps, %ld cell "
                "steps, %.1f s wall\n",
                t,
                total_passes,
                unconverged,
                cells.Steps(),
                wall.count());
    if (cc_paraview)
      current_oper.SaveParaView(0, t);
    return 0;
  }

  // Couple the collectors to the local model at one instant (the P2D is not connected yet).
  // OCV = 4.2 V, the same value the constant voltage coefficient used before.
  LocalModel local_model(cc_model, 4.2 / phi_scale, cc_R0, cc_j0);
  current_oper.SolveCoupled(local_model, cc_slope, cc_omega, cc_tol);
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
