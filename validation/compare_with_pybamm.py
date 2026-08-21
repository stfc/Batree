import matplotlib.pyplot as plt
import pybamm
import subprocess

# Close all figures
plt.close("all")

mfem_executable = "./../batree"

PLOT_MFEM  = True
PLOT_PYBAMM = True

# Standard plot settings for all plots.
def easyplot(ax, x, y, colour, linestyle, label, marker_indx = 1):

    ax.plot(
        x, y,
        linestyle,
        color=colour,
        linewidth=1.3,
        label=label,
        markevery=marker_indx
    )

# Function to run Batree and extract results from output.
def run_batree(sim_type, cell):

    print("Running... Batree | " + sim_type.ljust(4) + " | " + cell.ljust(8))
    cmd = [mfem_executable, "-m", sim_type, "-c", cell]
    result = subprocess.run(cmd, capture_output=True, text=True)

    times = []
    voltages = []

    for line in result.stdout.splitlines():
        parts = line.split()

        if len(parts) == 4 and parts[0].isdigit():
            _, time, voltage, _ = parts
            times.append(float(time))
            voltages.append(float(voltage))

    return times, voltages


# Function to run PyBAMM and extract results from output.
def run_pybamm(model, cell):

    print("Running... PyBaMM | " + str(model.__class__.__name__).ljust(4) + " | " + cell.ljust(8))
    solver = pybamm.IDAKLUSolver(options={"dt_max" : 1})
    sim = pybamm.Simulation(model, parameter_values=pybamm.ParameterValues(cell), solver=solver)
    soln = sim.solve([0, 3600])

    time = soln["Time [s]"].entries
    voltage = soln["Voltage [V]"].entries

    return time, voltage

# Running simulations and plotting results.
def run_and_plot(ax, cell, sim_type, pybamm_model, colour):

    if PLOT_MFEM:
        time, voltage = run_batree(sim_type, cell)
        easyplot(ax, time, voltage, colour, ".", f"{sim_type} (Batree)", 10)

    if PLOT_PYBAMM:
        time, voltage = run_pybamm(pybamm_model, cell)
        easyplot(ax, time, voltage, colour, "-", f"{sim_type} (PyBaMM)")


cells = ["Chen2020", "Ai2020"]
simulations = [
    ("SPM",  pybamm.lithium_ion.SPM(),  (197 / 255,  27 / 255, 138 / 255)),
    ("SPMe", pybamm.lithium_ion.SPMe(), ( 44 / 255, 127 / 255, 184 / 255)),
    ("DFN",  pybamm.lithium_ion.DFN(),  ( 49 / 255, 163 / 255,  84 / 255)),
]

fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True, constrained_layout=True)

for i, cell in enumerate(cells):
    ax = axs[i]
    for sim_type, pybamm_model, colour in simulations:
        run_and_plot(ax, cell, sim_type, pybamm_model, colour)
    ax.set_title(cell)
    ax.set_xlim(0, 3600)
    ax.set_xticks(range(0, 3601, 400))
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Voltage [V]")
    ax.grid(True, which="major", linestyle="--", alpha=0.4)
    ax.legend(fontsize="small", loc="lower left")

plt.savefig("compared_with_pybamm.png", dpi=300)
