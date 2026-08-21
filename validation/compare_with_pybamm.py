import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pybamm
import subprocess

# Close all figures
plt.close("all")

batree_executable = "./../batree"

PLOT_BATREE = True
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
def run_batree(sim_type, cell, c_rate):

    print("Running... Batree | " + sim_type.ljust(4) + " | " + cell.ljust(8) + " | " + c_rate + "C")
    cmd = [batree_executable, "-m", sim_type, "-c", cell, "-cr", c_rate, "-tf", str(3600 / float(c_rate))]
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
def run_pybamm(model, cell, c_rate):

    print("Running... PyBaMM | " + str(model.__class__.__name__).ljust(4) + " | " + cell.ljust(8) + " | " + c_rate + "C")
    solver = pybamm.IDAKLUSolver(options={"dt_max" : 1})
    sim = pybamm.Simulation(model, parameter_values=pybamm.ParameterValues(cell), C_rate=float(c_rate), solver=solver)
    soln = sim.solve([0, 3600 / float(c_rate)])

    time = soln["Time [s]"].entries
    voltage = soln["Voltage [V]"].entries

    return time, voltage

# Running simulations and plotting results.
def run_and_plot(ax, sim_type, pybamm_model, cell, c_rate, colour):

    if PLOT_BATREE:
        time, voltage = run_batree(sim_type, cell, c_rate)
        easyplot(ax, time, voltage, colour, ".", f"{sim_type} (Batree)", 20)

    if PLOT_PYBAMM:
        time, voltage = run_pybamm(pybamm_model, cell, c_rate)
        easyplot(ax, time, voltage, colour, "-", f"{sim_type} (PyBaMM)")


cells = ["Chen2020", "Ai2020"]
models = [
    ("SPM",  pybamm.lithium_ion.SPM(),  ["#cbc9e2", "#bdd7e7", "#bae4b3"]),
    ("SPMe", pybamm.lithium_ion.SPMe(), ["#9e9ac8", "#6baed6", "#74c476"]),
    ("DFN",  pybamm.lithium_ion.DFN(),  ["#6a51a3", "#2171b5", "#238b45"]),
]
c_rates = ["0.7", "1", "2"]

fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True, constrained_layout=True)

for i, cell in enumerate(cells):
    ax = axs[i]
    for sim_type, pybamm_model, colours in models:
        for colour, c_rate in zip(colours, c_rates):
            run_and_plot(ax, sim_type, pybamm_model, cell, c_rate, colour)
    ax.set_title(cell)
    ax.set_xlim(0, 5400)
    ax.set_ylim(2.4, 4.2)
    ax.set_xticks(range(0, 5401, 600))
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Voltage [V]")
    ax.grid(True, which="major", linestyle="--", alpha=0.4)

method_legend = fig.legend(
    handles=[
        Line2D([], [], color="black", linewidth=1.3, marker=".", linestyle="None", label="Batree"),
        Line2D([], [], color="black", linewidth=1.3, linestyle="-", label="PyBaMM"),
    ],
    title="Implementation",
    loc="upper left",
    bbox_to_anchor=(1.01, 0.90),
)

rate_legend = fig.legend(
    handles=[
        Line2D([], [], color=colour, linewidth=1.3, label=f"{c_rate}C")
        for c_rate, colour in zip(c_rates, ["#6a51a3", "#2171b5", "#238b45"])
    ],
    title="Discharge rate",
    loc="upper left",
    bbox_to_anchor=(1.01, 0.75),
)

model_legend = fig.legend(
    handles=[
        Line2D([], [], color=colour, linewidth=1.3, label=sim_type)
        for (sim_type, _, _), colour in zip(models, ["#cccccc", "#969696", "#525252"])
    ],
    title="Model",
    loc="upper left",
    bbox_to_anchor=(1.01, 0.55),
)

plt.savefig("compared_with_pybamm.png", dpi=300, bbox_inches="tight")
