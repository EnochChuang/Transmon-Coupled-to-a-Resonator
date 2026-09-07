import numpy as np
import matplotlib.pyplot as plt
import scqubits as scq
import qutip as qt
import os

Nr = 10   # resonator photon number truncation: |0>, ..., |9>
Nq = 3    # transmon levels: |0>, |1>, |2>
g = 0.05 
EC = 0.3
ng = 0.0
N_cut = 15     # charge states: n = -15, ..., 15
n_values = np.arange(-N_cut, N_cut + 1)
dim = len(n_values)
ratio = 50.0  # EJ/EC 的值
EJ = ratio * EC

# Part 1A: Bare Transmon
def transmon_hamiltonian(EJ,ng):
    H_charge = np.diag(
        4 * EC * (n_values - ng)**2
    )
    H_josephson = np.zeros((dim, dim))
    for i in range(dim - 1):
        H_josephson[i, i + 1] = -EJ / 2
        H_josephson[i + 1, i] = -EJ / 2

    H = H_charge + H_josephson

    return qt.Qobj(H)

H = transmon_hamiltonian(EJ, ng)

eigenenergies, eigenstates = H.eigenstates()
f01 = eigenenergies[1] - eigenenergies[0]
f12 = eigenenergies[2] - eigenenergies[1]
wq = f01
wr = f01 # resonator frequency
alpha = f12 - f01
fq_list = np.linspace(
    wq - 0.5,
    wq + 0.5,
    200
)
print("\n=== Part 1A: Bare Transmon Parameters ===")
print(f"f01 = {f01:.3f} GHz")
print(f"f12 = {f12:.3f} GHz")
print(f"alpha = {alpha:.3f} GHz")
# qubit-resonator detuning
Delta_qr = wq - wr
print(f"Delta_qr = {Delta_qr:.3f} GHz")

# Resonator annihilation operator
a = qt.tensor(
    qt.destroy(Nr),
    qt.qeye(Nq)
)

# Transmon annihilation operator
b = qt.tensor(
    qt.qeye(Nr),
    qt.destroy(Nq)
)

# Number operators
n_r = a.dag() * a
n_q = b.dag() * b

# Part 1B: Coupled Transmon-Resonator
def coupled_H(wq, wr, alpha, g):
    # Resonator Hamiltonian
    H_r = wr * n_r

    # Transmon Hamiltonian
    # Duffing oscillator approximation
    H_q = (
        wq * n_q
        + alpha / 2
        * b.dag() * b.dag() * b * b
    )

    # interaction Hamiltonian
    H_int = g * (a.dag() * b + a * b.dag())

    # Total coupled Hamiltonian
    H_total_GHz = H_r + H_q + H_int

    return H_total_GHz

num_eigenlevels = 3
energy_levels = np.zeros((num_eigenlevels, len(fq_list)))
for i, fq in enumerate(fq_list):
    H_total = coupled_H(fq, wr, alpha, g)
    energies = H_total.eigenenergies()
    energies = energies - energies[0]
    energy_levels[:, i] = np.real(energies[:num_eigenlevels])

plt.figure(figsize=(8, 5))
colors = ["#1B9431", "#E67E22"]
for level in range(1, num_eigenlevels):
    plt.plot(
        fq_list,
        energy_levels[level],
        color=colors[level - 1],
        label=f"Level {level}"
    )
plt.axvline(
    wr,
    color="gray",
    linestyle="--",
    label=r"$f_q=f_r$"
)
plt.xlabel(r"Transmon frequency $f_q$ [GHz]")
plt.ylabel("Energy relative to ground state [GHz]")
plt.title("Duffing Model: Transmon–Resonator Avoided Crossing")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig(
    "part1b_duffing_avoided_crossing.svg",
    format="svg",
    bbox_inches="tight"
)

# Part 1C: Building transmon-resonator by scqubits
transmon = scq.Transmon(EJ=EJ, EC=EC, ng=ng, truncated_dim=Nq, ncut=40)
resonator = scq.Oscillator(E_osc=wr, truncated_dim=Nr)
hilbertspace = scq.HilbertSpace([transmon, resonator])
hilbertspace.add_interaction(
    g=0.05,
    op1=(transmon.n_operator(), transmon),
    op2=(
        resonator.creation_operator()
        + resonator.annihilation_operator(),
        resonator
    ),
    add_hc=False,
    id_str="capacitive_coupling"
)

EJ_values = np.linspace(12.0, 18.0, 301)
f01_values = []
energy_levels = []

for EJ in EJ_values:
    transmon.EJ = EJ

    bare_evals = transmon.eigenvals(evals_count=3)
    f01 = bare_evals[1] - bare_evals[0]
    f01_values.append(f01)

    dressed_evals = hilbertspace.eigenvals(evals_count=3)
    dressed_evals = dressed_evals - dressed_evals[0]
    energy_levels.append(dressed_evals)

f01_values = np.array(f01_values)
energy_levels = np.array(energy_levels)

gap_values = energy_levels[:, 2] - energy_levels[:, 1]
gap_idx = np.argmin(gap_values)

print("\n=== Part 1C: Transmon–Resonator Avoided Crossing ===")
print(f"f01 at minimum gap = {f01_values[gap_idx]:.6f} GHz")
print(f"Minimum gap = {gap_values[gap_idx] * 1000:.6f} MHz")

plt.figure(figsize=(8, 5))
plt.plot(
    f01_values,
    energy_levels[:, 1],
    color="#1B9431",
    label="Dressed level 1"
)
plt.plot(
    f01_values,
    energy_levels[:, 2],
    color="#E67E22",
    label="Dressed level 2"
)
plt.axvline(
    wr,
    color="gray",
    linestyle="--",
    label=r"$f_{01}=f_r$"
)
plt.xlabel(r"Transmon frequency $f_{01}$ [GHz]")
plt.ylabel("Energy relative to ground state [GHz]")
plt.title("Full scqubits Model: Transmon–Resonator Avoided Crossing")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(
    "part1c_scqubits_avoided_crossing.svg",
    format="svg",
    bbox_inches="tight"
)
print("Figures are saved in:", os.getcwd())
plt.show()