# Transmon Coupled to a Resonator

This project investigates the energy spectrum of a Transmon qubit capacitively coupled to a microwave resonator. Its main goals are:

- to construct the bare Transmon Hamiltonian and calculate $f_{01}$, $f_{12}$, and the anharmonicity $\alpha$;
- to model the coupled system using the Duffing approximation and the rotating-wave approximation;
- to observe the avoided crossing produced when the Transmon frequency approaches the resonator frequency;
- to compare the simplified Duffing model with a numerical `scqubits` model using the full Transmon Hamiltonian.

---

## Part 1A: Bare Transmon

First, we have to know the $f_{01}$ (qubit frequency), $f_{12}$ and $\alpha$ of a transmon qubit. 

Thus, we will first build the Hamiltonian of a transmon qubit to solve for the parameters mentioned above.

Its Hamiltonian is

$$\hat{H_{\mathrm{T}}}=4E_C(\hat{n}-n_g)^2-E_J\cos(\hat{\phi}).$$

The physical quantities are:

- $E_C$: charging energy
- $E_J$: Josephson energy
- $n_g$: offset charge
- $\hat{n}$: Cooper-pair number operator
- $\hat{\phi}$: superconducting phase operator

The phase and Cooper-pair number operators are conjugate variables satisfying

$$[\hat{\phi},\hat{n}]=i.$$

The charging-energy term is

$$\hat{H}_C=4E_C(\hat{n}-n_g)^2,$$

while the Josephson-energy term is

$$\hat{H}_J=-E_J\cos(\hat{\phi}).$$

Here, we choose the charge basis to construct the matrix of the Hamiltonian, since

$$\hat{n} \lvert n \rangle=n \vert n \rangle.$$

Thus,

$$\langle m\vert H_C\vert n \rangle=4E_C(n-n_g)^2\delta_{mn},$$

the charging-energy contribution is diagonal.

For the Josephson contribution, since

$$\cos \hat{\phi} = \frac{e^{i\hat{\phi}}+e^{-i\hat{\phi}}}{2}$$

and

$$e^{i\hat{\phi}}\lvert n \rangle = \lvert n+1 \rangle, e^{-i\hat{\phi}}\lvert n \rangle = \lvert n-1 \rangle$$

so

$$H_J\lvert n\rangle = -\frac{E_J}{2} (\lvert n+1 \rangle + \lvert n-1 \rangle).$$

Hence

$$\langle n+1\vert H_J\vert n \rangle=\langle n\vert H_J\vert n+1 \rangle = -\frac{E_J}{2},$$

which is off-diagonal.

The complete Hamiltonian of transmon is

$$ H_{\mathrm{T}} =
\begin{pmatrix}
4E_C(n_1-n_g)^2 & -\dfrac{E_J}{2} & 0 & \cdots \\
-\dfrac{E_J}{2} & 4E_C(n_2-n_g)^2 & -\dfrac{E_J}{2} & \cdots \\
0 & -\dfrac{E_J}{2} & 4E_C(n_3-n_g)^2 & \cdots \\
\vdots & \vdots & \vdots & \ddots
\end{pmatrix}.
$$

So after solving $H_{\mathrm{T}}\lvert n\rangle =E_n\lvert n\rangle$, we can define

$$f_{01}=E_1-E_0,$$

$$f_{12}=E_2-E_1,$$

$$\alpha = f_{12}-f_{01}.$$

So if $E_J$, $E_C$ and $n_g$ are given, we can know $H_{\mathrm{T}}$, so we can solve the eigenvalue problem and get $f_{01}$ (qubit frequency), $f_{12}$ and $\alpha$ of a transmon qubit.


### Numerical implementation

The charge basis is truncated to

$n=-N_{\mathrm{cut}},\ldots,N_{\mathrm{cut}},$

with $N_{\mathrm{cut}}=15$. Therefore, the numerical Hamiltonian has dimension

$2N_{\mathrm{cut}}+1=31.$

The parameters used in Part 1A are:

| Parameter | Value |
|---|---:|
| $E_C$ | $0.3$ GHz |
| $E_J/E_C$ | $50$ |
| $E_J$ | $15$ GHz |
| $n_g$ | $0$ |
| $N_{\mathrm{cut}}$ | $15$ |

The charging term is constructed as a diagonal NumPy matrix. The Josephson term is constructed by placing $-E_J/2$ directly above and below the main diagonal. The resulting matrix is converted into a QuTiP \`Qobj\` and numerically diagonalized.

The program then calculates

$f_{01}=E_1-E_0,$

$f_{12}=E_2-E_1,$

and

$\alpha=f_{12}-f_{01}.$

For the parameters above, the program obtains approximately

| Quantity | Numerical result |
|---|---:|
| $f_{01}$ | $5.683$ GHz |
| $f_{12}$ | $5.338$ GHz |
| $\alpha$ | $-0.345$ GHz |

The resonator frequency is initially set equal to the calculated Transmon frequency,

$f_r=f_{01},$

so the uncoupled Transmon and resonator begin at resonance.

---

## Part 1B: Duffing Transmon–Resonator Model

We can do a Taylor expansion of $\cos \hat{\phi}$ on $H_{\mathrm{T}}$ to simplify the Hamiltonian, we then have

$$ H_{\mathrm{T}} \approx 4E_C\hat{n}^{2} + \frac{E_J}{2}\hat{\phi}^{2} - \frac{E_J}{24}\hat{\phi}^{4},$$

note that we have dropped the constant term.

### Ladder-operator representation

We introduce the Transmon annihilation and creation operators $b$ and
$b^\dagger$ through

$$\hat{\phi}=\phi_{\mathrm{zpf}}(b+b^\dagger),$$

$$\hat{n}=i n_{\mathrm{zpf}}(b^\dagger-b),$$

where

$$\phi_{\mathrm{zpf}}=
\left(\frac{2E_C}{E_J}\right)^{1/4},$$

$$n_{\mathrm{zpf}}=
\left(\frac{E_J}{32E_C}\right)^{1/4}.
$$

These definitions satisfy

$$[\hat{\phi},\hat{n}]=i$$

because

$$[b,b^\dagger]=1.$$

Substituting these operators into the quadratic part of the Hamiltonian
gives

$$4E_C\hat{n}^2+\frac{E_J}{2}\hat{\phi}^2=
\sqrt{8E_JE_C}
\left(
b^\dagger b+\frac{1}{2}
\right).
$$

Therefore, the harmonic plasma frequency is

$$
f_p=\sqrt{8E_JE_C},
$$

when all energies are expressed in frequency units.

The quartic contribution becomes

$$-\frac{E_J}{24}\hat{\phi}^4=
-\frac{E_C}{12}(b+b^\dagger)^4,
$$

since

$$
\phi_{\mathrm{zpf}}^4=
\frac{2E_C}{E_J}.
$$

Under the rotating-wave approximation, only the terms that preserve the
number of Transmon excitations are retained:

$$
(b+b^\dagger)^4
\longrightarrow
6b^{\dagger 2}b^2+12b^\dagger b+3.
$$

Thus,

$$
-\frac{E_C}{12}(b+b^\dagger)^4
\longrightarrow
-\frac{E_C}{2}b^{\dagger 2}b^2
-E_Cb^\dagger b
-\frac{E_C}{4}.
$$

After dropping another constant term, the low-energy Transmon
Hamiltonian becomes

$$
H_q=
\left(
\sqrt{8E_JE_C}-E_C
\right)b^\dagger b
-\frac{E_C}{2}b^{\dagger 2}b^2.
$$

Defining

$$
f_q\approx\sqrt{8E_JE_C}-E_C
$$

and

$$
\alpha\approx-E_C,
$$

we obtain the Duffing-oscillator Hamiltonian

$$
\boxed{
H_q=
f_qb^\dagger b
+
\frac{\alpha}{2}b^{\dagger 2}b^2
}.
$$

Here, $f_q$ is the approximate $0\rightarrow1$ transition frequency and
$\alpha$ is the Transmon anharmonicity.

### Energy levels of the Duffing oscillator

The number states satisfy

$$
b^\dagger b\lvert m\rangle=m\lvert m\rangle.
$$

Also,

$$
b^{\dagger 2}b^2\lvert m\rangle=
m(m-1)\lvert m\rangle.
$$

Therefore,

$$
H_q\lvert m\rangle=
\left[
mf_q+\frac{\alpha}{2}m(m-1)
\right]\lvert m\rangle,
$$

and the energy of the $m$-th level is

$$
\boxed{
E_m=
mf_q+\frac{\alpha}{2}m(m-1)
}.
$$

For the three lowest levels,

$$
E_0=0,
$$

$$
E_1=f_q,
$$

$$
E_2=2f_q+\alpha.
$$

Consequently,

$$
f_{01}=E_1-E_0=f_q,
$$

$$
f_{12}=E_2-E_1=f_q+\alpha,
$$

and hence

$$
\alpha=f_{12}-f_{01}.
$$

### Resonator Hamiltonian

A microwave resonator can be modeled as an LC harmonic oscillator:

$$
H_r=
\frac{\hat{Q}_r^2}{2C_r}
+
\frac{\hat{\Phi}_r^2}{2L_r},
$$

where $\hat{Q}_r$ is the charge operator and $\hat{\Phi}_r$ is the flux
operator of the resonator.

Its angular frequency is

$$
\omega_r=\frac{1}{\sqrt{L_rC_r}}.
$$

After quantization, its Hamiltonian becomes

$$
H_r=
\hbar\omega_r
\left(
a^\dagger a+\frac{1}{2}
\right),
$$

where $a^\dagger$ and $a$ create and annihilate one resonator photon.

After dropping the zero-point energy and expressing the Hamiltonian in
frequency units, we obtain

$$
\boxed{
H_r=f_ra^\dagger a
}.
$$

### Transmon–resonator coupling

Capacitive coupling connects the electric field of the resonator to the
charge of the Transmon. Its basic form can be written as

$$
H_{\mathrm{int}}\propto\hat{n}\hat{Q}_r.
$$

After quantization, the charge operators can be expressed as

$$
\hat{n}\propto i(b^\dagger-b),
\qquad
\hat{Q}_r\propto i(a^\dagger-a).
$$

After an equivalent redefinition of the phases of the ladder operators,
the interaction Hamiltonian can be written as

$$
H_{\mathrm{int}}=
g(a+a^\dagger)(b+b^\dagger).
$$

Expanding it gives

$$
H_{\mathrm{int}}=
g\left(
ab+ab^\dagger+a^\dagger b+a^\dagger b^\dagger
\right).
$$

To determine which terms are important, we move to the interaction
picture. The ladder operators evolve as

$$
a(t)=ae^{-if_rt},
\qquad
a^\dagger(t)=a^\dagger e^{if_rt},
$$

$$
b(t)=be^{-if_qt},
\qquad
b^\dagger(t)=b^\dagger e^{if_qt}.
$$

Therefore,

$$
ab\longrightarrow ab\,e^{-i(f_r+f_q)t},
$$

$$
a^\dagger b^\dagger
\longrightarrow
a^\dagger b^\dagger e^{i(f_r+f_q)t},
$$

$$
a^\dagger b
\longrightarrow
a^\dagger b\,e^{i(f_r-f_q)t},
$$

$$
ab^\dagger
\longrightarrow
ab^\dagger e^{-i(f_r-f_q)t}.
$$

When the Transmon and resonator are close to resonance,

$$
f_q\approx f_r,
$$

the terms $ab$ and $a^\dagger b^\dagger$ oscillate rapidly at approximately
$f_r+f_q$ and average to zero. The terms $a^\dagger b$ and $ab^\dagger$
vary slowly because they depend on the detuning

$$
\Delta=f_q-f_r.
$$

The rotating-wave approximation therefore gives

$$
\boxed{
H_{\mathrm{int}}
\approx
g(a^\dagger b+ab^\dagger)
}.
$$

The term $a^\dagger b$ converts a qubit excitation into a resonator photon,
whereas $ab^\dagger$ converts a resonator photon into a qubit excitation.
Therefore, the interaction exchanges an excitation between the Transmon
and the resonator while conserving the total number of excitations.

The complete coupled Hamiltonian is

$$
\boxed{
H=
f_qb^\dagger b
+
\frac{\alpha}{2}b^{\dagger 2}b^2
+
f_ra^\dagger a
+
g(a^\dagger b+ab^\dagger)
}.
$$

This is the Duffing version of the Jaynes–Cummings model.


### Numerical implementation of the Duffing model

The program represents the resonator and Transmon in a composite Hilbert space. The truncation dimensions are

| Subsystem | Truncation |
|---|---:|
| Resonator | $N_r=10$ photon states |
| Transmon | $N_q=3$ energy levels |

The operators acting on the composite Hilbert space are

$a=a_r\otimes I_q,$

$b=I_r\otimes b_q.$

Their number operators are

$n_r=a^\dagger a,$

$n_q=b^\dagger b.$

The Hamiltonian implemented in the function \`coupled_H\` is

$H_r=f_rn_r,$

$H_q=f_qn_q+\frac{\alpha}{2}b^{\dagger2}b^2,$

$H_{\mathrm{int}}=g(a^\dagger b+ab^\dagger),$

and

$H_{\mathrm{total}}=H_r+H_q+H_{\mathrm{int}}.$

The coupling strength is set to

$g=0.05\ \mathrm{GHz}.$

To reveal the avoided crossing, the simplified Transmon frequency is swept directly over the interval

$f_q\in[f_{01}-0.5,\ f_{01}+0.5]\ \mathrm{GHz},$

using 200 points, while the resonator frequency $f_r$ remains fixed.

At every value of $f_q$, the total Hamiltonian is diagonalized. The ground-state energy is subtracted from all eigenenergies:

$E_j^{\mathrm{relative}}=E_j-E_0.$

Only the first and second excited dressed levels are plotted. Far from resonance, these levels resemble the uncoupled Transmon and resonator states. Near $f_q=f_r$, the coupling mixes the two states and produces an avoided crossing.

---

## Part 1C: Transmon–Resonator Avoided Crossing

We now use the coupled Transmon–resonator Hamiltonian derived in Part 1B:

$$H=f_qb^\dagger b+\frac{\alpha}{2}b^{\dagger 2}b^2+f_ra^\dagger a+g(a^\dagger b+ab^\dagger).$$

The total excitation-number operator is

$$\hat{N}=a^\dagger a+b^\dagger b.$$

Under the rotating-wave approximation,

$$[H,\hat{N}]=0.$$

Therefore, the Hamiltonian can be separated into subspaces with fixed
total excitation number.

### Single-excitation subspace

In the single-excitation subspace,

$$N=1,$$

there are only two bare basis states:

$$\lvert g,1\rangle,\qquad \lvert e,0\rangle.$$

Their physical meanings are:

- $\lvert g,1\rangle$: the Transmon is in its ground state and the resonator contains one photon.
- $\lvert e,0\rangle$: the Transmon is in its first excited state and the resonator contains no photon.

The Hamiltonian matrix in this basis is

$$H_1=\left(\begin{array}{cc}
\langle g,1\vert H\vert g,1\rangle &
\langle g,1\vert H\vert e,0\rangle \\
\langle e,0\vert H\vert g,1\rangle &
\langle e,0\vert H\vert e,0\rangle
\end{array}\right).$$

Now, let's calculate each element.

First, we have to know that the ladder operators satisfy

$$a\lvert n_r\rangle=\sqrt{n_r}\lvert n_r-1\rangle,\qquad a^\dagger\lvert n_r\rangle=\sqrt{n_r+1}\lvert n_r+1\rangle,$$

and

$$b\lvert m_q\rangle=\sqrt{m_q}\lvert m_q-1\rangle,\qquad b^\dagger\lvert m_q\rangle=\sqrt{m_q+1}\lvert m_q+1\rangle.$$

#### 1. First diagonal element

For $\lvert g,1\rangle$, the Transmon excitation number is zero and the resonator photon number is one. Therefore,

$$\langle g,1\vert H\vert g,1\rangle=f_q(0)+\frac{\alpha}{2}(0)(-1)+f_r(1)+0=f_r.$$

Thus,

$$\boxed{\langle g,1\vert H\vert g,1\rangle=f_r.}$$

#### 2. Second diagonal element

For $\lvert e,0\rangle$, the Transmon excitation number is one and the resonator photon number is zero. Therefore,

$$\langle e,0\vert H\vert e,0\rangle=f_q(1)+\frac{\alpha}{2}(1)(1-1)+f_r(0)+0=f_q.$$

Thus,

$$\boxed{\langle e,0\vert H\vert e,0\rangle=f_q.}$$

The interaction contributes zero to both diagonal elements because it transfers an excitation between the two systems and therefore changes the original bare state.

#### 3. Upper-right off-diagonal element

The upper-right element is

$$\langle g,1\vert H\vert e,0\rangle=\langle g,1\vert H_{\mathrm{int}}\vert e,0\rangle.$$

Using

$$H_{\mathrm{int}}=g(a^\dagger b+ab^\dagger),$$

we obtain

$$\langle g,1\vert H_{\mathrm{int}}\vert e,0\rangle=g\left[\langle1\vert a^\dagger\vert0\rangle\langle g\vert b\vert e\rangle+\langle1\vert a\vert0\rangle\langle g\vert b^\dagger\vert e\rangle\right].$$

The required matrix elements are

$$\langle1\vert a^\dagger\vert0\rangle=1,\qquad \langle g\vert b\vert e\rangle=1,$$

while

$$\langle1\vert a\vert0\rangle=0,\qquad \langle g\vert b^\dagger\vert e\rangle=0.$$

Therefore,

$$\langle g,1\vert H_{\mathrm{int}}\vert e,0\rangle=g\left[(1)(1)+(0)(0)\right]=g.$$

Thus,

$$\boxed{\langle g,1\vert H\vert e,0\rangle=g.}$$

#### 4. Lower-left off-diagonal element

Similarly,

$$\langle e,0\vert H\vert g,1\rangle=\langle e,0\vert H_{\mathrm{int}}\vert g,1\rangle.$$

Expanding the interaction gives

$$\langle e,0\vert H_{\mathrm{int}}\vert g,1\rangle=g\left[\langle0\vert a^\dagger\vert1\rangle\langle e\vert b\vert g\rangle+\langle0\vert a\vert1\rangle\langle e\vert b^\dagger\vert g\rangle\right].$$

Here,

$$\langle0\vert a^\dagger\vert1\rangle=0,\qquad \langle e\vert b\vert g\rangle=0,$$

and

$$\langle0\vert a\vert1\rangle=1,\qquad \langle e\vert b^\dagger\vert g\rangle=1.$$

Therefore,

$$\langle e,0\vert H_{\mathrm{int}}\vert g,1\rangle=g\left[(0)(0)+(1)(1)\right]=g.$$

Thus,

$$\boxed{\langle e,0\vert H\vert g,1\rangle=g.}$$

Combining the four matrix elements gives

$$H_1=\left(\begin{array}{cc}
f_r & g \\
g & f_q
\end{array}\right).$$

### Eigenenergies

The eigenenergies satisfy

$$\det(H_1-EI)=0.$$

Therefore,

$$\det(H_1-EI)=
\left|\begin{array}{cc}
f_r-E & g \\
g & f_q-E
\end{array}\right|
=0.$$

which gives

$$(f_r-E)(f_q-E)-g^2=0.$$

After expansion,

$$E^2-(f_r+f_q)E+f_rf_q-g^2=0.$$

The two eigenenergies are

$$E_\pm=\frac{f_r+f_q}{2}\pm\frac{1}{2}\sqrt{(f_q-f_r)^2+4g^2}.$$

Defining the detuning as

$$\Delta=f_q-f_r,$$

we obtain

$$E_\pm=\frac{f_r+f_q}{2}\pm\frac{1}{2}\sqrt{\Delta^2+4g^2}.$$

The separation between the two levels is

$$E_+-E_-=\sqrt{\Delta^2+4g^2}.$$

### Avoided crossing

Without coupling, $g=0$, the two bare energy levels cross when

$$f_q=f_r.$$

When $g\neq0$, the energy separation never becomes zero. Its minimum
occurs at resonance:

$$\Delta=0.$$

At resonance, the minimum gap is

$$E_+-E_-=2g.$$

Therefore, the two energy levels repel each other and form an avoided
crossing.

### Dressed states at resonance

At resonance, let

$$f_q=f_r=f_0.$$

The Hamiltonian becomes

$$H_1=\left(\begin{array}{cc}
f_0 & g \\
g & f_0
\end{array}\right).$$

Its eigenenergies are

$$E_+=f_0+g,\qquad E_-=f_0-g.$$

The corresponding dressed states are

$$\lvert+\rangle=\frac{\lvert g,1\rangle+\lvert e,0\rangle}{\sqrt{2}},$$

and

$$\lvert-\rangle=\frac{\lvert g,1\rangle-\lvert e,0\rangle}{\sqrt{2}}.$$

At exact resonance, both dressed states are equal mixtures of the two
bare states:

$$\left|\langle g,1\vert+\rangle\right|^2=\left|\langle e,0\vert+\rangle\right|^2=\frac{1}{2},$$

$$\left|\langle g,1\vert-\rangle\right|^2=\left|\langle e,0\vert-\rangle\right|^2=\frac{1}{2}.$$

The dressed states are exactly $50/50$ mixtures only at exact resonance.

### Connection to the numerical model

In the numerical calculation, $E_J$ is varied to tune the Transmon
frequency $f_q$ through the fixed resonator frequency $f_r$. The
eigenenergies of the complete coupled Hamiltonian are calculated at each
value of $E_J$.

For the simplified interaction

$$H_{\mathrm{int}}=g(a^\dagger b+ab^\dagger),$$

the minimum gap is

$$\text{minimum gap}=2g.$$

However, the `scqubits` model uses

$$H_{\mathrm{int}}=g\hat{n}(a+a^\dagger).$$

Therefore, the effective coupling between $\lvert e,0\rangle$ and
$\lvert g,1\rangle$ is

$$g_{\mathrm{eff}}=\left|\langle g,1\vert g\hat{n}(a+a^\dagger)\vert e,0\rangle\right|.$$

Evaluating the resonator matrix element gives

$$g_{\mathrm{eff}}=g\left|\langle0\vert\hat{n}\vert1\rangle\right|.$$

Thus, the numerical minimum gap is approximately

$$\text{minimum gap}\approx2g_{\mathrm{eff}}=2g\left|\langle0\vert\hat{n}\vert1\rangle\right|.$$

This is why the minimum gap obtained from the full `scqubits` model does
not necessarily equal exactly $2g$, where $g$ is the coupling coefficient
entered in the code.


### Full numerical validation with scqubits

The simplified Duffing calculation treats $f_q$, $\alpha$, and $g$ as effective parameters. Part 1C validates the avoided crossing using the full Transmon Hamiltonian provided by \`scqubits\`.

The numerical subsystems are constructed as

~~~python
transmon = scq.Transmon(
    EJ=EJ,
    EC=EC,
    ng=ng,
    ncut=40,
    truncated_dim=Nq
)

resonator = scq.Oscillator(
    E_osc=wr,
    truncated_dim=Nr
)
~~~

The two subsystems are combined into one Hilbert space:

~~~python
hilbertspace = scq.HilbertSpace([transmon, resonator])
~~~

The program adds capacitive coupling through the Transmon charge operator and the resonator field quadrature:

~~~python
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
~~~

Thus, the interaction used in this calculation is

$H_{\mathrm{int}}=g\hat n(a+a^\dagger).$

Unlike the simplified Duffing model, the effective coupling is determined by the charge matrix element

$g_{\mathrm{eff}}=g\left|\langle0|\hat n|1\rangle\right|.$

The Josephson energy is swept over

$E_J\in[12,18]\ \mathrm{GHz}$

using 300 points. At each value of $E_J$, the program performs the following steps:

1. Update the Transmon Josephson energy.
2. Diagonalize the bare Transmon Hamiltonian and calculate its $f_{01}$.
3. Diagonalize the complete Transmon–resonator Hamiltonian.
4. Subtract the dressed ground-state energy.
5. Store the first two excited dressed energies.

The gap between the two relevant dressed levels is calculated from

$\mathrm{gap}(E_J)=E_2(E_J)-E_1(E_J).$

The avoided-crossing point is identified numerically using

~~~python
gap_idx = np.argmin(gap_values)
~~~

For the parameters used in the program, the result is approximately

| Quantity | Numerical result |
|---|---:|
| Bare Transmon frequency at the minimum gap | $5.6806$ GHz |
| Resonator frequency | $5.683$ GHz |
| Minimum dressed-state gap | $108.8$ MHz |

The small difference between the bare Transmon frequency and the resonator frequency at the numerical minimum is caused by the finite parameter grid and by energy shifts from the full interaction.

Part 1B and Part 1C therefore demonstrate the same physical phenomenon at two different levels of modeling:

| Part | Model | Tuned parameter | Coupling |
|---|---|---|---|
| 1B | Duffing approximation | Effective frequency $f_q$ | $g(a^\dagger b+ab^\dagger)$ |
| 1C | Full \`scqubits\` Transmon | Josephson energy $E_J$ | $g\hat n(a+a^\dagger)$ |

The agreement between the two calculations confirms that the observed spectral splitting is produced by Transmon–resonator hybridization.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/EnochChuang/Transmon-Coupled-to-a-Resonator.git
```

Enter the project folder:

```bash
cd Transmon-Coupled-to-a-Resonator
```

Install the required Python packages:

```bash
pip install -r requirements.txt
```

The main dependencies are:

- NumPy
- Matplotlib
- QuTiP
- scqubits

---

## Usage

After adding the simulation script to the repository, run it with

```bash
python "github 2(2).py"
```

Parts 1A–1C print the bare Transmon parameters and the minimum avoided-crossing gap, and generate two avoided-crossing figures: one from the Duffing model and one from the full `scqubits` model.
