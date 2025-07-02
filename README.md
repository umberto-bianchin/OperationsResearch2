# Operations Research 2 – Traveling Salesman Problem (TSP)

**Course:** Operations Research 2
**Professor:** Matteo Fischetti [matteo.fischetti@unipd.it]
**University:** University of Padua – Department of Information Engineering
**Academic Year:** 2024–2025

**Authors:**
* Umberto Bianchin – [umberto.bianchin@studenti.unipd.it] 
* Francesco De Nicola – [francesco.denicola@studenti.unipd.it]

---

## 🎯 Project Overview

This project presents a complete experimental framework for solving the **Traveling Salesman Problem (TSP)** using multiple algorithmic families:

* **Heuristics:** Nearest Neighbor, Extra-Mileage
* **Metaheuristics:** GRASP, VNS, Tabu Search, Genetic Algorithm
* **Exact Methods:** Benders Decomposition, Branch-and-Cut with CPLEX and Concorde
* **Matheuristics:** Hard Fixing, Local Branching

Each method is implemented in C and shares unified data structures, solution trackers, and input parsers. We also provide a **performance profiler** that benchmarks all algorithms developed by prof.
Domenico Salvagnin from the University of Padua called perfprof.py.

## ⚙️ Usage

## Requirements

- **C Compiler**: GCC or MSVC
- **Python 3** with `matplotlib` and `numpy`
- **IBM ILOG CPLEX** (academic license available)
- **Concorde TSP Solver** for enhanced exact solving

## 📦 Command Line Options

To run the program, you can pass various options like this:

```bash
./main -file <filename> -a <algorithm> -t <time_limit> ...
```

The following options are supported (aliases are accepted as shown):

### 📂 Input

* `-file`, `-input`, `-f <filename>`
  Specifies the `.tsp` input file containing the TSP instance (TSPLIB format).

* `-n`, `-nodes <int>`
  Generates a random instance with the given number of nodes.

### 🧠 Algorithm Selection

* `-a`, `-algorithm <char>`
  Selects the algorithm to use. Valid codes:

  * `N`: Nearest Neighbor
  * `E`: Extra-Mileage
  * `V`: Variable Neighborhood Search (VNS)
  * `T`: Tabu Search
  * `G`: GRASP
  * `Z`: Genetic Algorithm
  * `B`: Benders Decomposition
  * `C`: Branch-and-Cut
  * `H`: Hard Fixing
  * `L`: Local Branching

### ⚙️ General Configuration

* `-t`, `-time_limit <float>`
  Sets the maximum running time in seconds.

* `-seed <int>`
  Seed value for the random number generator.

* `-i <0|1>`
  Enables (`1`) or disables (`0`) integer edge costs.

* `-r <char>`
  Sets the running mode:

  * `n`: normal
  * `b`: benchmark
  * `t`: test

### 🔧 Algorithm-Specific Parameters

* `-kick <int>`
  Number of perturbation moves (e.g., in VNS or GRASP).

* `-kopt <int>`
  Value of `k` in the k-OPT move (must be ≥ 3 for VNS).

* `-alpha <int>`
  GRASP randomness parameter (0–100).

* `-minc <int>`
  Number of candidate nodes considered in GRASP (MIN\_COSTS).

* `-maxt <int>`
  Maximum tabu tenure in Tabu Search.

* `-mint <int>`
  Minimum tabu tenure in Tabu Search.

* `-stept <int>`
  Tenure increment step in Tabu Search.

* `-warmup <0|1>`
  Enables warm-starting with a heuristic solution.

* `-posting <0|1>`
  Enables posting heuristic solutions during Branch-and-Cut.

* `-concorde <0|1>`
  Enables use of Concorde for fractional cut separation.

* `-depth <int>`
  Maximum node depth for posting heuristic solutions.

* `-fixedprob <0|1>`
  In Hard Fixing, whether to use a fixed probability (`1`) or decreasing schedule (`0`).

* `-probability <int>`
  Edge-fixing probability for Hard Fixing (up to 90).

* `-klocal <int>`
  Local Branching radius (number of edges allowed to change).

* `-cdepth <int>`
  CPLEX node limit for each iteration of matheuristics.

* `-population <int>`
  Population size for the Genetic Algorithm.

* `-generation <int>`
  Number of children generated per generation in GA.

### 🆘 Help

* `-help`, `--help`
  Displays the help message and usage instructions.

---

## 📊 Performance Profiling

Use the C implementation to populate CSV logs and `perfprof.py` (Python) to visualize:

```bash
python3 perfprof.py results.csv output.pdf
```

---

## 📌 Installation – IBM ILOG CPLEX

See full installation instructions in the report (Appendix A), or:

* **Windows:** Add CPLEX and CONCERT paths to the `Path` system variable.
* **macOS:** Export `CPLEX_HOME` and extend `PATH`.
* **VS Code:** Use `c_cpp_properties.json` for header inclusion.

---

## 🧪 Evaluation & Results

We tested all algorithms on synthetic TSP instances of 300 to 1000 nodes with various time limits. The **best-performing** methods (in terms of cost and time) were:

* 🥇 **Hard Fixing (Matheuristic)**
* 🥈 **Variable Neighborhood Search (Metaheuristic)**
* 🥉 **Nearest Neighbor + 2-OPT (Heuristic)**

---

## 📎 Links

🔗 [GitHub Repository](https://github.com/umberto-bianchin/OperationsResearch2)

📄 [📄 Download the full report (PDF)](./Operations_Research_2.pdf) for in-depth explanations, comparisons, and performance plots.
