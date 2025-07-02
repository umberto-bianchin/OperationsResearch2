# Operations Research 2 – Traveling Salesman Problem (TSP)

**Course:** Operations Research 2
**Professor:** [Matteo Fischetti] (matteo.fischetti@unipd.it)
**University:** University of Padua – Department of Information Engineering
**Academic Year:** 2024–2025

**Authors:**
* Umberto Bianchin – [umberto.bianchin@studenti.unipd.it] (umberto.bianchin@studenti.unipd.it)
* Francesco De Nicola – [francesco.denicola@studenti.unipd.it] (francesco.denicola@studenti.unipd.it)

---

## 🎯 Project Overview

This project presents a complete experimental framework for solving the **Traveling Salesman Problem (TSP)** using multiple algorithmic families:

* **Heuristics:** Nearest Neighbor, Extra-Mileage
* **Metaheuristics:** GRASP, VNS, Tabu Search, Genetic Algorithm
* **Exact Methods:** Benders Decomposition, Branch-and-Cut with CPLEX and Concorde
* **Matheuristics:** Hard Fixing, Local Branching

Each method is implemented in C and shares unified data structures, solution trackers, and input parsers. We also provide a **performance profiler** that benchmarks all algorithms using Dolan–Moré profiles generated via a Python script.

---

## 📁 Project Structure
OperationsResearch2/
├── CMakeLists.txt        
├── README.md           
├── include/           # Header files
│   ├── benders.h
│   ├── data_struct.h
│   ├── heuristics.h
│   ├── matheuristics.h
│   ├── metaheuristics.h
│   ├── parsers.h
│   ├── utils.h
│   └── exact.h
├── src/               # Source code
│   ├── benders.c
│   ├── data_struct.c
│   ├── heuristics.c
│   ├── main.c
│   ├── matheuristics.c
│   ├── metaheuristics.c
│   ├── parsers.c
│   ├── utils.c
│   └── exact.c
└── results/           # Experimental results and CSV files

---

## ⚙️ Usage

## Requirements

- **C Compiler**: GCC or MSVC
- **Python 3** with `matplotlib` and `numpy`
- **IBM ILOG CPLEX** (academic license available)
- **Concorde TSP Solver** for enhanced exact solving
  
### 🔧 Compilation

```bash
make
```

---

### 🚀 Run Example

To run the solver on an instance:

```bash
./tsp -file tsp_instance.tsp -a V -i 1 -t 60 -seed 42 -r n
```

### 🧭 Command Line Options

| Option           | Description                                                                                                                                                     |                   |                                      |
| ---------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------- | ------------------------------------ |
| `-file <file>`   | Input `.tsp` file (TSPLIB format)                                                                                                                               |                   |                                      |
| `-n <nodes>`     | Number of randomly generated nodes                                                                                                                              |                   |                                      |
| `-a <algorithm>` | Algorithm: `N` (NN), `E` (EM), `V` (VNS), `T` (Tabu), `G` (GRASP), `Z` (Genetic), `B` (Benders), `C` (Branch-and-Cut), `H` (Hard Fixing), `L` (Local Branching) |                   |                                      |
| `-t <seconds>`   | Time limit                                                                                                                                                      |                   |                                      |
| \`-i \[0         | 1]\`                                                                                                                                                            | Use integer costs |                                      |
| `-seed <int>`    | Seed for random number generation                                                                                                                               |                   |                                      |
| \`-r \[n         | b                                                                                                                                                               | t]\`              | Run mode: normal, benchmark, or test |
| `-help`          | Show help message with full parameter list                                                                                                                      |                   |                                      |

**Example with metaheuristics tuning:**

```bash
./main -n 1000 -a V -t 60 -seed 123 -i 1 -kopt 5 -kick 7 -r b
```

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

## 📚 References

This project implements techniques described in:

* Dantzig, Fulkerson & Johnson (1954) – Cutting-plane approach
* Lin & Kernighan (1973) – k-OPT heuristics
* Glover (1990s) – Tabu Search
* Fischetti & Lodi (2003) – Local Branching
* Applegate et al. – Concorde TSP Solver
* IBM ILOG CPLEX – ILP Solver (used for all exact methods)

---

## 🧪 Evaluation & Results

We tested all algorithms on synthetic TSP instances of 300 to 1000 nodes with various time limits. The **best-performing** methods (in terms of cost and time) were:

* 🥇 **Hard Fixing (Matheuristic)**
* 🥈 **Variable Neighborhood Search (Metaheuristic)**
* 🥉 **Nearest Neighbor + 2-OPT (Heuristic)**

---

## 📎 Links

🔗 [GitHub Repository (if applicable)](https://github.com/umberto-bianchin/OperationsResearch2)

📄 Full report: See `Operations_Research_2.pdf` for in-depth explanations, comparisons, and performance plots.
