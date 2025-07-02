# Operations Research 2 – Traveling Salesman Problem (TSP)

**Course:** Operations Research 2\
**Professor:** Matteo Fischetti [matteo.fischetti@unipd.it]\
**University:** University of Padua – Department of Information Engineering\
**Academic Year:** 2024–2025\
**Authors:**
* Umberto Bianchin – [umberto.bianchin@studenti.unipd.it] 
* Francesco De Nicola – [francesco.denicola@studenti.unipd.it]

---

## Project Overview

This project presents a complete experimental framework for solving the **Traveling Salesman Problem (TSP)** using multiple algorithmic families:

* **Heuristics:** Nearest Neighbor, Extra-Mileage
* **Metaheuristics:** GRASP, VNS, Tabu Search, Genetic Algorithm
* **Exact Methods:** Benders Decomposition, Branch-and-Cut with CPLEX and Concorde
* **Matheuristics:** Hard Fixing, Local Branching

Each method is implemented in C and shares unified data structures, solution trackers, and input parsers. We also provide a **performance profiler** that benchmarks all algorithms developed by prof.
Domenico Salvagnin from the University of Padua called perfprof.py.

## Requirements

- **C Compiler**: GCC or MSVC
- **Python 3** with `matplotlib` and `numpy`
- **IBM ILOG CPLEX** 
- **Concorde TSP Solver** for enhanced exact solving

--- 

## Command Line Options

To run the program, you can pass various options like this:

```bash
./main -file <filename> -a <algorithm> -t <time_limit> ...
```

The following options are supported (aliases are accepted as shown):

### Input

* `-file`, `-input`, `-f <filename>`
  Specifies the `.tsp` input file containing the TSP instance (TSPLIB format).

* `-n`, `-nodes <int>`
  Generates a random instance with the given number of nodes.

### Algorithm Selection

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

### General Configuration

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

### Algorithm-Specific Parameters

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

### Help

* `-help`, `--help`
  Displays the help message and usage instructions.

---

## Performance Profiling

Use the C implementation to populate CSV logs and `perfprof.py` (Python) to visualize:

```bash
python3 perfprof.py results.csv output.pdf
```

---

## Installation – IBM ILOG CPLEX

See full installation instructions in the report (Appendix A), or:

* **Windows:** Add CPLEX and CONCERT paths to the `Path` system variable.
* **macOS:** Export `CPLEX_HOME` and extend `PATH`.
* **VS Code:** Use `c_cpp_properties.json` for header inclusion.

---

## References

1. D. Applegate et al. *Concorde TSP Solver*. [http://www.math.uwaterloo.ca/tsp/concorde.html](http://www.math.uwaterloo.ca/tsp/concorde.html). Accessed: May 30, 2025.

2. D. Applegate, R. Bixby, V. Chvátal, W. Cook. *The Traveling Salesman Problem: A Computational Study*. Princeton University Press, 2006.

3. J.F. Benders. "Partitioning procedures for solving mixed-variables programming problems." *Numerische Mathematik*, 4 (1962/63), pp. 238–252. [http://eudml.org/doc/131533](http://eudml.org/doc/131533)

4. E. Carrizosa et al. "Counterfactual Optimization for Fault Prevention in Complex Wind Energy Systems." *European Journal of Operational Research*, preprint, May 2025.

5. G.A. Croes. "A method for solving traveling-salesman problems." *Operations Research*, 6(6), 1958, pp. 791–812.

6. G. Dantzig, R. Fulkerson, S. Johnson. "Solution of a large-scale traveling-salesman problem." *Journal of the Operations Research Society of America*, 2(4), 1954, pp. 393–410.

7. D. Davendra (ed.). *Traveling Salesman Problem: Theory and Applications*. IntechOpen, 2010. ISBN: 978-953-307-426-9. [https://doi.org/10.5772/38734](https://doi.org/10.5772/38734)

8. E.D. Dolan, J.J. Moré. "Benchmarking optimization software with performance profiles." *Mathematical Programming*, 91(2), 2002, pp. 201–213.

9. T. Feo, M.G.C. Resende. "Greedy Randomized Adaptive Search Procedures." *Journal of Global Optimization*, 6, 1995, pp. 109–133. [https://doi.org/10.1007/BF01096763](https://doi.org/10.1007/BF01096763)

10. M. Fischetti, A. Lodi. "Local branching." *Mathematical Programming*, 98, 2003, pp. 23–47. [https://doi.org/10.1007/s10107-003-0395-5](https://doi.org/10.1007/s10107-003-0395-5)

11. F. Glover, M. Laguna. *Tabu Search I*. Interfaces, 1(3), 1999, pp. 190–206. [https://doi.org/10.1287/ijoc.1.3.190](https://doi.org/10.1287/ijoc.1.3.190)

12. W\.R. Hamilton. "Account of the Icosian Calculus." *Proceedings of the Royal Irish Academy*, 6, 1858, pp. 415–416.

13. IBM Corporation. *IBM ILOG CPLEX Optimization Studio*. [https://www.ibm.com/products/ilog-cplex-optimization-studio](https://www.ibm.com/products/ilog-cplex-optimization-studio). Accessed: May 30, 2025.

14. S. Lin. "Computer solutions of the traveling salesman problem." *Bell System Technical Journal*, 44(10), 1965, pp. 2245–2269.

15. S. Lin, B.W. Kernighan. "An Effective Heuristic Algorithm for the Traveling-Salesman Problem." *Operations Research*, 21(2), 1973, pp. 498–516.

16. N. Mladenović, P. Hansen. "Variable neighborhood search." *Computers & Operations Research*, 24(11), 1997, pp. 1097–1100. [https://doi.org/10.1016/S0305-0548(97)00031-2](https://doi.org/10.1016/S0305-0548%2897%2900031-2)

17. J.V. Potvin. "Genetic Algorithms for the Traveling Salesman Problem." *Annals of Operations Research*, 63, 1996, pp. 337–370. [https://doi.org/10.1007/BF02125403](https://doi.org/10.1007/BF02125403)

18. G. Reinelt. "TSPLIB – A Traveling Salesman Problem Library." *INFORMS Journal on Computing*, 3(4), 1991, pp. 376–384. [https://doi.org/10.1287/ijoc.3.4.376](https://doi.org/10.1287/ijoc.3.4.376)

19. D.J. Rosenkrantz, R.E. Stearns, P.M. Lewis II. "An Analysis of Several Heuristics for the Traveling Salesman Problem." *SIAM Journal on Computing*, 6(3), 1977, pp. 563–581. [https://doi.org/10.1137/0206041](https://doi.org/10.1137/0206041)

---

## Evaluation & Results

We tested all algorithms on synthetic TSP instances of 300 to 1000 nodes with various time limits. The **best-performing** methods (in terms of cost and time) were:

* **Hard Fixing (Matheuristic)**
* **Variable Neighborhood Search (Metaheuristic)**
* **Nearest Neighbor + 2-OPT (Heuristic)**

---

## Links

[GitHub Repository](https://github.com/umberto-bianchin/OperationsResearch2)

[Download the full report (PDF)](./Operations_Research_2.pdf) for in-depth explanations, comparisons, and performance plots.
