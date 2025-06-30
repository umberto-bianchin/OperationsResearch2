Performance Profile by D. Salvagnin (2016)
Internal use only, not to be distributed   

Dolan ED, Moŕe J (2002) Benchmarking optimization software with performance profiles. Mathematical Programming 91(2):201–213

@article{dolan2002benchmarking,
  title={Benchmarking optimization software with performance profiles},
  author={Dolan, E. D. and Mor{\'e}, J.J.},
  journal={Mathematical Programming},
  volume={91},
  number={2},
  pages={201--213},
  year={2002},
  publisher={Springer}
}


Usage:
python3 utils/perfprof.py -D ',' -T 100 -S 0.001 -M 1.015 results/results_test.csv results/pp_test.pdf -P "Example" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.01 -M 1.2 results/results_heuristics.csv results/pp_heuristics.pdf -P "Heuristics Algorithm" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 1 -M 5 results/results_exact.csv results/pp_exact.pdf -P "Exact Algorithms" -X 'Time Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.1 -M 3 results/results_C.csv results/ppC.pdf -P "Branch & Cut" -X 'Time Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.1 -M 1.13 results/results_B.csv results/ppB.pdf -P "Benders" -X 'Time Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.1 -M 1.6 results/results_C_depth.csv results/ppC_depth.pdf -P "Posting Depth Tuning" -X 'Time Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.001 -M 1.01 results/results_H_0.csv results/ppH_0.pdf -P "Hard Fixing" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.01 -M 1.1 results/results_H_100.csv results/ppH_100.pdf -P "Hard Fixing" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.01 -M 1.1 results/results_H_1000.csv results/ppH_1000.pdf -P "Hard Fixing" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.01 -M 1.025 results/results_H.csv results/ppH.pdf -P "Hard Fixing" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.01 -M 1.04 results/results_L_5000.csv results/ppL_5000.pdf -P "Local Branching" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.01 -M 1.025 results/results_L_7000.csv results/ppL_7000.pdf -P "Local Branching" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.01 -M 1.025 results/results_L_9000.csv results/ppL_9000.pdf -P "Local Branching" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.01 -M 1.02 results/results_L_20000.csv results/ppL_20000.pdf -P "Local Branching" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.01 -M 1.01 results/results_L_depth.csv results/ppL_DEPTH.pdf -P "Local Branching" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.01 -M 1.04 results/results_L_best.csv results/ppL.pdf -P "Local Branching" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.1 -M 1.45 results/results_Pold.csv results/ppGA.pdf -P "Genetic Algorithm" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.01 -M 1.07 results/results_GAlow.csv results/ppGA_low.pdf -P "Genetic Algorithm" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.001 -M 1.01 results/results_math.csv results/pp_math.pdf -P "Metaheuristic Algorithms" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.001 -M 1.1 results/results_heuristic.csv results/pp_heuristic.pdf -P "Metaheuristic Algorithms" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.001 -M 1.014 results/results_T.csv results/ppT.pdf -P "Tabu Search" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.001 -M 1.014 results/results_V.csv results/ppV.pdf -P "Variable Neighborhood Search" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.001 -M 1.032 results/results_G.csv results/ppG.pdf -P "GRASP" -X 'Cost Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.1 -M 1.26 results/results.csv results/ppB.pdf -P "Benders" -X 'Time Ratio'
python3 utils/perfprof.py -D ',' -T 100 -S 0.1 -M 1.07 results/results.csv results/ppTop.pdf -P "Best Algorithms" -X 'Cost Ratio'













Parameters:

self.parser.add_option("-D", "--delimiter", dest="delimiter", default=None, help="delimiter for input files")
self.parser.add_option("-M", "--maxratio", dest="maxratio", default=4, type=int, help="maxratio for perf. profile")
self.parser.add_option("-S", "--shift", dest="shift", default=0, type=float, help="shift for data")
self.parser.add_option("-L", "--logplot", dest="logplot", action="store_true", default=False, help="log scale for x")
self.parser.add_option("-T", "--timelimit", dest="timelimit", default=1e99, type=float, help="time limit for runs")
self.parser.add_option("-P", "--plot-title", dest="plottitle", default=None, help="plot title")
self.parser.add_option("-X", "--x-label", dest="xlabel", default='Time Ratio', help="x axis label")
self.parser.add_option("-B", "--bw", dest="bw", action="store_true", default=False, help="plot B/W")


