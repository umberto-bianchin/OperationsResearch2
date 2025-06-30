import subprocess
from datetime import datetime
import threading
import os

tabu_commands = [
    "./tsp -n 1000 -t 60 -r b -a t -mint 300 -maxt 700 -stept 50",
    "./tsp -n 1000 -t 60 -r b -a t -mint 300 -maxt 700 -stept 100",
    "./tsp -n 1000 -t 60 -r b -a t -mint 500 -maxt 700 -stept 50",
    "./tsp -n 1000 -t 60 -r b -a t -mint 500 -maxt 700 -stept 100",
    "./tsp -n 1000 -t 60 -r b -a t -mint 700 -maxt 900 -stept 50",
    "./tsp -n 1000 -t 60 -r b -a t -mint 700 -maxt 900 -stept 100",
]

grasp_commands = [
    "./tsp -n 1000 -t 60 -r b -a g -alpha 10 -minc 7",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 5 -minc 3",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 5 -minc 7",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 3 -minc 3",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 3 -minc 7",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 2 -minc 3",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 2 -minc 7",
]

variable_commands = [
    "./tsp -n 1000 -t 60 -r b -a v -kick 5 -kopt 3",
    "./tsp -n 1000 -t 60 -r b -a v -kick 3 -kopt 5",
    "./tsp -n 1000 -t 60 -r b -a v -kick 5 -kopt 5",
    "./tsp -n 1000 -t 60 -r b -a v -kick 7 -kopt 3",
    "./tsp -n 1000 -t 60 -r b -a v -kick 7 -kopt 5",
]

best_heu_commands = [
    "./tsp -n 1000 -t 60 -r b -a v -kick 7 -kopt 5",
    "./tsp -n 1000 -t 60 -r b -a t -mint 700 -maxt 900 -stept 50",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 2 -minc 3"  
]

branch_commands = [
    "./tsp -r t -a c -n 300 -t 60 -warmup 0 -posting 0 -concorde 0",
    "./tsp -r t -a c -n 300 -t 60 -warmup 1 -posting 0 -concorde 0",
    "./tsp -r t -a c -n 300 -t 60 -warmup 0 -posting 1 -depth 100 -concorde 0",
    "./tsp -r t -a c -n 300 -t 60 -warmup 0 -posting 0 -concorde 1", 
    "./tsp -r t -a c -n 300 -t 60 -warmup 1 -posting 0 -concorde 1 ",
    "./tsp -r t -a c -n 300 -t 60 -warmup 1 -posting 1 -depth 100 -concorde 0",
    "./tsp -r t -a c -n 300 -t 60 -warmup 0 -posting 1 -depth 100 -concorde 1",
    "./tsp -r t -a c -n 300 -t 60 -warmup 1 -posting 1 -depth 100 -concorde 1",
]

benders_command = [
    "./tsp -r t -a b -n 300 -t 60 -warmup 0 -posting 0 -concorde 0",
    "./tsp -r t -a b -n 300 -t 60 -warmup 1 -posting 0 -concorde 0",
]

posting_commands = [
    "./tsp -r t -a c -n 300 -t 60 -warmup 0 -posting 1 -depth 10 -concorde 0",
    "./tsp -r t -a c -n 300 -t 60 -warmup 0 -posting 1 -depth 50 -concorde 0",
    "./tsp -r t -a c -n 300 -t 60 -warmup 0 -posting 1 -depth 100 -concorde 0",
    "./tsp -r t -a c -n 300 -t 60 -warmup 0 -posting 1 -depth 1000 -concorde 0"
]

cplex_best_commands = [
    "./tsp -r t -a b -n 300 -t 60 -warmup 1 -posting 0 -concorde 0",
    "./tsp -r t -a c -n 300 -t 60 -warmup 1 -posting 1 -depth 100 -concorde 1"
]

hard_fixing_commands_0 = [
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 0 -cdepth 0",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 1 -probability 20 -cdepth 0",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 1 -probability 40 -cdepth 0",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 1 -probability 50 -cdepth 0",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 1 -probability 80 -cdepth 0",
]

hard_fixing_commands_100 = [
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 0 -cdepth 100",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 1 -probability 20 -cdepth 100",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 1 -probability 40 -cdepth 100",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 1 -probability 50 -cdepth 100",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 1 -probability 80 -cdepth 100",
]

hard_fixing_commands_1000 = [
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 0 -cdepth 1000",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 1 -probability 20 -cdepth 1000",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 1 -probability 40 -cdepth 1000",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 1 -probability 50 -cdepth 1000",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 1 -probability 80 -cdepth 1000",
]

hard_fixing_commands = [
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 0 -cdepth 1000",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 0 -cdepth 2000",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 0 -cdepth 3000",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 0 -cdepth 5000",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 0 -cdepth 7000",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 0 -cdepth 9000",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 0 -cdepth 20000",
]

local_branching_commands_5000 = [
    #"./tsp -r b -a l -n 1000 -t 180 -klocal 10 -cdepth 5000",
    #"./tsp -r b -a l -n 1000 -t 180 -klocal 20 -cdepth 5000",
    #"./tsp -r b -a l -n 1000 -t 180 -klocal 30 -cdepth 5000",
    "./tsp -r b -a l -n 1000 -t 180 -klocal 50 -cdepth 5000",
    "./tsp -r b -a l -n 1000 -t 180 -klocal 70 -cdepth 5000",
]

local_branching_commands_7000 = [
    #"./tsp -r b -a l -n 1000 -t 180 -klocal 10 -cdepth 7000",
    #"./tsp -r b -a l -n 1000 -t 180 -klocal 20 -cdepth 7000",
    #"./tsp -r b -a l -n 1000 -t 180 -klocal 30 -cdepth 7000",
    "./tsp -r b -a l -n 1000 -t 180 -klocal 50 -cdepth 7000",
    "./tsp -r b -a l -n 1000 -t 180 -klocal 70 -cdepth 7000",
]

local_branching_commands_9000 = [
    "./tsp -r b -a l -n 1000 -t 180 -klocal 10 -cdepth 9000",
    "./tsp -r b -a l -n 1000 -t 180 -klocal 20 -cdepth 9000",
    "./tsp -r b -a l -n 1000 -t 180 -klocal 30 -cdepth 9000",
]

local_branching_commands = [
    "./tsp -r b -a l -n 1000 -t 180 -klocal 10 -cdepth 9000",
    "./tsp -r b -a l -n 1000 -t 180 -klocal 20 -cdepth 9000",
    "./tsp -r b -a l -n 1000 -t 180 -klocal 30 -cdepth 9000",
    "./tsp -r b -a l -n 1000 -t 180 -klocal 50 -cdepth 9000",
    "./tsp -r b -a l -n 1000 -t 180 -klocal 70 -cdepth 9000",
]

local_branching_commands_20000 = [
    #"./tsp -r b -a l -n 1000 -t 180 -klocal 10 -cdepth 20000",
    #"./tsp -r b -a l -n 1000 -t 180 -klocal 20 -cdepth 20000",
    #"./tsp -r b -a l -n 1000 -t 180 -klocal 30 -cdepth 20000",
    "./tsp -r b -a l -n 1000 -t 180 -klocal 50 -cdepth 20000",
    "./tsp -r b -a l -n 1000 -t 180 -klocal 70 -cdepth 20000",
]

genetic_commands = [
    #"./tsp -r b -a p -n 1000 -t 60 -population 500 -generation 20",
    #"./tsp -r b -a p -n 1000 -t 60 -population 500 -generation 50",
    #"./tsp -r b -a p -n 1000 -t 60 -population 1000 -generation 20",
    #"./tsp -r b -a p -n 1000 -t 60 -population 1000 -generation 50",
    #"./tsp -r b -a p -n 1000 -t 60 -population 1000 -generation 100",
    "./tsp -r b -a p -n 1000 -t 60 -population 200 -generation 20",
    "./tsp -r b -a p -n 1000 -t 60 -population 200 -generation 50",
    "./tsp -r b -a p -n 1000 -t 60 -population 50 -generation 20",
    "./tsp -r b -a p -n 1000 -t 60 -population 25 -generation 10",
    "./tsp -r b -a p -n 1000 -t 60 -population 10 -generation 5",
]

best_commands = [
    #"./tsp -r b -a n -n 1000 -t 180",
    #"./tsp -r b -a v -n 1000 -t 180 -kick 7 -kopt 5",
    "./tsp -r b -a h -n 1000 -t 180 -fixedprob 0 -cdepth 7000",
]

def run_command(command):
    try:
        print(f"Executing: {command}")
        subprocess.run(command, shell=True, check=True)
        print(f"Completed: {command}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while executing {command}: {e}")


print(f"Script started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

'''
threads = []
for command in best_commands:
    thread = threading.Thread(target=run_command, args=(command,))
    threads.append(thread)
    thread.start()
'''
for command in best_commands:
    print(f"Executing: {command}")
    subprocess.run(command, shell=True, check=True)
    print(f"Completed: {command}")

'''
src = "./results/results_L.csv"
dst = "./results/results_L_5.csv"
if os.path.exists(src):
    os.rename(src, dst)
    print(f"Rinominato {src} → {dst}")
else:
    print(f"Attenzione: file '{src}' non trovato, impossibile rinominare.")

for command in local_branching_commands_7000:
    print(f"Executing: {command}")
    subprocess.run(command, shell=True, check=True)
    print(f"Completed: {command}")

src = "./results/results_L.csv"
dst = "./results/results_L_7.csv"
if os.path.exists(src):
    os.rename(src, dst)
    print(f"Rinominato {src} → {dst}")
else:
    print(f"Attenzione: file '{src}' non trovato, impossibile rinominare.")

for command in local_branching_commands_20000:
    print(f"Executing: {command}")
    subprocess.run(command, shell=True, check=True)
    print(f"Completed: {command}")

src = "./results/results_L.csv"
dst = "./results/results_L_20.csv"
if os.path.exists(src):
    os.rename(src, dst)
    print(f"Rinominato {src} → {dst}")
else:
    print(f"Attenzione: file '{src}' non trovato, impossibile rinominare.")
'''
'''
for thread in threads:
    thread.join()
'''
print(f"Script completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

233863.6375
232995.8414
230850.5722
227828.8611
232587.7922
232676.5585
228497.4295
232336.5610
233097.7538
231739.0618
229638.8822
231336.9911
228190.9787
229701.2904
229048.7810
234813.9059
232587.4142
232604.5057
232934.6679
231680.9228
229361.2368
233810.0243
231238.3399
232403.7551
234554.4706
234978.5645
231079.2201
230512.3593
233027.6034
230204.4778
230349.6246
235459.5160
239681.9061
232674.7494
232114.3898
231919.3137
230573.5891
232953.4879
233783.4078
231515.1579
231588.3945
230453.1125
234466.6246
231206.1087
230700.1572
235189.5948
234002.1790
228170.9100
232163.6320