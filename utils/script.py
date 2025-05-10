import subprocess
from datetime import datetime
import threading

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
    "./tsp -r t -a b -n 300 -t 60 -warmup 1 -posting 0 -concorde 01",
    "./tsp -r t -a c -n 300 -t 60 -warmup 1 -posting 1 -depth 100 -concorde 1"
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
for command in cplex_best_commands:
    print(f"Executing: {command}")
    subprocess.run(command, shell=True, check=True)
    print(f"Completed: {command}")
'''
for thread in threads:
    thread.join()
'''
print(f"Script completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
