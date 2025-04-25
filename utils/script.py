import subprocess
from datetime import datetime
import threading

commands = [
    "./tsp -n 1000 -t 60 -r b -a v -kick 1 -kopt 3",
    "./tsp -n 1000 -t 60 -r b -a v -kick 5 -kopt 3",
    "./tsp -n 1000 -t 60 -r b -a v -kick 1 -kopt 5",
    "./tsp -n 1000 -t 60 -r b -a v -kick 3 -kopt 5",
    "./tsp -n 1000 -t 60 -r b -a v -kick 1 -kopt 7",
    "./tsp -n 1000 -t 60 -r b -a t -mint 10 -maxt 50 -stept 5",
    "./tsp -n 1000 -t 60 -r b -a t -mint 25 -maxt 155 -stept 10",     
    "./tsp -n 1000 -t 60 -r b -a t -mint 50 -maxt 250 -stept 25",
    "./tsp -n 1000 -t 60 -r b -a t -mint 100 -maxt 500 -stept 50",
    "./tsp -n 1000 -t 60 -r b -a t -mint 200 -maxt 500 -stept 50",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 20 -minc 3",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 50 -minc 3",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 100 -minc 3",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 20 -minc 7",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 50 -minc 7"    
]

best_commands = [
    "./tsp -n 1000 -t 60 -r b -a v -kick 5 -kopt 3",
    "./tsp -n 1000 -t 60 -r b -a t -mint 200 -maxt 500 -stept 50",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 20 -minc 3"  
]

cplex_commands = [
    "./tsp -r t -a b -n 200 -t 60 -warmup 0",
    "./tsp -r t -a b -n 200 -t 60 -warmup 1",
    "./tsp -r t -a c -n 200 -t 60 -warmup 0 -posting 0 -concorde 0",
    "./tsp -r t -a c -n 200 -t 60 -warmup 1 -posting 0 -concorde 0",
    #"./tsp -r t -a c -n 200 -t 60 -warmup 0 -posting 1 -concorde 0 -depth 3",
    #"./tsp -r t -a c -n 200 -t 60 -warmup 0 -posting 1 -concorde 0 -depth 5",
    "./tsp -r t -a c -n 200 -t 60 -warmup 0 -posting 1 -concorde 0 -depth 10",
    "./tsp -r t -a c -n 200 -t 60 -warmup 0 -posting 0 -concorde 1", 
    "./tsp -r t -a c -n 200 -t 60 -warmup 1 -posting 0  -concorde 1 ",
    "./tsp -r t -a c -n 200 -t 60 -warmup 1 -posting 1 -concorde 0 -depth 10",
    "./tsp -r t -a c -n 200 -t 60 -warmup 0 -posting 1 -concorde 1 -depth 10",
    "./tsp -r t -a c -n 200 -t 60 -warmup 1 -posting  1 -concorde 1 -depth 10",
]

cplex_best_commands = [
    "./tsp -r t -a b -n 200 -t 60 -warmup 1",
    "./tsp -r t -a c -n 200 -t 60 -warmup 0 -posting 0 -concorde 0"
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
