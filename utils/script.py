import subprocess
from datetime import datetime

commands = [
    "./tsp -n 1000 -t 60 -r b -a v -kick 1 -kopt 3",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a v -kick 5 -kopt 3",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a v -kick 1 -kopt 5",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a v -kick 3 -kopt 5",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a v -kick 1 -kopt 7",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a t -mint 10 -maxt 50 -stept 5",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a t -mint 25 -maxt 155 -stept 10",     
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a t -mint 50 -maxt 250 -stept 25",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a t -mint 100 -maxt 500 -stept 50",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 20 -minc 3",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 50 -minc 3",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 100 -minc 3",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 20 -minc 7",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 50 -minc 7"    
]

best_commands = [
    "./tsp -n 1000 -t 60 -r b -a v -kick 3 -kopt 5",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a t -mint 100 -maxt 500 -stept 50",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a g -alpha 20 -minc 3"  
]

print(f"Script started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

for command in commands:
    print(f"Executing: {command}")
    subprocess.run(command, shell=True, check=True)
    print(f"Completed: {command}")

print(f"Script completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
