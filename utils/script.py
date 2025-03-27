import subprocess

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
    "./tsp -n 1000 -t 60 -r b -a t -maxt 50 -mint 10 -stept 5",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a t -maxt 155 -mint 25 -stept 10",     
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a t -maxt 250 -mint 50 -stept 25",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a t -maxt 500 -mint 100 -stept 50",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a t -maxt 155 -mint 25 -stept 10",
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

for command in commands:
    print(f"Executing: {commands}")
    subprocess.run(commands, shell=True, check=True)
    print(f"Completed: {commands}")