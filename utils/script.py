import subprocess

commands = [
    "./tsp -n 1000 -t 60 -r b -a v -kick 5",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a v -kick 10",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a t -maxt 500 -mint 100 -stept 50",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a t -maxt 50 -mint 10 -stept 5",
    "sleep 2",
    "./tsp -n 1000 -t 60 -r b -a t -maxt 5000 -mint 1000 -stept 100"
]

for command in commands:
    print(f"Executing: {commands}")
    subprocess.run(commands, shell=True, check=True)
    print(f"Completed: {commands}")