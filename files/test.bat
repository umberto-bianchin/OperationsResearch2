@echo off
echi Start FOR cicle
for /L %%i in (1, 1, 10) do (
	echo "tsp -n 1000 -time_limit 60 -seed %%i > run_n1000_t60_seed%%i.log"
)
echo End FOR cicle
pause
000