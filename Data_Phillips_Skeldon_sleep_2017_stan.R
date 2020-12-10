# Phillips Skeldon sleep 2017 stan data

# initial values
y0 = c(0.6,-6.7,14.6,0.35,-0.86,-0.59)

## all time is in hours

# time step 
time_step = 0.05

# time values (not including time 0)
ts = seq(time_step,48, time_step)  

# number of points (not including time 0)
T_n = length(ts)

state_test = readRDS("L:/Lab_JamesR/lachlanW/sleep model/Phillips2007SleepStan/state_tester.rds")

prop_test =  sum(state_test==1)/length(state_test)