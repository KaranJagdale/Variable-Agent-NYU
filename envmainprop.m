
n_s =8; %Number of Stations
n_b = 2; %Number of bus, 1 bus = 2 modules
n_st = 3; %Number of states 
n_a = 4; % Number of actions bus can take
%a_par = rand(1,n_s); %These are Ps values for the stops
almu = 2/n_s;
%a_par = normrnd(almu, almu/10, 1, n_s); 
%arr_par = rand(1,n_s)/60*2; % Assuming on 90 passengers arrive in 30 minutes
arrmu = 0.015/3; arrsigma = arrmu/10;
%arr_par = normrnd(arrmu, arrsigma, 1,n_s);
dis_stp = 300*(rand(1,n_s) + 1); %Distance between stops distributed between 300 - 600 meters
v_bus = 20*5/18; % Speed of bus 20 Km/h
w_wait = 2.1; w_walk = 2.2;% weights of walk time and wait time in final cost
cap_bus = 50;
unit_cap = cap_bus/2;
v_pas = 5.4*5/18; %Passenger speed in Km/h
t_bo = 5; %boarding time per passenger in seconds
t_al = 2;
fixdt = 30; %fixed time lost per stop
splittime = 1;

