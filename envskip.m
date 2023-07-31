n_b = 12;
hwmult = 1;
betmult = 0;
n_s = 20;
n_st = 3; 
%disp(arrmu*n_s/2*hwt_i)
state0 = zeros(n_st,n_b);
state0(n_st,:) = ones(1,n_b)*0;
state0(1,:) = ones(1,n_b)*n_s;
nb_start = n_b;
tot_mod = nb_start;
