n_b = 4;
hwmult = 1;
betmult = 1;
n_s = 8;
n_st = 3; 
%disp(arrmu*n_s/2*hwt_i)
state0 = zeros(n_st,n_b);
state0(n_st,:) = ones(1,n_b)*0;
state0(1,:) = ones(1,n_b);
