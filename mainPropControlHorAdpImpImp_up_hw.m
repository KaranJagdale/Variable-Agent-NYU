%%Code to simulate the variable agent algorithm
rng(21) %95
close all
%CHECK Reward_t function before running this code specifically check the
%gain values fPas1, fPas2
%This code is using the most recent headway function that is headwaycuu
%this code also accounts the leftover passengers if the bus skips the stop
%but at this point of time we are assming that the bus not skipping stops
%successively
%dbstop if naninf
nom = true;%change the boolian value to obtain performance of nominal or proposed policy

tic;
simtime = 3; %simulation time in hours
n_s = 20; %Number of Stations
n_b = 24; %Number of bus, 1 bus = 2 modules modules starts in the detached state

nb_start = n_b;
tot_mod = nb_start;
n_st = 3; %Number of states 
n_a = 5; % Number of actions bus can take
%a_par = rand(1,n_s); %These are Ps values for the stops
hor = 1;
fPas1 =0; fPas2 =1; %nominal 0, 1.5
jP1 = 1; jP2 = 1; %reward function join parameters
almu = 2/n_s;
a_par = normrnd(almu, almu/10, 1, n_s); 
%arr_par = rand(1,n_s)/60*2; % Assuming on 90 passengers arrive in 30 minutes
hdem = 1500;
arrmu = hdem/3600/n_s; arrsigma = arrmu/10;
arr_par = normrnd(arrmu, arrsigma, 1,n_s);
atstop = zeros(1,n_b); %These are the flags which will be 1  if hte particular bus is waiting to leave a stop
% if the corresponding bus is at stop and 0 it it is on the road
fixdt = 20; %fixed time lost per stop
%dis_stp = 333*(rand(1,n_s) + 1); %Distance between stops distributed between 333 - 666 meters
dis_stp = normrnd(400, 40, 1,n_s);
v_bus = 20*5/18; % Speed of bus 20 Km/h
w_wait = 2.1; w_walk = 2.2;% weights of walk time and wait time in final cost
lapass = zeros(1,n_b);
lpass = zeros(1,n_s);
w_pass = zeros(n_b,n_s); %to store the number of walking passengers to a particular stop if the modules are leaving in the joined state, first parameter should be 2*n_b
%at max there can be 4 types of passengers walking towards a particular
%stop      
Hway = zeros(1,n_b);
tw_pass = zeros(2*n_b,n_s); %time to reach the next stop
cap_bus = 80;
unit_cap = cap_bus/2;
v_pas = 4.5*5/18; %Passenger speed in Km/h
t_bo = 4; %boarding time per passenger in seconds
t_al = 3; %Alighting time per passenger in seconds
scamcount = 0;
gam_k = 2; gam_theta = 5.1; %parameters of Gamma rv
gamma = 1.5;
bet = 0.6;
hwt_i = (sum(dis_stp)/v_bus + n_s*fixdt)/(n_b - (t_al + t_bo)*arrmu*n_s);
hwt = hwt_i;
hwA = hwt*ones(1, n_b); %array to store the headways of the buses
cumtime = zeros(1,n_b);
splittime = 1; %splittime being 1 second -- Introducing this constant to avoid bus overtake
stload = floor(arrmu*n_s/2*hwt_i); 
state0 = zeros(n_st,n_b);
%state0(n_st,:) = ones(1,n_b); %uncomment if the modules start in the
%attached state
state0(1,:) = ones(1,n_b)*n_s;
state0(2,:) = stload;
state = state0;

t_nxt = zeros(1,n_b); %this variable stores the time bus will require to reach the upcoming stop
for i=1:n_b  %At start each bus leaves with an interval of 90 seconds
    t_nxt(i) = hwt_i*(i-1);% + dis_stp(1)/v_bus; as the bus are starting 
end
%disp(t_nxt)
l_action = zeros(n_st,n_b); %last action taken by th bus
papbcum = 0;
papbskip = 0;
count = 1; %This will account the number of rounds
%consec_t = dis_stp/v_bus;  %We are assuming that, in the first round the passengers begin to arrive at the stop after the bus leaves the previous stop 
gencount = 1; %delierately initialized gencount to 1 to syncit with Time array
skcount = 0;
stcount = 0;
sjcount = 0;
spcount = 0;
snbcount = 0;
time = 0;
Time = 0;
Time1 = 0;
Time2 = 0;
Time3 = 0;
Time4 = 0;
State1 = state(:,1);
State2 = state(:,1);
State3 = state(:,2);
loc1 = state(1,1); loc2 = state(1,1); loc3 = state(1,2); loc4 = state(1,2);
locf1 = loc1;
reachtime = zeros(1,n_s);
leavetime = zeros(1,n_s);
nbtime = zeros(1,n_s); %array to capture the time at which a bus leaves from the s^th stop
pa = stload*n_b;
Pa = stload*n_b;
Hw = 0;
pe_cum = 0; %exiting passengers
Pe_cum = 0;
Pd_cum = 0;
Pa_cum = stload*n_b;
Pb_cum = stload*n_b;
Pw_cum = 0;
pw_cum = 0;
pd_cum = 0; %cumulative deboarding passengers
pb_cum = stload*n_b; %cumulative boarding passengers
pa_cum = stload*n_b; %cumulative arriving passengers
ord = 1:n_b; %will keep track of the order of the buses. This is basically 
next_stop = false;
rotCounted = true;
rotCB= false;
%the mapping between the state index and the bus number
%Creating arriving passenger array 
% pa_pre = poissrnd(arrmu*n_s, 1, simtime*3600);
Pa_pre_cum = zeros(1,simtime*3600+1);
% pa_pre_cum = stload*n_b;
Time_pre = 0:simtime*3600;

%generating arrivals for each stop
tflst = sum(dis_stp)/v_bus + (t_al + t_bo)*hwt_i*sum(arr_par) + fixdt*n_s;
tfstops = zeros(1,n_s);
for i = 1:n_s -1
    tfstops(i+1) = tfstops(i) + dis_stp(i)/v_bus + (t_al + t_bo)*hwt_i*arr_par(i) + fixdt;
end
Pa_all_cum = zeros(n_s, simtime*3600+1);
nwpasscount = 0;
streachct = 0;
skgstct = 0;
rewc = 0;
for i = 1:simtime*3600
    for st = 1:n_s
        if i < tflst
            if tfstops(st) < hwt_i
                newpass = poissrnd(arr_par(st));
            else
                if tfstops(st) - i <= hwt_i
                    newpass = poissrnd(arr_par(st));
                else
                    newpass = 0;
                end
            end
        else
            newpass = poissrnd(arr_par(st));
        end
        Pa_all_cum(st, i+1) = Pa_all_cum(st, i) + newpass;
    end
end
padec = zeros(n_s, simtime*3600+1);
for i = 1: 3600*simtime + 1
    Pa_pre_cum(i) = sum(Pa_all_cum(:,i)) + stload*n_b;
end

rotcount = 0;
ns_reach = 0; %# of times the bus reached the state n_s
stcn = 0; spcn = 0; skcn = 0; jncn = 0; nbcn = 0; %action counters in particular time interval
pVeh = 0; %number of passengers in vehicle
pWait = []; %number of wiating passengers
locs = cell(1,tot_mod); %to store the module location 
for i=1:tot_mod
    locs{i} = [];
end

mod_time = cell(1, tot_mod); %to store module time
for i=1:tot_mod
    mod_time{i} = [];
end
stopTimeA = 0; %average of stop times
while true
    if floor(pe_cum) ~= pe_cum
        disp('pe_cum is not more int')
        break
    end
    pVeh = [pVeh sum(state(2,:))];
    %calculating the value of target headway
    hwt_i = (sum(dis_stp)/v_bus + n_s*fixdt)/(n_b - (t_al + t_bo)*arrmu*n_s);
    hwt = hwt_i;

    if state(1,1) == n_s %count encodes whether the leading module has covered a circle or not
        ns_reach = ns_reach + 1;
    end
    if ns_reach > 1 %this condition means that the first bus has completed a circle
        count = 0;
    end
%     if count == 0
%         fprintf('count 0 at gencount : %i \n', gencount)
%         break
%     end
    gencount = gencount + 1;
%     if gencount < 4
%         disp(state)
%     end
    %disp(state)
    n_b = size(state,2);
    hw = headwaycuuu(state, t_nxt, atstop, dis_stp, v_bus, n_b,count,ord);
    %disp(hw)


    [M,im] = min(t_nxt);  %This step finds which bus reaches/leaves the next stop first
    time = time + M;
    Time = [Time time];
    if state(1, 1) == n_s-1 && next_stop && im == 1 && atstop(1) == 0
        rotcount = rotcount + 1;
        rotCounted = false;
    end
    if rotcount == 2 && ~rotCounted
        sind = gencount;
        stime = time;
        etime = stime + 1*60*60;
        rotCounted = true;
        rotCB = true;
        inSt = locs{1}(end);
    end
    if rotCB
        if time > stime + 500 && locs{1}(end) == inSt
            cycT = (time - stime)/60;
            rotCB = false;
        end
    end
    %finding other indices with the minimum value
    imin = find(t_nxt == M);
    im = min(imin);
    t_nxt_sim = t_nxt; % This is to pass in the rerwsum function
    implus = iplus(im,n_b);
    implusp = iplus(implus,n_b);
    for i=1:n_b
            if i ~= im
                t_nxt(i) = t_nxt(i) - t_nxt(im);  %nxt for im will be calculated as the stop time at the stop
            end       
    end

    for i= 1:nb_start
        for j = 1:n_s
            if w_pass(i,j) ~=0
                tw_pass(i,j) = tw_pass(i,j) - t_nxt(im);
                if tw_pass(i,j) < 0 || tw_pass(i,j) == 0
                    pe_cum = pe_cum + w_pass(i,j);
                    w_pass(i,j) = 0;
                    tw_pass(i,j) = 0;
                end
            end
        end
    end

    if atstop(im) == 0
        
       % hway = hw(im);
        if im ==1 && count == 1
            hway = hwt;
        else
            hway = hwA(im);
        end
        %hway = 
        
        if state(1,im) == n_s %Because our stops are circular n_s->1
            state(1,im) = 1;
        else
            state(1,im) = state(1,im) + 1; %state(1,im) is storing the last bus-stop number the bus has visited
        end
        opt_act = rewsum_t(hor, state, t_nxt_sim, atstop, t_al, t_bo, arr_par, a_par, dis_stp, im, hway, v_bus, v_pas, lpass, lapass, cap_bus, hwt,ord, count, gencount, reachtime, leavetime, time,Hway, Pa_all_cum, fPas1, fPas2, jP1, jP2, ex_wt);
        %opt_act = rewsum_expt(hor, state, t_nxt_sim, atstop, t_al, t_bo, arr_par, a_par, dis_stp, im, hway, v_bus, v_pas, lpass, lapass, cap_bus, hwt,ord, count, gencount);

        optSeq = opt_act;
        bstop = state(1,im);
        opt_act = str2num(opt_act(1));
        
        imr = opt_act + 1; %just a shift as the code below is written in terms of 
        if nom
            imr = 1;
        end
        b_st = [];
        ex_wt = 0; %extra wait the bus has to do because of already reched buses
        %following lines find the max waiting time among the buses already
        %at the stop where the vurrent bus (im) is reached.
        for i=1:n_b
            if state(1,i) == bstop && atstop(i) == 1
                b_st = [b_st i];
                if t_nxt(i) > ex_wt
                    ex_wt = t_nxt(i);
                end
            end
        end
%         if ex_wt > 0 && imr == 2
%             disp('heya')
%             disp(state)
%             fprintf('im : %i \n', im)
%         end
%         if ex_wt > 0 && state(3,im) == 0
%             imr = 1;
%         end
        
        if state(3,im) == 1 
            if imr == 3
                action  = 'split';
                spcount = spcount +1;
                n_b = n_b +1; %The number of modules increase due to splitting
                %compute the number of passengers getting down at the current stop
                pds_c = binornd((state(2,im)-lapass(im)),a_par(state(1,im))) + lapass(im);
    
                stmin = iminus(state(1,im),n_s);
                for iw=1:nb_start
                    if w_pass(iw,stmin) == 0
                        break
                    end
                end
                
                if prod(w_pass(:,stmin)) ~=0
                    disp('damn... this is not working')
                    break
                end
                w_pass(iw,stmin) = lapass(im); %If the deboarding passengers are more than unit capacity then the left over passengers gets priority for deboarding
                if w_pass(iw,stmin) ~=0    
                    tw_pass(iw,stmin) = dis_stp(stmin)/v_pas;
                end
                
                load_r = pds_c;  %floor is not necessary but having it will not affect as well
                lapassf = 0;
                if state(2,im) - load_r > unit_cap
                    load_r = state(2,im) - unit_cap;
                end
                
                if load_r > unit_cap
                    lapassf = load_r - unit_cap;
                    load_r = unit_cap;
                    pds_c = load_r;    %as all passengers cannot deboard in this case               
                end
                pdsr = pds_c;
                pd_cum = pd_cum + pdsr;
                pe_cum = pe_cum + pdsr - lapass(im);
    %             if gencount == 167
    %                 disp(load_r)
    %             end
                load_f = state(2,im) - load_r;
    
                state_i = zeros(n_st,n_b);
                for i=1:im-1
                        state_i(:,i) = state(:,i);
                end
                state_i(:,im) = [state(1,im);load_f;0];
                state_i(:,im+1) = [state(1,im);load_r;0];
                for i=im+2:n_b
                        state_i(:,i) = state(:,i-1);
                end
         
                state = state_i;
    
                %Updating atstop
                atstop_i = zeros(1,n_b);
                for i=1:im-1
                    atstop_i(i) = atstop(i);
                end
                atstop_i(im) = 0;
                atstop_i(im + 1) = 1; %as rear module has stopped at the stop
                for i = im + 2 : n_b
                    atstop_i(i) = atstop(i-1);
                end
                atstop = atstop_i;
    
                %Updating lapass
    
                lapass_i = zeros(1,n_b);
                for i=1:im-1
                    lapass_i(i) = lapass(i);
                end
                lapass_i(im) = lapassf;
                lapass_i(im+1) = 0;
                for i = im+2:n_b
                    lapass_i(i) = lapass(i-1);
                end
                lapass = lapass_i; 
                
                state(2,im+1) = state(2,im+1) - pdsr; %This will be zero because the way we have distributed the passengers           
                f_ind = floor(reachtime(state(1,im))); s_ind = floor(time);
                if s_ind >= size(Pa_all_cum,2)
                    s_ind = size(Pa_all_cum,2)-1;
                end
                padec(state(1,im), f_ind+1:s_ind+1) = 1;
                pa = Pa_all_cum(state(1,im),s_ind+1) - Pa_all_cum(state(1,im),f_ind+1);
                pa_cum = pa_cum + pa;
                
                pbsr = min(pa + lpass(state(1,im)), unit_cap - state(2,im+1));
                pb_cum = pb_cum + pbsr;
                lpass(state(1,im)) = pa + lpass(state(1,im)) - pbsr; %assuming the leftover passengers gets priority boarding and they are never more than unitcap
                papbcum = papbcum + pa-pbsr;
    
                state(2,im+1) = state(2,im+1) + pbsr;
                stopTime = max(t_al*pdsr , t_bo*pbsr) + fixdt + ex_wt + splittime;
                tspent = stopTime;
                stopTimeA = stopTimeA + stopTime;
                %Updating t_nxt as the number of buses has changed
                t_nxt_i = zeros(1,n_b);
                for i = 1:im-1
                    t_nxt_i(i) = t_nxt(i);
                end
                for i = im+2:n_b
                    t_nxt_i(i) = t_nxt(i-1);
                end

                %updating hwA
                hwA_i = zeros(1,n_b);
                for i = 1:im-1
                    hwA_i(i) = hwA(i);
                end
                for i = im+2:n_b
                    hwA_i(i) = hwA(i-1);
                end
                hwA = hwA_i;


                t_nxt_i(im) = ex_wt + 1; %adding 1 to ensure the module behind leaves later than the module in the front
                atstop(im) = 1;
                %t_nxt_i(im) = ct +  ex_wt + splittime; %As the leading module does not stop.
                t_nxt_i(im + 1) = tspent; %time to leave the stop
                
                
                t_nxt = t_nxt_i;
                
            elseif imr == 1
                action = 'stop';
                stcount = stcount + 1;
                
                pdsr = lapass(im) + binornd((state(2,im) - lapass(im)),a_par(state(1,im)));
                
                pd_cum = pd_cum + pdsr;
                pe_cum  = pe_cum + pdsr - lapass(im);
                stmin = iminus(state(1,im),n_s);
                for iw=1:nb_start
                    if w_pass(iw,stmin) == 0
                        break
                    end
                end
                if prod(w_pass(:,stmin)) ~=0
                    disp('damn... this is not working')
                    break
                end
                w_pass(iw,stmin) = lapass(im); %If the deboarding passengers is more than unit capacity then the left over passengers gets priority for deboarding
                if w_pass(iw,stmin) ~=0
                    tw_pass(iw,stmin) = dis_stp(stmin)/v_pas;
                end
%
                lapass(im) = 0;
                
                state(2,im) = state(2,im) - pdsr;

                f_ind= floor(reachtime(state(1,im))); s_ind = floor(time);
                if s_ind >= size(Pa_all_cum,2)
                    s_ind = size(Pa_all_cum,2)-1;
                end
                padec(state(1,im), f_ind+1:s_ind+1) = 1;
                pa = Pa_all_cum(state(1,im),s_ind+1) - Pa_all_cum(state(1,im),f_ind+1);
                pa_cum = pa_cum + pa;

                pbsr = min(pa + lpass(state(1,im)), cap_bus - state(2,im));
                pb_cum = pb_cum + pbsr;
                lpass(state(1,im)) = pa + lpass(state(1,im)) - pbsr;
                papbcum = papbcum + pa-pbsr;
                state(2,im) = state(2,im) + pbsr;
                stopTime = max(t_al*pdsr , t_bo*pbsr) + fixdt + ex_wt;
                
                stopTimeA = stopTimeA + stopTime;
                tspent = stopTime;
                %Updating t_nxt, atstop, 
                atstop(im) = 1;
                t_nxt(im) = tspent;
            else
                pw = binornd((state(2,im) - lapass(im)),a_par(state(1,im))); %Number of passengers have to walk. Can modify this line to incorporate multiple skip leftover passengers
                pw_cum = pw_cum + pw;
                lapass(im) = pw + lapass(im);
                
                f_ind= floor(reachtime(state(1,im))); s_ind = floor(time);
                if s_ind>= size(Pa_all_cum,2)
                    s_ind = size(Pa_all_cum,2)-1;
                end
                padec(state(1,im), f_ind+1:s_ind+1) = 1;
                pa = Pa_all_cum(state(1,im), s_ind+1) - Pa_all_cum(state(1,im),f_ind+1);
                pa_cum = pa_cum + pa;
                
                lpass(state(1,im)) = pa + lpass(state(1,im));
                papbskip = papbskip + pa;
                skcount = skcount + 1;
                action = 'skip';
                tspent = 1 + ex_wt;           
                atstop(im) = 1;
                t_nxt(im) = tspent; %adding 1 to ensure the module reaching late leave late
                
            end
        

    
        else  %if the bus is splitted already then the competition will be in stop and skip
            if imr == 1
                action = 'stop';
                stcount = stcount + 1;
                pdsr = lapass(im) + binornd((state(2,im) - lapass(im)),a_par(state(1,im)));

                pd_cum = pd_cum + pdsr;
                pe_cum  = pe_cum + pdsr - lapass(im);
                stmin = iminus(state(1,im),n_s);
                for iw=1:nb_start
                    if w_pass(iw,stmin) == 0
                        break
                    end
                end
                if prod(w_pass(:,stmin)) ~=0
                    disp('damn... this is not working')
                    break
                end
                w_pass(iw,stmin) = lapass(im); %If the deboarding passengers is more than unit capacity then the left over passengers gets priority for deboarding
                if w_pass(iw,stmin) ~=0
                    tw_pass(iw,stmin) = dis_stp(stmin)/v_pas;
                end
%
                lapass(im) = 0;
                
                state(2,im) = state(2,im) - pdsr;

                f_ind= floor(reachtime(state(1,im))); s_ind = floor(time);
                if s_ind >= size(Pa_all_cum,2)
                    s_ind = size(Pa_all_cum,2)-1;
                end
                padec(state(1,im), f_ind+1:s_ind+1) = 1;
                pa = Pa_all_cum(state(1,im),s_ind+1) - Pa_all_cum(state(1,im),f_ind+1);
                pa_cum = pa_cum + pa;
                pbsr = min(pa + lpass(state(1,im)), unit_cap - state(2,im));
                pb_cum = pb_cum + pbsr;
                lpass(state(1,im)) = pa + lpass(state(1,im)) - pbsr;
                papbcum = papbcum + pa-pbsr;
                state(2,im) = state(2,im) + pbsr;
                stopTime =  max(t_al*pdsr , t_bo*pbsr) + fixdt + ex_wt;
                stopTimeA = stopTimeA + stopTime;
%                 if gencount == 534
%                     fprintf('tspent ; %f, pdsr : %f, pbsr : %f, sind : %i, find : %i \n', stopTime, pdsr, pbsr, s_ind, f_ind)
%                 end
                tspent = stopTime;
                %Updating t_nxt, atstop, 
                atstop(im) = 1;
                t_nxt(im) = tspent;

                 
            else  %If skip action is taken then the time spent at the stop will be zero
                pw = binornd((state(2,im) - lapass(im)),a_par(state(1,im))); %Number of passengers have to walk. Can modify this line to incorporate multiple skip leftover passengers
                pw_cum = pw_cum + pw;
                lapass(im) = pw + lapass(im);
                
                f_ind= floor(reachtime(state(1,im))); s_ind = floor(time);
                if s_ind>= size(Pa_all_cum,2)
                    s_ind = size(Pa_all_cum,2)-1;
                end
                padec(state(1,im), f_ind+1:s_ind+1) = 1;
                pa = Pa_all_cum(state(1,im), s_ind+1) - Pa_all_cum(state(1,im),f_ind+1);
                pa_cum = pa_cum + pa;
                
                lpass(state(1,im)) = pa + lpass(state(1,im));
                papbskip = papbskip + pa;
                skcount = skcount + 1;
                action = 'skip';
                
                tspent = 1 + ex_wt;            
                atstop(im) = 1;
                t_nxt(im) = tspent; %adding 1 to ensure the module reaching late leaves late 
            end           
        end 
        reachtime(bstop) = time;
     
    
    else %If the bus is not at the stop we will evaluate join possibility
        
        if im == 1 && count == 1
            hway = hwt;
        else
            hway = time - leavetime(state(1,im));
        end
        
        opt_act = rewsum_t(hor, state, t_nxt_sim, atstop, t_al, t_bo, arr_par, a_par, dis_stp, im, hway, v_bus, v_pas, lpass, lapass, cap_bus, hwt,ord, count, gencount, reachtime, leavetime, time,Hway, Pa_all_cum, fPas1, fPas2, jP1, jP2, ex_wt);
        %opt_act = rewsum_expt(hor, state, t_nxt_sim, atstop, t_al, t_bo, arr_par, a_par, dis_stp, im, hway, v_bus, v_pas, lpass, lapass, cap_bus,hwt,ord, count, gencount);
%         if gencount == 1783
%             fprintf('optimal sequence : %s \n', opt_act)
%         end
        optSeq = opt_act;
        %fprintf('opt_act : %s on module : %i \n \n', opt_act, im)
        
        opt_act = str2num(opt_act(1));
        imr = opt_act + 1;
        if state(1,im) == state(1,implus) && atstop(implus) == 1 && t_nxt(implus) < 31 && state(3,im) == 0 && state(3,implus) == 0
            %fprintf('t_nxt(implus) : %i, gencount : %i \n', t_nxt(implus), gencount)
            imr = 4;
        end
        if nom
            imr = 5;
        end
        
        if state(3,im) == 0 && state(3,implus) == 0 && imr == 4 %Join have positive reward and it is indeed possible to join
            r_jn = 1;
        else
            r_jn = -1; %This is dummy value if the statement in 'if' does not executte then join is anyway not feasible
        end
        %r_jn = -1;
        leavetime(state(1,im)) = time;
        
        joincond =  r_jn > 0;%  && ((atstop(implus) == 0 && state(1,implus) ~= state(1, im)) || state(1, implus) == state(1,im)); 
            if joincond
                %fprintf('join happening at time : %f, gencount %i for module : %i \n', time, gencount, im)
                %disp(table(state', atstop', t_nxt'))
                action = 'join';
                sjcount = sjcount + 1;
                if atstop(implus) == 0
                    tspent = hw(implus);  %Front bus has to wait for this time
                else
                    tspent = t_nxt(implus);
                end
                leavetime(state(1,im)) = time + tspent; %in the case of join the module will leave after additional tspent time

                load_rj = state(2,implus); % this value will requre to caculate lapass
                stop_rj = state(1,implus);
                lap_ac = 0;
                if stop_rj <= state(1,im)
                    bskip = stop_rj+1:state(1,im);
                else
                    bskip = stop_rj+1:n_s;
                    bskip = [bskip 1:state(1,im)];
                end
                for st = bskip
                    lap = floor((load_rj - lapass(implus) - lap_ac)*a_par(st));
                    lap_ac = lap_ac + lap;
                end
                
                if im == n_b                    
                    n_b = n_b -1;
                    state_i = zeros(n_st,n_b);
                    for i=2:im-1
                        state_i(:,i-1) = state(:,i);
                    end
                    state_i(:,n_b) = [state(1,im);state(2,im)+state(2,implus);1];
                    state = state_i;
                    %updating t_nxt
                    t_nxt_i = zeros(1,n_b);
                    for i=2:n_b
                        t_nxt_i(i-1) = t_nxt(i);
                    end
                    t_nxt_i(n_b) = tspent;%t_nxt(n_b + 1) + tspent;
                    %t_nxt_i(n_b) =  tspent;
                    t_nxt = t_nxt_i;
                    %updating hwA
                    hwA_i = zeros(1,n_b);
                    for i=2:n_b
                        hwA_i(i-1) = hwA(i);
                    end
                    hwA = hwA_i;
                    
                    %updating atstop

                    atstop_i = zeros(1,n_b);

                    for i=2:n_b
                        atstop_i(i-1) = atstop(i);
                    end
                    atstop_i(n_b) = 1;
                    atstop = atstop_i;

                    %updating lapass
                    %as the rear bus has to skip the stop state(1,im) we
                    %have to add people willing to deboard in the lapass 
                    %pw = binornd((load_rj-lapass(1)),a_par(state(1,n_b)));
                    lapass_i = zeros(1,n_b);
                    for i =2:n_b
                        lapass_i(i-1) = lapass(i);
                    end
                    lapass_i(n_b) = lapass(1)+ lapass(n_b+1) + lap_ac;% + pw;
                    lapass = lapass_i;

                    %In this case the bus 2 will correspond to state(:,1)
                    ord_i = ord;
                    ord(1) = ord_i(tot_mod);
                    ord(2:end) = ord_i(1:tot_mod-1);
                else
                    n_b = n_b -1;
    
                    %Updating state
                    state_i = zeros(n_st,n_b);
                    for i=1:im-1
                        state_i(:,i) = state(:,i);
                    end
                    state_i(:,im) = [state(1,im);state(2,im)+state(2,implus);1];
                    for i=im+1:n_b
                        state_i(:,i) = state(:,i+1);
                    end
                    state = state_i;

                    %Updating t_nxt
                    
                    t_nxt_i = zeros(1,n_b);
                    for i = 1:im-1
                        t_nxt_i(i) = t_nxt(i);
                    end
                    
                    t_nxt_i(im) = tspent;%t_nxt(im) + tspent;
                    %t_nxt_i(im) = tspent;
                    for i =im+1:n_b
                        t_nxt_i(i) = t_nxt(i+1);
                    end
                    t_nxt = t_nxt_i;

                    %updating hwA
                    hwA_i = zeros(1,n_b);
                    for i=1:im-1
                        hwA_i(i) = hwA(i);
                    end
                    for i=im+1:n_b
                        hwA_i(i)= hwA(i+1);
                    end
                    hwA= hwA_i;

                    %Updaing atstop
                    atstop_i = zeros(1,n_b);
                    for i = 1:im-1
                        atstop_i(i) = atstop(i);
                    end
    
                    atstop_i(im) = 1;  %this was 0 before
    
                    for i=im+1:n_b
                        atstop_i(i) = atstop(i+1);
                    end
                    atstop = atstop_i;

                    %updating lapass
                    %pw = binornd((load_rj-lapass(implus)),a_par(state(1,im)));
                    
                    lapass_i = zeros(1,n_b);
                    for i=1:im-1
                        lapass_i(i) = lapass(i);
                    end
                    lapass_i(im) = lapass(im) + lapass(im+1) + lap_ac;% + pw;
                    for i=im+1:n_b
                        lapass_i(i) = lapass(i+1);
                    end
                    lapass = lapass_i;
                end   
        %If join is not feasible  
            else %modefy this to avoid overtake
                hwA(im) = time - nbtime(state(1,im));
                nbtime(state(1,im)) = time;
                snbcount = snbcount + 1;
                action = 'nextbs';              
                atstop(im) = 0;
                next_stop = true;
                eps_tr = gamrnd(gam_k, gam_theta) - gam_k*gam_theta;
                ct = dis_stp(state(1,im))/v_bus + eps_tr;
                mxtraveltime = 0;
                for bus = 1:n_b %this takes care of overtaking
                    if state(1,bus) == state(1,im) && atstop(bus) == 0
                        if t_nxt(bus) > mxtraveltime
                            mxtraveltime = t_nxt(bus);
                        end
                    end
                end

                ct =  max(ct, mxtraveltime) + 0.1; %adding 1 to ensure the module behind reaches late than the module in the front
                tspent = ct;
                t_nxt(im) = tspent;
            end          
    end
%     if ex_wt > 0 && imr == 2
%         fprintf('skip happened with ex_wt > 0 at gencount : %i \n', gencount)
%     end
    
    if time < simtime*3600
        if int32(time) == 0
            if Pa_pre_cum(1) < pb_cum 
                disp('khatra1')
                disp(pb_cum)
                disp(Pa_pre_cum(int16(time)))
            
            end
        else
            if Pa_pre_cum(ceil(time)) < pb_cum 
                disp('khatra2')
                disp(gencount)
                disp(pb_cum)
                disp(Pa_pre_cum(int16(time)))          
            end
        end
        if pa_cum < pb_cum
            disp('khatra1')
        end
    end 
    %Now we construct arrays that consist the bus locations at different
    %time instants
    Pa = [Pa pa];
%    Hw = [Hw hw(im)];
    Pa_cum = [Pa_cum pa_cum];
    Pd_cum = [Pd_cum pd_cum];
    Pb_cum = [Pb_cum pb_cum];
    Pw_cum = [Pw_cum pw_cum];
    Pe_cum = [Pe_cum pe_cum];
    mod_loc = zeros(1,tot_mod);
    atstopM = zeros(1,tot_mod);
    i = 1;
    for j = 1:size(state,2)
        if state(3,j) == 1
            mod_loc(i) = state(1,j);
            mod_loc(i+1) = state(1,j);
            atstopM(i) = atstop(j);
            atstopM(i+1) = atstop(j);
            i = i + 2;
        else
            mod_loc(i) = state(1,j);
            atstopM(i) = atstop(j);
            i = i + 1;
        end
    end
            
    for i=1:tot_mod
        if size(locs{i},2) > 0
            if atstopM(ord(i)) == 1 || locs{i}(size(locs{i},2)) ~= mod_loc(ord(i))
                locs{i} = [locs{i} mod_loc(ord(i))];
                mod_time{i} = [mod_time{i} time];
            end
        else
            locs{i} = [locs{i} mod_loc(ord(i))];
            mod_time{i} = [mod_time{i} time];
        end
    end


%     if n_b == 2
%         bus_loc = [state(1,1) state(1,1) state(1,2) state(1,2)];
%         atstopf = [atstop(1) atstop(1) atstop(2) atstop(2)];
%     elseif n_b == 3
%         if state(3,1) == 1
%             bus_loc = [state(1,1) state(1,1) state(1,2) state(1,3)];
%             atstopf = [atstop(1) atstop(1) atstop(2) atstop(3)];
%         elseif state(3,2) == 1
%             bus_loc = [state(1,1) state(1,2) state(1,2) state(1,3)];
%             atstopf = [atstop(1) atstop(2) atstop(2) atstop(3)];
%         elseif state(3,3) == 1
%             bus_loc = [state(1,1) state(1,2) state(1,3) state(1,3)];
%             atstopf = [atstop(1) atstop(2) atstop(3) atstop(3)];
%         end
%     elseif n_b == 4
%         bus_loc = state(1,:);
%         atstopf = atstop;
%     end
%     locf1 = [locf1 bus_loc(ord(1))];
     
%     if atstopf(ord(1)) == 1 || bus_loc(ord (1)) ~= loc1(size(loc1,2))
%         loc1 = [loc1 bus_loc(ord(1))];
%         Time1 = [Time1 time];
%     end
%     if atstopf(ord(2)) == 1 || bus_loc(ord(2)) ~= loc2(size(loc2,2))
%         loc2 = [loc2 bus_loc(ord(2))];
%         Time2 = [Time2 time];
%     end
% 
%     if atstopf(ord(3)) == 1 || bus_loc(ord(3)) ~= loc3(size(loc3,2))
%         loc3 = [loc3 bus_loc(ord(3))];
%         Time3 = [Time3 time];
%     end
%     if atstopf(ord(4)) == 1 || bus_loc(ord(4)) ~= loc4(size(loc4,2))
%         loc4 = [loc4 bus_loc(ord(4))];
%         Time4 = [Time4 time];
%     end

    if time > simtime*60*60
        break
    end
%     if time > 1500 && time < 1600
%         disp(gencount)
%         disp(1:24)
%         disp(state)
%         disp(atstop)
%         disp(t_nxt)
%     end
    %count = count + 1;
%     if gencount == 30
%         fprintf('pbs : %f \n', pbsr)
%         fprintf('pds : %f \n', pdsr)
%         fprintf('hway :%f \n',hway)
%     end
%     if gencount < 10
%         disp(table(state', atstop', t_nxt'))
%     end
    
%     if time < 1300
%         if imr == 1
%             stcn = stcn + 1;
%         elseif imr == 2
%             skcn = skcn + 1;
%         elseif imr == 3
%             spcn  = spcn + 1;
%         elseif imr == 4
%             jncn = jncn + 1;
%         elseif imr == 5
%             nbcn = nbcn + 1;
%         end
%     end
%       if gencount > 500 && gencount < 536
%           if im == 2
%               fprintf('module 2 gencount : %i \n', gencount)
%           end
%       end
%          if gencount < 10
%              fprintf('count : %i \n', count)
%          end
%          if gencount ==  415
%             fprintf('%s action on module %d \n', action, im)
%             fprintf('action comb : %s \n', optSeq)
%             fprintf('module %d spends time %f with extra time %f \n', im, tspent, ex_wt)
%             fprintf('state at t = %d \n', time)
%             fprintf('new added time is %d \n',M)
%             %fprintf('pdsr : %i, pbsr : %i, ex_wt : %i, load : %i, a_par : %i, lapass(im) : %i \n', pdsr, pbsr, ex_wt, state(2,im), a_par(state(1,im)), lapass(im))
%             disp(1:n_b)
%             disp(state)
%             disp(ord)
%             disp(atstop)
%             disp('tnxt after')    
%             disp(t_nxt)
%             disp(lapass)
%             disp(gencount)
%         end

        if pd_cum < pe_cum
            disp('scam')
            disp(tw_pass)
            disp(w_pass)
            break
            
        end
        if pa_cum < pb_cum
            disp('fuck my life')
        end
end 
i = 1;
while Time(i) < etime
    i = i+1;    
end
eind = i-1;

for i=3:size(Time,2)-1
    pWait = [pWait Pa_pre_cum(round(Time(i))) - Pb_cum(i)];
end
stopTimeA = stopTimeA/(stcount + spcount);

%finding the waiting and in-vehicle passengers
PInVeh = Pb_cum - Pd_cum;

% PWait = zeros(1,size(Pb_cum,2));
% for i = 1: size(Pb_cum,2)
%     PWait(i) = Pa_pre_cum(round(Time(i))) - Pb_cum;
% end
%now finding the area under the curves
interv = sind:eind;
Pb_cum_int = trapz(Time(interv), Pb_cum(interv));
Pd_cum_int = trapz(Time(interv), Pd_cum(interv));
Pa_cum_int = trapz(Time(interv), Pa_cum(interv));
Pw_cum_int = trapz(Time(interv), Pw_cum(interv));
Pe_cum_int = trapz(Time(interv), Pe_cum(interv));
Pa_pre_cum_int = trapz(Pa_pre_cum(ceil(stime):ceil(etime)));

avg_inveh = (Pb_cum_int - Pd_cum_int)/(Pb_cum(eind)- Pb_cum(sind));
avg_wait = (Pa_cum_int - Pb_cum_int)/Pa_cum(size(Pa_cum,2));
avg_walk = Pw_cum_int/Pd_cum(size(Pd_cum,2)); %this formulation is wrong
avg_walk1 = (Pd_cum_int - Pe_cum_int)/(Pd_cum(eind) - Pb_cum(sind));
avg_wait1 = (Pa_pre_cum_int - Pb_cum_int)/(Pa_pre_cum(ceil(etime)) - Pa_pre_cum(ceil(stime)));
fprintf('Average in-vehicle time (m) : %f \n', avg_inveh/60)
fprintf('Average waiting time (m) : %f \n', avg_wait1/60)
fprintf('Average walk time (m) : %f \n', avg_walk1/60)
fprintf('Policy cost (m) : %f \n', (w_wait*avg_wait1 + avg_inveh + w_walk*avg_walk1)/60)
% fprintf('final time : %f hr \n', time/3600)
% fprintf('join count : %i \n', sjcount)
% fprintf('skip count : %i \n', skcount)
% fprintf('stop count : %i \n', stcount)
% fprintf('sind : %i, eind : %i \n', sind, eind)
% fprintf('stime : %f, etime : %f \n', stime, etime)
% fprintf('Average stopping time : %f \n', stopTimeA)
% fprintf('Code run time : %f \n', toc)
% fprintf('Cycle time : %f \n', cycT)
fprintf('pb_cum_int : %f, pd_cum_int : %f, boardings : %f \n', Pb_cum_int, Pd_cum_int, Pb_cum(eind)- Pb_cum(sind))
%preparing to create plots
ilist = cell(1,tot_mod);
ilist_c = cell(1,tot_mod);
% for i=1:tot_mod
%     ilist{i}= 0;
% end
seInd = zeros(2,tot_mod);
% for i=1:tot_mod
%     for j = 1: size(locs{i},2) - 1
%         if mod_time{i}(j) > stime
%             seInd(1,i) = j;
%             break
%         end
for i=1:tot_mod
    [~, seInd(1,i)] = min(abs(mod_time{i} - stime));
    [~, seInd(2,i)] = min(abs(mod_time{i} - etime));
    ilist{i} = [ilist{i} seInd(1,i)];
    ilist_c{i} = [ilist_c{i} 1];
end

for i=1:tot_mod
    for j = 1: size(locs{i},2) - 1
        %if locs{i}(j) == n_s && locs{i}(j+1) == 1
        if locs{i}(j) > locs{i}(j+1) && mod_time{i}(j) >= stime && mod_time{i}(j) <= etime
            ilist{i} = [ilist{i} j];
        end
        if locs{i}(j) > locs{i}(j+1)
            ilist_c{i} = [ilist_c{i} j];
        end
    end
end
for i=1:tot_mod
    ilist{i} = [ilist{i} seInd(2,i)];
end
colA = ["#FF0000",'#FFFF00','#00EAFF','#AA00FF','#FF7F00','#BFFF00','#0095FF','#FF00AA','#FFD400','#6AFF00','#0040FF','#EDB9B9','#B9D7ED','#E7E9B9','#DCB9ED','#B9EDE0','#8F2323','#23628F','#8F6A23','#6B238F','#4F8F23','#000000','#737373','#CCCCCC'];
figure(10) %plotting only the evaluation period
hold on
for i = 1: tot_mod
    for j = 1:size(ilist{i},2)-1
        i1 = ilist{i}(j)+1; i2 = ilist{i}(j+1);
        plot(mod_time{i}(i1:i2), locs{i}(i1:i2), 'Color', colA(i), LineWidth=1.2)
    end
end
xlabel('Time (s)', FontSize=12)
ylabel('Stops', FontSize=12)
yticks(1:n_s)
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';
%plotting for all time
figure(11)
hold on
for i = 1: tot_mod
    for j = 1:size(ilist_c{i},2)-1
        i1 = ilist_c{i}(j)+1; i2 = ilist_c{i}(j+1);
        plot(mod_time{i}(i1:i2), locs{i}(i1:i2), 'Color', colA(i), LineWidth=1.2)
    end
end
ax.FontSize = 12;
yticks(1:n_s)
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';
set(gca,"FontSize",12)



figure(6)
%plot(Time, Pa_cum)
hold on
plot(Time_pre, Pa_pre_cum, LineWidth=1.2)
plot(Time, Pb_cum, LineWidth=1.2)
plot(Time, Pd_cum, LineWidth=1.2) 
plot(Time, Pe_cum, LineWidth=1.2)
plot(Time(3:size(Time,2)-1), pWait, LineWidth=1.2)
plot(Time, PInVeh, LineWidth=1.2)
xline(stime, '--')
xline(etime, '--')
xlabel('Time (s)', 'FontSize',12)
ylabel('Cumulative count', FontSize=12)
legend( 'arrival','boardings','alightings','exiting', 'waiting', 'in-vehicle', 'Location','northwest', 'FontSize', 11)

%ylim([0 5000])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

figure(8)
plot(Time, pVeh)
hold on
plot(Time(3:size(Time,2)-1), pWait)
legend('In vehicle passengers', 'waiting passengers')
grid on

figure(9)
plot(Time(sind:eind), pVeh(sind:eind))
hold on
plot(Time(sind:eind), pWait(sind:eind))

% disp(stime)
% disp(etime)
% 
t_sta_comb = 4500; t_sto_comb = 7500;
[~, i_sta] = min(abs(t_sta_comb - Time));
[~, i_sto] = min(abs(t_sto_comb - Time));
P_veh_p = PInVeh(i_sta:i_sto); P_wait_p = pWait(i_sta:i_sto);
P_walk_p = Pd_cum(i_sta:i_sto) - Pe_cum(i_sta:i_sto);
t_p = Time(i_sta:i_sto);
figure(12)
plot(t_p, P_veh_p)
hold on
plot(t_p, P_wait_p)
grid on