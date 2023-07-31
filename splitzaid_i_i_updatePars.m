%%Code to simulate the variable agent algorithm
rng(95) %13, 20, 23
close all
%This code is using the most recent headway function that is headwaycuu
%this code also accounts the leftover passengers if the bus skips the stop
%but at this point of time we are assming that the bus not skipping stops
%successively
%dbstop if naninf
simtime = 3; %simulation time in hours
n_s = 20; %Number of Stations
n_b = 12; %Number of bus, 1 bus = 2 modules
tot_mod = 2*n_b;
n_st = 3; %Number of states 
n_a = 4; % Number of actions bus can take
%a_par = rand(1,n_s); %These are Ps values for the stops
almu = 2/n_s;
a_par = normrnd(almu, almu/10, 1, n_s);
%arr_par = rand(1,n_s)/60*2; % Assuming on 90 passengers arrive in 30 minutes
hdem = 1500;
arrmu = hdem/3600/n_s; arrsigma = arrmu/10;
arr_par = normrnd(arrmu, arrsigma, 1,n_s);
atstop = zeros(1,n_b); %These are the flags which will be 1
% if the corresponding bus is at stop and 0 it it is on the road
fixdt = 20; %fixed time lost per stop
dis_stp = 333*(rand(1,n_s) + 1); %Distance between stops distributed between 300 - 600 meters
dis_stp = normrnd(400, 40, 1,n_s);
v_bus = 20*5/18; % Speed of bus 20 Km/h
w_wait = 2.1; w_walk = 2.2;% weights of walk time and wait time in final cost
lapass = zeros(1,n_b);
lpass = zeros(1,n_s);
w_pass = zeros(4,n_s); %to store the number of walking passengers to a particular stop
%at max there can be 4 types of passengers walking towards a particular
%stop
tw_pass = zeros(4,n_s); %time to reach the next stop
cap_bus = 40;
unit_cap = cap_bus/2;
v_pas = 4.5*5/18; %Passenger speed in Km/h
t_bo = 4; %boarding time per passenger in seconds
t_al = 3; %Alighting time per passenger in seconds
scamcount = 0;
gam_k = 2; gam_theta = 5.1; %parameters of Gamma rv
gamma = 1.5;
bet = 0.15;
hwt_i = (sum(dis_stp)/v_bus + n_s*fixdt)/(n_b - (t_al + t_bo)*arrmu*n_s);
hwt = hwt_i;
cumtime = zeros(1,n_b);
splittime = 1; %splittime being 1 second -- Introducing this constant to avoid bus overtake
stload = floor(arrmu*n_s/2*hwt_i); 
% state0 = zeros(n_st,n_b);
% state0(n_st,:) = ones(1,n_b);
% state0(1,:) = ones(1,n_b);

%state0(2,:) = stload;
state0 = zeros(n_st,n_b);
state0(n_st,:) = ones(1,n_b); %uncomment if the modules start in the
%attached state
state0(1,:) = ones(1,n_b)*n_s;
state0(2,:) = stload;
state = state0;
stop = false;
split = false;
join = false;
justdeb = false;
waiting = false;
nextbus = false;

%Busses will depart with headay of 90 seconds
%We will work on the basis of time a bus required to reach the next stop.
%It will we stored in the array n_xt. Note that this array will change as
%the code proceeds and its size will also change.
t_nxt = zeros(1,n_b); %this variable stores the time bus will require to reach the upcoming stop
for i=1:n_b  %At start each bus leaves with an interval of 90 seconds
    t_nxt(i) = hwt_i*(i-1) ;%+ dis_stp(1)/v_bus;
end
%disp(t_nxt)
l_action = zeros(n_st,n_b); %last action taken by th bus
papbcum = 0;
papbskip = 0;
count = 1; %This will account the number of rounds
%consec_t = dis_stp/v_bus;  %We are assuming that, in the first round the passengers begin to arrive at the stop after the bus leaves the previous stop 
gencount = 0;
skcount = 0;
stcount = 0;
sjcount = 0;
spcount = 0;
snbcount = 0;
alcount = 0;
scamnb = 0;
T = 0;
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
ord = 1:tot_mod; %will keep track of the order of the buses. This is basically 
%the mapping between the state index and the bus number
%Creating arriving passenger array 
% pa_pre = poissrnd(arrmu*n_s, 1, simtime*3600);
Pa_pre_cum = zeros(1,simtime*3600+1);
% pa_pre_cum = stload*n_b;
Time_pre = 0:simtime*3600;
% pa_pre_cum(1) = pa_pre_cum;
% for i = 1:size(pa_pre,2)
%     pa_pre_cum = pa_pre_cum + pa_pre(i);
%     Pa_pre_cum(i+1) = pa_pre_cum;
% end
%generating arrivals for each stop
tflst = sum(dis_stp)/v_bus + (t_al + t_bo)*hwt_i*sum(arr_par) + fixdt*n_s;
tfstops = zeros(1,n_s);
for i = 1:n_s -1
    tfstops(i+1) = tfstops(i) + dis_stp(i)/v_bus + (t_al + t_bo)*hwt_i*arr_par(i) + fixdt;
end
Pa_all_cum = zeros(n_s, simtime*3600+1);
nwpasscount = 0;
% for i = 1:n_s
%     %passnger coming to ith stop in one target headway
%     passhwt = poissrnd(arr_par(i)*hwt_i);
%     for k = 1:n_s-1
%         for j=floor(tfstops(k))+1:floor(tfstops(k+1))
%             if Pa_all_cum(i,j) >= passhwt && i > k
%                 newpass = 0;
%                 nwpasscount = nwpasscount + 1;
%             else
%                 newpass = poissrnd(arr_par(i));
%             end
%             Pa_all_cum(i,j+1) = Pa_all_cum(i,j) + newpass;
%         end
%         for j = floor(tfstops(k+1))+1: simtime*3600
%             Pa_all_cum(i,j+1) = Pa_all_cum(i,j) + poissrnd(arr_par(i));
%         end
%     end
% end

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
singcount = 0;
rotcount = 0;
locs = cell(1,tot_mod); %to store the module location 
for i=1:tot_mod
    locs{i} = [];
end
mod_time = cell(1, tot_mod); %to store module time
for i=1:tot_mod
    mod_time{i} = [];
end
ns_reach = 0;
pVeh = 0;
while true
%     if n_b > 2
%         disp('attention')
%     end
    %array for identifying front and rear module of a bus
    
%     fr = zeros(1,n_b);
%     if n_b == 4
%         fr(1) = 1; fr(3) = 1;
%     elseif n_b == 3
%         if state(3,1) == 1  
%             fr(2) = 1;
%         else
%             fr(1) = 1;
%         end
%     end
    
    pVeh = [pVeh sum(state(2,:))];
    if state(1,1) == n_s %count contains whether the leading module has covered a circle or not
        ns_reach = ns_reach + 1;
    end
    if ns_reach > 1
        count = 0;
    end
    gencount = gencount + 1;


    [M,im] = min(t_nxt);  %This step finds which bus reaches/leaves the next stop first
    time = time + M;
    
    if state(1, 1) == n_s && im == 1
        rotcount = rotcount + 1;
    end
    if rotcount == 4
        sind = gencount;
        stime = time;
        etime = stime + 1*60*60;
    end

    if time == inf
        disp('stopping')
        break
    end
    T = [T M];
    %finding other indices with the minimum value
    imin = find(t_nxt == M);
    im = min(imin);
    fr = fR(state, im);
    implus = iplus(im,n_b);
    implusp = iplus(implus,n_b);
    for i=1:n_b
            if i ~= im
                t_nxt(i) = t_nxt(i) - t_nxt(im);  %nxt for im will be calculated as the stop time at the stop
            end       
    end

    for i= 1:4
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
    ex_wt = 0;
    if atstop(im) == 0

        if state(1,im) == n_s %Because our stops are circular n_s->1
            state(1,im) = 1;
        else
            state(1,im) = state(1,im) + 1; %state(1,im) is storing the last bus-stop number the bus has visited
        end

        if im ==1 && count == 1
            hway = hwt;
        else
            hway = time - reachtime(state(1,im));
        end
        
       
        acap = cap_bus - state(2,im);
    
        
        st_reach = state(1,im); %because we already have updated the stop in the begining of the loop
        b_st = [];
        ex_wt = 0; %extra wait the bus has to do because of already reched buses
        %following lines find the max waiting time among the buses already
        %at the stop where the current bus (im) is reached.
        for i=1:n_b
            if state(1,i) == st_reach && atstop(i) == 1 && fr ~= 0 
                b_st = [b_st i];
                if t_nxt(i) > ex_wt && t_nxt(i) ~= inf
                    ex_wt = t_nxt(i);
                end
            end
        end

    
        if state(3,im) == 1
%             if time < 12000 && time > 11000
%                 fprintf('hw : %f \n', hway)
%                 fprintf('hwt : %f \n', gamma*hwt)
%                 fprintf('gencount, time : %i, %f \n', gencount, time)
%             end
            if hway > gamma*hwt
                if im ==1
                    fprintf('split 1 gencount : %i \n', gencount)
                end
                action  = 'split';
                split = true;
                spcount = spcount +1;
                n_b = n_b +1; %The number of modules increase due to splitting
                %compute th number of passengers getting down at the current stop
                %this is not how zaid is doiing things i mean the
                %deboarding passengers count as done following
                pds_c = binornd((state(2,im)-lapass(im)),a_par(state(1,im))) + lapass(im);
                if pds_c>= state(2,im)
                    alcount = alcount + 1;
                end
                %Assuming that only the passengers getting down at the current stop
                %are in the rear module
                stmin = iminus(state(1,im),n_s);
                for iw=1:4
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
                load = load_r;
                pdsr = pds_c;
                pd_cum = pd_cum + pdsr;
                pe_cum = pe_cum + pdsr - lapass(im);
 
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
                
           
                f_ind= floor(reachtime(state(1,im))); s_ind = floor(time);
                if s_ind >= size(Pa_all_cum,2)
                        s_ind = size(Pa_all_cum,2)-1;
                end
                padec(state(1,im), f_ind+1:s_ind+1) = 1;
                pa = Pa_all_cum(state(1,im),s_ind+1) - Pa_all_cum(state(1,im),f_ind+1);
                pa_cum = pa_cum + pa;
                
                pbsr = min(pa + lpass(state(1,im)), unit_cap - state(2,im+1));
              
                reachtime(state(1,im)) = time;
                pb_cum = pb_cum + pbsr;
                lpass(state(1,im)) = pa + lpass(state(1,im)) - pbsr; %assuming the leftover passengers gets priority boarding and they are never more than unitcap
                papbcum = papbcum + pa-pbsr;
    %             if gencount == 498
    %                 disp(state(2,im+1))
    %                 disp(pbsr)
    %                 fprintf('bpass  = %f \n',bpass)
    %                 fprintf('headway = %f \n',hw(im))
    %             end
                state(2,im+1) = state(2,im+1) + pbsr;
                tspent = max(t_al*pdsr, t_bo*pbsr) + fixdt + ex_wt + splittime;
               
                
                %Updating t_nxt as the number of buses has changed
                t_nxt_i = zeros(1,n_b);
                for i = 1:im-1
                    t_nxt_i(i) = t_nxt(i);
                end
                ct = dis_stp(state(1,im))/v_bus + gamrnd(gam_k, gam_theta) - gam_k*gam_theta;  %this tspent can be the t_nxt two lines above

                t_nxt_i(im) = ct + ex_wt + splittime; %As the leading module does not stop.
                t_nxt_i(im + 1) = tspent; %time to leave the stop
                
                for i = im+2:n_b
                    t_nxt_i(i) = t_nxt(i-1);
                end
                t_nxt = t_nxt_i;
            else
                action = 'stop';
                stop = true;
                stcount = stcount + 1;
                load = state(2,im);
                pdsr = lapass(im) + binornd((state(2,im) - lapass(im)),a_par(state(1,im)));
                pd_cum = pd_cum + pdsr;
                pe_cum  = pe_cum + pdsr - lapass(im);
                stmin = iminus(state(1,im),n_s);
                for iw=1:4
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

                lapass(im) = 0;
                
                state(2,im) = state(2,im) - pdsr;
%                 if gencount < 10
%                     disp(state)
%                 end
                
                
                f_ind= floor(reachtime(state(1,im))); s_ind = floor(time);
                if s_ind >= size(Pa_all_cum,2)
                    s_ind = size(Pa_all_cum,2)-1;
                end
                padec(state(1,im), f_ind+1:s_ind+1) = 1;
                pa = Pa_all_cum(state(1,im),s_ind+1) - Pa_all_cum(state(1,im),f_ind+1);
                pa_cum = pa_cum + pa;
                reachtime(state(1,im)) = time;
%                 if gencount == 71
%                     disp(lpass)
%                     disp(pa)
%                 end
                
                pbsr = min(pa + lpass(state(1,im)), cap_bus - state(2,im));
                pb_cum = pb_cum + pbsr;
                lpass(state(1,im)) = pa + lpass(state(1,im)) - pbsr;
                papbcum = papbcum + pa-pbsr;
                state(2,im) = state(2,im) + pbsr;
                
                tspent = max(t_al*pdsr , t_bo*pbsr) + fixdt + ex_wt;
                %Updating t_nxt, atstop, 
                atstop(im) = 1;
                t_nxt(im) = tspent;
                %Not updating the headway as we are neglecting the time
                %spent at the stop in the headway calculation
            end


        
    
        else  %if the bus is splitted already then the competition will be in stop and skip
            singcount = singcount + 1;
            if fr == 0
                action = 'just deboard';
                %disp(action)
                load = state(2,im);
                justdeb = true;
                stcount = stcount + 1;
                pdsr = lapass(im) + binornd((state(2,im) - lapass(im)),a_par(state(1,im)));
                pd_cum = pd_cum + pdsr;
                pe_cum  = pe_cum + pdsr - lapass(im);
                stmin = iminus(state(1,im),n_s);
                for iw=1:4
                    if w_pass(iw,stmin) == 0
                        break
                    end
                end
                if prod(w_pass(:,stmin)) ~=0
                    disp('damn... this is not working')
                    break
                end
                w_pass(iw,stmin) = lapass(im); %If the deboarding passengers is more than unit capacity then the left over passengers gets priority for deboarding
                if w_pass(iw,stmin) ~=0 %here we are assuming that the walking passengers are walking towards the last stop but inreality there can be some walking passengers who want to walk to two or more stops before but they are rare so ignoring that
                    tw_pass(iw,stmin) = dis_stp(stmin)/v_pas;
                end

                lapass(im) = 0;
                
                state(2,im) = state(2,im) - pdsr;
%                 if gencount < 10
%                     disp(state)
%                 end
                
                
                f_ind= floor(reachtime(state(1,im))); s_ind = floor(time);
                if s_ind >= size(Pa_all_cum,2)
                    s_ind = size(Pa_all_cum,2)-1;
                end
                padec(state(1,im), f_ind+1:s_ind+1) = 1;
                pa = Pa_all_cum(state(1,im),s_ind+1) - Pa_all_cum(state(1,im),f_ind+1);
                pa_cum = pa_cum + pa;
                reachtime(state(1,im)) = time;
                pbsr = 0;
                lpass(state(1,im)) = pa + lpass(state(1,im)) - pbsr;
                tspent = t_al*pdsr  + fixdt + ex_wt;
                %Updating t_nxt, atstop, 
                atstop(im) = 1;
                t_nxt(im) = tspent;
                %Not updating the headway as we are neglecting the time
                %spent at the stop in the headway calculation
            else
                action = 'stop';
                %disp(action)
                stop = true;
                stcount = stcount + 1;
                load = state(2,im);
                pdsr = lapass(im) + binornd((state(2,im) - lapass(im)),a_par(state(1,im)));
                pd_cum = pd_cum + pdsr;
                pe_cum  = pe_cum + pdsr - lapass(im);
                stmin = iminus(state(1,im),n_s);
                for iw=1:4
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

                lapass(im) = 0;
                
                state(2,im) = state(2,im) - pdsr;

                
                f_ind= floor(reachtime(state(1,im))); s_ind = floor(time);
                if s_ind >= size(Pa_all_cum,2)
                    s_ind = size(Pa_all_cum,2)-1;
                end
                padec(state(1,im), f_ind+1:s_ind+1) = 1;
                pa = Pa_all_cum(state(1,im),s_ind+1) - Pa_all_cum(state(1,im),f_ind+1);
                pa_cum = pa_cum + pa;
                reachtime(state(1,im)) = time;
%                 if gencount == 71
%                     disp(lpass)
%                     disp(pa)
%                 end
                
                pbsr = min(pa + lpass(state(1,im)), unit_cap - state(2,im));
                pb_cum = pb_cum + pbsr;
                lpass(state(1,im)) = pa + lpass(state(1,im)) - pbsr;
                papbcum = papbcum + pa-pbsr;
                state(2,im) = state(2,im) + pbsr;
                
                tspent = max(t_al*pdsr , t_bo*pbsr) + fixdt + ex_wt;
                %Updating t_nxt, atstop, 
                atstop(im) = 1;
                t_nxt(im) = tspent;

                 
            
            %reachtime(state(1,im)) = time;
            end 
        end
     
    
    else %If the bus is not at the stop we will evalueate join possibility
        
        if im == 1 && count == 1
            hway=hwt_i;                                                                                                                                                                                         y = hwt;
        else
            hway = time - leavetime(state(1,im));
        end
        leavetime(state(1,im)) = time;         
        
        if state(3,im) == 0  
            if fr == 1
                t_nxt(im) = inf;
                action  = 'wait';
                wait = true;
                waiting = true;
                %write code so that the module waits here till the rear
                %module comes
            else %join imth module with the module in front
                if state(1,im) == state(1,im -1) && atstop(im-1) == 1
                    sjcount = sjcount + 1;                    
                    action = 'join';
                    join = true;
                    %disp(action)
                    n_b = n_b -1;
    
                    %Updating state
                    state_i = zeros(n_st,n_b);
                    for i=1:im-2
                        state_i(:,i) = state(:,i);
                    end
                    state_i(:,im-1) = [state(1,im-1);state(2,im)+state(2,im-1);1];
                    for i=im:n_b
                        state_i(:,i) = state(:,i+1);
                    end
                    state = state_i;

                    %Updating t_nxt
                    
                    t_nxt_i = zeros(1,n_b);
                    for i = 1:im-2
                        t_nxt_i(i) = t_nxt(i);
                    end
                    
                    t_nxt_i(im-1) = 0;
    
                    for i =im:n_b
                        t_nxt_i(i) = t_nxt(i+1);
                    end
                    t_nxt = t_nxt_i;

                    %Updaing atstop
                    atstop_i = zeros(1,n_b);
                    for i = 1:im-2
                        atstop_i(i) = atstop(i);
                    end
    
                    atstop_i(im-1) = 1;  %this was 0 before
    
                    for i=im:n_b
                        atstop_i(i) = atstop(i+1);
                    end
                    atstop = atstop_i;

                    %updating lapass
                    lapass_i = zeros(1,n_b);
                    for i=1:im-2
                        lapass_i(i) = lapass(i);
                    end
                    lapass_i(im-1) = lapass(im) + lapass(im-1);
                    for i=im:n_b
                        lapass_i(i) = lapass(i+1);
                    end
                    lapass = lapass_i;
                    waiting = false;
                else
                    scamnb = scamnb + 1;
                    snbcount = snbcount + 1;
                    action = 'nextbs';
                    nextbus = true;
                    atstop(im) = 0;
               
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

                    ct =  max(ct, mxtraveltime) + 1; %adding 1 to ensure the module behind reaches late than the module in the front
                    tspent = ct;
                    t_nxt(im) = tspent;
                end
            end                
        else
                snbcount = snbcount + 1;
                action = 'nextbs';
                nextbus = true;
                tspent = dis_stp(state(1,im))/v_bus + gamrnd(gam_k, gam_theta) - gam_k*gam_theta;
                t_nxt(im) = tspent;
                atstop(im) = 0;         
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
    atstopt = atstop; %atstop temp
    if nextbus
        atstopt(im) = 1;
    end
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

    Time = [Time time];
    if time > simtime*60*60
        break
    end
    if time > 250 && time < 270
        fprintf('action : %s \n', action)
        fprintf('gencount : %i, count : %i \n', gencount,count)
        fprintf('im : %i, hway : %f \n', im, hway)
        fprintf('time : %f \n', time)

    end
% 
    if time > 250 && time < 270
        fprintf('gamma*hwt : %f, hw : %f \n', gamma*hwt, hway)
        fprintf('%s action on module %d \n', action, im)
        fprintf('module %d spends time %f with extra time %f \n', im, tspent, ex_wt)
        fprintf('state at t = %d \n', time)
        fprintf('new added time is %d \n',M)
%         fprintf('No. of boarding passengers : %i \n', pbsr)
%         fprintf('lpass : %i\n', lpass(1))
%         fprintf('No. of deboarding passengers : %i \n',pdsr)
        disp(state)       
        %disp(ord)
        disp(atstop)
        disp('tnxt after')
        disp(t_nxt)
        disp(gencount)
    end

        if pd_cum < pe_cum
            disp('scam')
            disp(tw_pass)
            disp(w_pass)
            break
            
        end
        if pa_cum < pb_cum
            disp('fuck my life')
        end
        stop = false;
        split = false;
        join = false;
        justdeb = false;
        nextbus = false;       
end 
timin = inf;
for i = size(Time,2)
        if abs(etime - Time(i)) < timin
            timin = abs(etime - Time(i));
            eind = i;
        end
end
i = 1;
while Time(i) < etime
    i = i+1;    
end
eind = i-1;

pWait = [];
for i=3:size(Time,2)-1
    pWait = [pWait Pa_pre_cum(round(Time(i))) - Pb_cum(i)];
end
%in-vehicle passengers
PInVeh = Pb_cum - Pd_cum;

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
fprintf('No. of passenger boarding : %i \n', Pb_cum(eind)- Pb_cum(sind))

% fprintf('final time : %f hr \n', time/3600)
% fprintf('join count : %i \n', sjcount)
% fprintf('skip count : %i \n', skcount)
% fprintf('stop count : %i \n', stcount)
fprintf('pb_cum_int : %f, pd_cum_int : %f, boardings : %f \n', Pb_cum_int, Pd_cum_int, Pb_cum(eind)- Pb_cum(sind))
Pa_pre_cum1 = Pa_pre_cum; 

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
xlabel('Time (s)')
ylabel('Stops')
yticks(1:n_s)
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';
ax.FontSize = 12;
%yticks(1:n_s)
%plotting for all time
figure(11)
hold on
for i = 1: tot_mod
    for j = 1:size(ilist_c{i},2)-1
        i1 = ilist_c{i}(j)+1; i2 = ilist_c{i}(j+1);
        plot(mod_time{i}(i1:i2), locs{i}(i1:i2), 'Color', colA(i), LineWidth=1.2)
    end
end
yticks(1:n_s)
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';


figure(6)
%plot(Time, Pa_cum)
hold on
plot(Time_pre, Pa_pre_cum, LineWidth=1.2)
plot(Time, Pb_cum, LineWidth=1.2)
plot(Time, Pd_cum, LineWidth=1.2) 
plot(Time, Pe_cum, LineWidth=1.2)
plot(Time(3:size(Time,2)-1), pWait, LineWidth=1.2)
plot(Time, Pb_cum - Pd_cum, LineWidth=1.2)
xline(stime, '--')
xline(etime, '--')
xlabel('Time (s)', 'FontSize',12)
ylabel('Cumulative count', FontSize=12)
legend( 'arrival','boardings','alightings','exiting', 'waiting', 'in-vehicle', 'Location','northwest', 'FontSize', 11)

ylim([0 5000])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% figure(8)
% plot(Time, pVeh)
% hold on
% plot(Time(3:size(Time,2)-1), pWait)
% legend('In vehicle passengers', 'waiting passengers')
% grid on
% disp(stime)
% disp(etime)
t_sta_comb = 4500; t_sto_comb = 7500;
[~, i_sta] = min(abs(t_sta_comb - Time));
[~, i_sto] = min(abs(t_sto_comb - Time));

P_veh_z = pVeh(i_sta:i_sto); P_wait_z = pWait(i_sta:i_sto);
t_z = Time(i_sta:i_sto);
figure(8)
plot(t_z, P_veh_z)
hold on
plot(t_z, P_wait_z)
grid on
