%check the value of nom variable before running the code
nsim = 100;
costHor1 = zeros(1,nsim);
jAve = 0; 
skAve = 0;
stAve = 0;
nbAve = 0;
spAve = 0;
nom = true;
for simn = 1:nsim
    rng(simn)    
    simtime = 3; %simulation time in hours
    n_s = 20; %Number of Stations
    n_b = 12*2; %Number of bus, 1 bus = 2 modules modules starts in the detached state
    nb_start = n_b;
    tot_mod = nb_start;
    n_st = 3; %Number of states 
    n_a = 5; % Number of actions bus can take
    %a_par = rand(1,n_s); %These are Ps values for the stops
    hor = 1;
    fPas1 = 0; fPas2 = 1; %reward function parameters
    jP1 = 1; jP2 = 1;
    almu = 2/n_s;
    a_par = normrnd(almu, almu/10, 1, n_s); 
    %arr_par = rand(1,n_s)/60*2; % Assuming on 90 passengers arrive in 30 minutes
    hdem = 1500;
    arrmu = hdem/3600/n_s; arrsigma = arrmu/10;
    arr_par = normrnd(arrmu, arrsigma, 1,n_s);
    atstop = zeros(1,n_b); %These are the flags which will be 1
    % if the corresponding bus is at stop and 0 it it is on the road
    fixdt = 30; %fixed time lost per stop
    %dis_stp = 333*(rand(1,n_s) + 1); %Distance between stops distributed between 300 - 600 meters
    dis_stp = normrnd(500, 50, 1,n_s);
    v_bus = 20*5/18; % Speed of bus 20 Km/h
    w_wait = 2.1; w_walk = 2.2;% weights of walk time and wait time in final cost
    lapass = zeros(1,n_b);
    lpass = zeros(1,n_s);
    w_pass = zeros(n_b,n_s); %to store the number of walking passengers to a particular stop if the modules are leaving in the joined state, first parameter should be 2*n_b
    %at max there can be 4 types of passengers walking towards a particular
    %stop      
    Hway = zeros(1,n_b);
    tw_pass = zeros(2*n_b,n_s); %time to reach the next stop
    cap_bus = 100;
    unit_cap = cap_bus/2;
    v_pas = 5.4*5/18; %Passenger speed in Km/h
    t_bo = 5; %boarding time per passenger in seconds
    t_al = 2; %Alighting time per passenger in seconds
    scamcount = 0;
    gam_k = 2; gam_theta = 2; %parameters of Gamma rv
    gamma = 1.5;
    bet = 0.6;
    hwt_i = (sum(dis_stp)/v_bus + n_s*fixdt)/(n_b - (t_al + t_bo)*arrmu*n_s);
    hwt = hwt_i;
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
    ord = [1 2 3 4]; %will keep track of the order of the buses. This is basically 
    next_stop = false;
    rotCounted = true;
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
    % 
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
    
    rotcount = 0;
    ns_reach = 0; %# of times the bus reached the state n_s
    stcn = 0; spcn = 0; skcn = 0; jncn = 0; nbcn = 0; %action counters in particular time interval
    pVeh = 0; %number of passengers in vehicle
    pWait = []; %number of wiating passengers
    
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
    %         fprintf('rotcount : %i , time : %f \n', rotcount, time)
    %         disp(state)
            rotCounted = false;
        end
        if rotcount == 2 && ~rotCounted
            %disp(time)
            %disp('hi')
            sind = gencount;
            stime = time;
            etime = stime + 1*60*60;
            rotCounted = true;
        end
    
        T = [T M];
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
                hway = time - reachtime(iplus(state(1,im),n_s));
            end
            
            if state(1,im) == n_s %Because our stops are circular n_s->1
                state(1,im) = 1;
            else
                state(1,im) = state(1,im) + 1; %state(1,im) is storing the last bus-stop number the bus has visited
            end
            opt_act = rewsum_t(hor, state, t_nxt_sim, atstop, t_al, t_bo, arr_par, a_par, dis_stp, im, hway, v_bus, v_pas, lpass, lapass, cap_bus, hwt,ord, count, gencount, reachtime, leavetime, time,Hway, Pa_all_cum, fPas1, fPas2, jP1, jP2);
            %opt_act = rewsum_expt(hor, state, t_nxt_sim, atstop, t_al, t_bo, arr_par, a_par, dis_stp, im, hway, v_bus, v_pas, lpass, lapass, cap_bus, hwt,ord, count, gencount);
    
            opt_a = opt_act;
            bstop = state(1,im);
    %         fprintf('rewsum1 im : %i \n', im)
    %         disp('t_nxt 1')
    %         disp(t_nxt)
            %disp('state passed in rewsum')
            %disp(state)
            %fprintf('opt_act1 : %s \n', opt_act1)
            %fprintf('opt_act : %s on module : %i \n \n', opt_act, im)
            opt_act = str2num(opt_act(1));
            rewo = Rewarduhwi(im, state, opt_act, a_par,arr_par, dis_stp, v_pas, ...
            hw(im),cap_bus,lpass(state(1,im)),0, lapass(im), hw(implus), t_bo, t_al,hw(implusp), hwt);
            rewc = rewc + rewo;
            imr = opt_act + 1; %just a shift as the code below is written in terms of imr
            if nom
                imr = 1;
            end
            st_reach = state(1,im); %because we already have updated the stop in the begining of the loop
            b_st = [];
            ex_wt = 0; %extra wait the bus has to do because of already reched buses
            %following lines find the max waiting time among the buses already
            %at the stop where the vurrent bus (im) is reached.
            for i=1:n_b
                if state(1,i) == st_reach && atstop(i) == 1
                    b_st = [b_st i];
                    if t_nxt(i) > ex_wt
                        ex_wt = t_nxt(i);
                    end
                end
            end
            
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
                    tspent = max(t_al*pdsr , t_bo*pbsr) + fixdt + ex_wt + splittime;
                   
                    
                    %Updating t_nxt as the number of buses has changed
                    t_nxt_i = zeros(1,n_b);
                    for i = 1:im-1
                        t_nxt_i(i) = t_nxt(i);
                    end
                    for i = im+2:n_b
                        t_nxt_i(i) = t_nxt(i-1);
                    end
                    
    %                 if ex_wt ~=0
    %                     t_nxt_i(im) = ex_wt + 1; %adding 1 to ensure the module behind leaves later than the module in the front
    %                     atstop(im) = 1;
    %                 else
    %                     t_nxt_i(im) =  ct +  splittime;
    %                 end
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
                    
                    tspent = max(t_al*pdsr , t_bo*pbsr) + fixdt + ex_wt;
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
    %                 if gencount == 2
    %                     fprintf('load : %i, lapass(im) : %i, a_par : %f, pdsr : %i \n', state(2,im), lapass(im), a_par(state(1,im)), pdsr)
    %                     fprintf('rand eval : %i \n', binornd((state(2,im) - lapass(im)),a_par(state(1,im))))
    %                 end
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
    %                 if gencount == 85
    %                     disp('here')
    %                     disp(pdsr)
    %                     disp(pbsr)
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
                    %Not updating the headway as we are neglecting the time
                    %spent at the stop in the headway calculation
    
                     
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
            
    %         disp('atstop')
    %         disp(atstop)
    %         fprintf('rewsum 2, im : %i \n', im)
    %         disp('t_nxt2')
    %         disp(t_nxt)
            %opt_act1 = rewsum(1, state, t_nxt_sim, atstop, t_al, t_bo, arr_par, a_par, dis_stp, im, hway, v_bus, v_pas, lpass, lapass, cap_bus,hwt,ord);
    
            opt_act = rewsum_t(hor, state, t_nxt_sim, atstop, t_al, t_bo, arr_par, a_par, dis_stp, im, hway, v_bus, v_pas, lpass, lapass, cap_bus, hwt,ord, count, gencount, reachtime, leavetime, time,Hway, Pa_all_cum, fPas1, fPas2, jP1, jP2);
            %opt_act = rewsum_expt(hor, state, t_nxt_sim, atstop, t_al, t_bo, arr_par, a_par, dis_stp, im, hway, v_bus, v_pas, lpass, lapass, cap_bus,hwt,ord, count, gencount);
    
            opt_a = opt_act;
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
            rewo = Rewarduhwi(im, state, opt_act, a_par,arr_par, dis_stp, v_pas, ...
            hw(im),cap_bus,lpass(state(1,im)),0, lapass(im), hw(implus), t_bo, t_al,hw(implusp), hwt);
            rewc = rewc + rewo;
            if state(3,im) == 0 && state(3,implus) == 0 && imr == 4 %Join have positive reward and it is indeed possible to join
                r_jn = 1;
            else
                r_jn = -1; %This is dummy value if the statement in 'if' does not executte then join is anyway not feasible
            end
            leavetime(state(1,im)) = time;
            %r_jn  = -1;
            joincond =  r_jn > 0; % && ((atstop(implus) == 0 && state(1,implus) ~= state(1, im)) || state(1, implus) == state(1,im)); 
                if joincond
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
                        %updating hw
     
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
    %                     ord_i = ord;
    %                     ord(1) = ord_i(tot_mod);
    %                     ord(2:end) = ord_i(1:tot_mod-1);
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
        
                        %Updating headway
                        %No need to update headway here
    %                     hw_i = zeros(1,n_b);
    %                     for i=1:im-1
    %                         hw_i(i) = hw(i);
    %                     end
    %                     
    %                     hw_i(im) = hw(im) + tspent;
    %                     
    %                     for i = im+1:n_b
    %                         hw_i(i) = hw(i+1);
    %                     end
    %                     hw = hw_i;
        
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
        
        if time < 1300
            if imr == 1
                stcn = stcn + 1;
            elseif imr == 2
                skcn = skcn + 1;
            elseif imr == 3
                spcn  = spcn + 1;
            elseif imr == 4
                jncn = jncn + 1;
            elseif imr == 5
                nbcn = nbcn + 1;
            end
        end
        %      if gencount < 3
        %         fprintf('%s action on module %d \n', action, im)
        %         fprintf('action comb : %s \n', opt_a)
        %         fprintf('module %d spends time %f with extra time %f \n', im, tspent, ex_wt)
        %         fprintf('state at t = %d \n', time)
        %         fprintf('new added time is %d \n',M)
        %         %fprintf('pdsr : %i, pbsr : %i, ex_wt : %i, load : %i, a_par : %i, lapass(im) : %i \n', pdsr, pbsr, ex_wt, state(2,im), a_par(state(1,im)), lapass(im))
        %         disp(state)
        %         disp(ord)
        %         disp(atstop)
        %         disp('tnxt after')    
        %         disp(t_nxt)
        %         disp(lapass)
        %         disp(gencount)
        %     end
    % %     if gencount < 10
    %         disp(tw_pass)
    %         disp(w_pass)
    %     end
    
    
        
         %disp(gencount)   
            %Currently we are formulating the problem as if the stop, skip and
            %split will happen just at the stop but the join will happen after the
            %passengers have boarded and deboarded at the stop.
            
            %For now we are going forward with formulation such that the join
            %action will be considered only if the stop action has taken by the bus
            %in the front. 
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
        pWait = [pWait Pa_pre_cum(ceil(Time(i))) - Pb_cum(i)];
    end
    jAve = jAve + sjcount;
    skAve = skAve + skcount;
    stAve = stAve + stcount;
    nbAve = nbAve + snbcount;
    spAve = spAve + spcount;
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
    %disp(fPas2/parsep +1)
    costHor1(simn) = (w_wait*avg_wait1 + avg_inveh + w_walk*avg_walk1)/60;
end
fprintf('average travel time obtained : %f \n', mean(costHor1))
fprintf('split action average : %i \n', round(spAve/nsim))
fprintf('join action average : %i \n', round(jAve/nsim))
fprintf('stop action average : %i \n', round(stAve/nsim))
fprintf('skip action average : %i \n', round(skAve/nsim))
fprintf('next bs action average : %i \n', round(nbAve/nsim))