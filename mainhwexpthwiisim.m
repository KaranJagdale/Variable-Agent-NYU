%Before running the code check state0, n_b, Cost array name, bet, alpha,
%hw_t

nsim = 100;
Costsplit = zeros(1,nsim);

for simno=1:nsim
    %%Code to simulate the variable agent algorithm
    rng(simno)
    disp(simno)
    %This code is using the most recent headway function that is headwaycuu
    %this code also accounts the leftover passengers if the bus skips the stop
    %but at this point of time we are assming that the bus not skipping stops
    %successively
    %dbstop if naninf
    simtime = 2; %simulation time in hours
    n_s = 8; %Number of Stations
    n_b = 2; %Number of bus, 1 bus = 2 modules
    n_st = 3; %Number of states 
    n_a = 4; % Number of actions bus can take
    %a_par = rand(1,n_s); %These are Ps values for the stops
    almu = 2/n_s;
    a_par = normrnd(almu, almu/10, 1, n_s); 
    %arr_par = rand(1,n_s)/60*2; % Assuming on 90 passengers arrive in 30 minutes
    arrmu = 0.015/3; arrsigma = arrmu/10;
    arr_par = normrnd(arrmu, arrsigma, 1,n_s);
    atstop = ones(1,n_b); %These are the flags which will be 1
    % if the corresponding bus is at stop and 0 it it is on the road
    fixdt = 30; %fixed time lost per stop
    dis_stp = 120*(rand(1,n_s) + 1); %Distance between stops distributed between 300 - 600 meters
    v_bus = 20*5/18; % Speed of bus 20 Km/h
    w_wait = 2.1; w_walk = 2.2;% weights of walk time and wait time in final cost
    lapass = zeros(1,n_b);
    lpass = zeros(1,n_s);
    w_pass = zeros(4,n_s); %to store the number of walking passengers to a particular stop
    %at max there can be 4 types of passengers walking towards a particular
    %stop
    tw_pass = zeros(4,n_s); %time to reach the next stop
    cap_bus = 50;
    unit_cap = cap_bus/2;
    v_pas = 5.4*5/18; %Passenger speed in Km/h
    t_bo = 5; %boarding time per passenger in seconds
    t_al = 2; %Alighting time per passenger in seconds
    scamcount = 0;
    gam_k = 2; gam_theta = 2; %parameters of Gamma rv
    gamma = 1.5;
    bet = 0.6;
    ifsp = 1; %ifsp = 1000 if doing normal
    hwt_i = (sum(dis_stp)/v_bus + n_s*fixdt)/(n_b - (t_al + t_bo)*arrmu*n_s);
    hwt = ifsp*hwt_i;
    cumtime = zeros(1,n_b);
    splittime = 1; %splittime being 1 second -- Introducing this constant to avoid bus overtake
    stload = floor(arrmu*n_s/2*hwt_i); 
    state0 = zeros(n_st,n_b);
    state0(n_st,:) = ones(1,n_b);
    state0(1,:) = ones(1,n_b);
    state0(2,:) = stload;
    state = state0;
    
   
    t_nxt = zeros(1,n_b); %this variable stores the time bus will require to reach the upcoming stop
    for i=1:n_b  %At start each bus leaves with an interval of hw_i seconds
        t_nxt(i) = hwt_i*(i-1) + dis_stp(1)/v_bus;
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
    %the mapping between the state index and the bus number
    %Creating arriving passenger array 
    % pa_pre = poissrnd(arrmu*n_s, 1, simtime*3600);
    Pa_pre_cum = zeros(1,simtime*3600+1);
    % pa_pre_cum = stload*n_b;
    Time_pre = 0:simtime*3600;
  
    %generating arrivals for each stop
    tflst = sum(dis_stp) + (t_al + t_bo)*hwt_i*sum(arr_par) + fixdt*n_s;
    tfstops = zeros(1,n_s);
    for i = 1:n_s -1
        tfstops(i+1) = tfstops(i) + dis_stp(i) + (t_al + t_bo)*hwt_i*arr_par(i) + fixdt;
    end
    Pa_all_cum = zeros(n_s, simtime*3600+1);
    nwpasscount = 0;
    for i = 1:n_s
        %passnger coming to ith stop in one target headway
        passhwt = poissrnd(arr_par(i)*hwt_i);
        for k = 1:n_s-1
            for j=floor(tfstops(k))+1:floor(tfstops(k+1))
                if Pa_all_cum(i,j) >= passhwt && i > k
                    newpass = 0;
                    nwpasscount = nwpasscount + 1;
                else
                    newpass = poissrnd(arr_par(i));
                end
                Pa_all_cum(i,j+1) = Pa_all_cum(i,j) + newpass;
            end
            for j = floor(tfstops(k+1))+1: simtime*3600
                Pa_all_cum(i,j+1) = Pa_all_cum(i,j) + poissrnd(arr_par(i));
            end
        end
    end
    padec = zeros(n_s, simtime*3600+1);
    for i = 1: 3600*simtime + 1
        Pa_pre_cum(i) = sum(Pa_all_cum(:,i)) + stload*n_b;
    end
    
    while true
        hwt_i = (sum(dis_stp)/v_bus + n_s*fixdt)/(n_b - (t_al + t_bo)*arrmu*n_s);
        hwt = ifsp*hwt_i; 
        if state(1,1) == 8 %count encodes whether the leading module has covered a circle or not
            count = 0;
        end
   
        gencount = gencount + 1;
  
        hw = headwaycuu(state, t_nxt, atstop, dis_stp, v_bus, n_b,count);
    
    
    
        [M,im] = min(t_nxt);  %This step finds which bus reaches/leaves the next stop first
        time = time + M;
        T = [T M];
        %finding other indices with the minimum value
        imin = find(t_nxt == M);
        im = min(imin);
        
        implus = iplus(im,n_b);
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
        
    
        if atstop(im) == 0
            
           % hway = hw(im);
    
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
            bstop = state(1,im);
            %disp(t_nxt)
            
            %nxt is the time to reach a stop when bus is not at the stop and it is
            %equal to the time to leave the stop when the bus  is as the stop. With
            %this definitions will have to update the way headway is calculated.
            
            %Here we will assume that the last bus starts from the first stop
            %before the first bus completes the circle. Which kind of makes sense
            %as otherwise it will mean that the number of buses is more than
            %required.
            %While computing the headway we are assuming the 
            %Computing headway is divided in four different cases depending on
            %whether the two busses are either on the road or at the stop
            
    
            
            %The current headway computation does not consider the time spent at
            %the bus-stop. Will improve this in future
        
            %Now computing rewards
            %Need to calculate  lapass, lpass
    
    
            
            %acap: accomodating capacity
            acap = cap_bus - state(2,im);
        
            %lpass: for now ignoring the left over passengers
        %     if l_action(im+1) ==1
        %         lpass = arr_par(im);
        %     end
            %Computing rewards
            %lpass = 0;
    
            
            %implus is the index of the bus behind im
            r_st = Rewardhwi(im, state, 0,hway,n_s, gamma, bet, hwt);
        
            r_sk = Rewardhwi(im, state, 1,hway,n_s, gamma, bet, hwt);
        
            r_sp = Rewardhwi(im, state, 2,hway,n_s, gamma, bet, hwt);
           % r_st = 1; r_sk = 0;
    %         fprintf('Stop rew %f \n', r_st)
    %         fprintf('Split rew %f \n', r_sp)
            
    %         if r_st ~= -r_sk && gencount < 13
    %             fprintf('Issue at count = %f \n', gencount)
    %             break
    %         end
            %Finding if some other bus/module is already at the stop. This
            %is to avoid overtaking.
            st_reach = state(1,im); %because we already have updated the stop in the begining of the loop
            b_st = [];
            ex_wt = 0; %extra wait the bus has to do because of already reched buses
            for i=1:n_b
                if state(1,i) == st_reach && atstop(i) == 1
                    b_st = [b_st i];
                    ex_wt = ex_wt + t_nxt(i);
                end
            end
    %         if gencount == 4 || gencount == 3
    %             disp('this is imp')
    %             disp(state)
    %             disp(b_st)
    %             fprintf('ex_wt %f', ex_wt)
    %             fprintf('st_reach %f',st_reach)
    %         end
        
            if state(3,im) == 1 %If the bus is joined then the split action will always have highest reward as we are using the AVERAGE in rewards
                
                action  = 'split';
                spcount = spcount +1;
                n_b = n_b +1; %The number of modules increase due to splitting
                %compute th number of passengers getting down at the current stop
                pds_c = binornd((state(2,im)-lapass(im)),a_par(state(1,im))) + lapass(im);
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
                %At this point we have distributed the load and updated the states
                %due to increase in the number of busses
        
                %Now the rear bus will stop at the stop and the front bus will
                %leave
                
                %pdsr = binornd(state(2,im)+state(2,im+1),a_par(state(1,im))); %all the passengers in the joint bus will be moved to the rear bus
    %             if pdsr > state(2,im+1)  %pdsr is random so can be more than the number of passengers in the rear bus
    %                 pdsr = state(2,im+1);
    %             end
                
                
    %             if gencount == 1
    %                 disp(state)
    %                 disp(pdsr)
    %             end
                %In the code we have distributed the load already so summing
                %loads of both rear and the front module in the above line
    %             if gencount == 2
    %                 disp(state(2,im))
    %                 disp(pdsr)
    %             end
                state(2,im+1) = state(2,im+1) - pdsr; %This will be zero because the way we have distributed the passengers
                
                %In Zaid's work he has modelled the split such as the rear module
                %stops at the stop and deboards the passengers but don't boards any
                %passenger so that the rear module can catch up with the front one.
                %Right now I am allowing the passengers to board to make the
                %solution more general
    
                %in hw we have used im although technically we should pass im+1
                %but, we haven't updated hw yet and headway for both im and
                %im+1 will be same as they were joined before
                
                
                
                %bpass = arr_par(state(1,im+1))*hway;
    %             if gencount < 10
    %                 disp(state)
    %             end
                
    %             f_ind = ceil(reachtime(state(1,im))); 
    %             s_ind = floor(time);
    %             if f_ind == 0
    %                 f_ind = f_ind + 1;
    %             end
    %             if 
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
    %             if gencount == 498
    %                 disp(state(2,im+1))
    %                 disp(pbsr)
    %                 fprintf('bpass  = %f \n',bpass)
    %                 fprintf('headway = %f \n',hw(im))
    %             end
                state(2,im+1) = state(2,im+1) + pbsr;
                tspent = t_al*pdsr + t_bo*pbsr + fixdt + ex_wt + splittime;
               
                
                %Updating t_nxt as the number of buses has changed
                t_nxt_i = zeros(1,n_b);
                for i = 1:im-1
                    t_nxt_i(i) = t_nxt(i);
                end
    
                t_nxt_i(im) = dis_stp(state(1,im))/v_bus + ex_wt + splittime; %As the leading module does not stop.
                t_nxt_i(im + 1) = tspent; %time to leave the stop
                
                for i = im+2:n_b
                    t_nxt_i(i) = t_nxt(i-1);
                end
                t_nxt = t_nxt_i;
    
    
            
        
            else  %if the bus is splitted already then the competition will be in stop and skip
                if r_st > r_sk
                    action = 'stop';
                    stcount = stcount + 1;
                    pdsr = lapass(im) + binornd((state(2,im) - lapass(im)),a_par(state(1,im)));
    %                 if gencount == 85
    %                     disp('hehe')
    %                     disp(lapass(im))
    %                     disp(state(2,im))
    %                     disp(binornd((state(2,im) - lapass(im)),a_par(state(1,im))))
    %                 end
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
                    
                    tspent = t_al*pdsr + t_bo*pbsr + fixdt + ex_wt;
                    %Updating t_nxt, atstop, 
                    atstop(im) = 1;
                    t_nxt(im) = tspent;
                    %Not updating the headway as we are neglecting the time
                    %spent at the stop in the headway calculation
    
                     
                else  %If skip action is taken then the time spent at the stop will be zero
                    pw = binornd((state(2,im) - lapass(im)),a_par(state(1,im)));% + lapass(im); %Number of passengers have to walk. Can modify this line to incorporate multiple skip leftover passengers
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
                    tspent = 0 + ex_wt;
                    if ex_wt ~=0
                        atstop(im) = 1;
                        t_nxt(im) = ex_wt;
                    else
                        t_nxt(im) = dis_stp(state(1,im))/v_bus + ex_wt;
                    end
                end
                
            end 
            reachtime(bstop) = time;
         
        
        else %If the bus is not at the stop we will evalueate join possibility
            
            if im == 1 && count == 1
                hway = hwt;
            else
                hway = time - leavetime(state(1,im));
            end
            hway = hw(im);
    %         if gencount == 2
    %             disp(r_jn)
    %             disp(state(3,im))
    %             disp(state(3,implus))
    %         end
            if state(3,im) == 0 && state(3,implus) == 0 %Join have positive reward and it is indeed possible to join
                r_jn = Rewardhwi(im, state, 3,hway,n_s, gamma, bet, hwt);
            else
                r_jn = -1; %This is dummy value if the statement in 'if' does not executte then join is anyway not feasible
            end
            leavetime(state(1,im)) = time;
            %r_jn  = -1;
            joincond =  r_jn > 0  && ((atstop(implus) == 0 && state(1,implus) ~= state(1, im)) || state(1, implus) == state(1,im)); 
                if joincond
                    action = 'join';
                    sjcount = sjcount + 1;
                    if atstop(implus) == 0
                        tspent = hw(implus);  %Front bus has to wait for this time
                    else
                        tspent = t_nxt(implus);
                    end
                    load_rj = state(2,implus); % this value will requre to caculate lapass
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
                        t_nxt_i(n_b) = t_nxt(n_b + 1) + tspent;
                        %t_nxt_i(n_b) =  tspent;
                        t_nxt = t_nxt_i;
                        %updating hw
                        %no need to update headway here
    
    %                     hw_i = zeros(1,n_b);
    % 
    %                     for i=2:n_b
    %                         hw_i(i-1) = hw(i);
    %                     end
    % 
    %                     hw_i(n_b) = hw(n_b+1) + tspent;
    % 
    %                     hw = hw_i;
    
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
                        pw = binornd((load_rj-lapass(1)),a_par(state(1,n_b)));
                        lapass_i = zeros(1,n_b);
                        for i =2:n_b
                            lapass_i(i-1) = lapass(i);
                        end
                        lapass_i(n_b) = lapass(1)+ lapass(n_b+1) + pw;
                        lapass = lapass_i;
    
                        %In this case the bus 2 will correspond to state(:,1)
                        ord_i = ord;
                        ord(1) = ord_i(4);
                        ord(2:4) = ord_i(1:3);
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
                        
                        t_nxt_i(im) = t_nxt(im) + tspent;
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
                        pw = binornd((load_rj-lapass(implus)),a_par(state(1,im)));
                        lapass_i = zeros(1,n_b);
                        for i=1:im-1
                            lapass_i(i) = lapass(i);
                        end
                        lapass_i(im) = lapass(im) + lapass(im+1) + pw;
                        for i=im+1:n_b
                            lapass_i(i) = lapass(i+1);
                        end
                        lapass = lapass_i;
                    end              
            %If join is not feasible  
                else
                    snbcount = snbcount + 1;
                    action = 'nextbs';
                    t_nxt(im) = dis_stp(state(1,im))/v_bus;
                    atstop(im) = 0;
                    tspent = dis_stp(state(1,im))/v_bus;  %this tspent can be the t_nxt two lines above
                end          
        end
        
        if time < simtime*3600
        
            if Pa_pre_cum(ceil(time)) < pb_cum 
                disp('khatra')
                disp(pb_cum)
                disp(Pa_pre_cum(int16(time)))
                %break
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
        if n_b == 2
            bus_loc = [state(1,1) state(1,1) state(1,2) state(1,2)];
            atstopf = [atstop(1) atstop(1) atstop(2) atstop(2)];
        elseif n_b == 3
            if state(3,1) == 1
                bus_loc = [state(1,1) state(1,1) state(1,2) state(1,3)];
                atstopf = [atstop(1) atstop(1) atstop(2) atstop(3)];
            elseif state(3,2) == 1
                bus_loc = [state(1,1) state(1,2) state(1,2) state(1,3)];
                atstopf = [atstop(1) atstop(2) atstop(2) atstop(3)];
            elseif state(3,3) == 1
                bus_loc = [state(1,1) state(1,2) state(1,3) state(1,3)];
                atstopf = [atstop(1) atstop(2) atstop(3) atstop(3)];
            end
        elseif n_b == 4
            bus_loc = state(1,:);
            atstopf = atstop;
        end
        locf1 = [locf1 bus_loc(ord(1))];
        Time = [Time time];
        if atstopf(ord(1)) == 1 || bus_loc(ord(1)) ~= loc1(size(loc1,2))
            loc1 = [loc1 bus_loc(ord(1))];
            Time1 = [Time1 time];
        end
        if atstopf(ord(2)) == 1 || bus_loc(ord(2)) ~= loc2(size(loc2,2))
            loc2 = [loc2 bus_loc(ord(2))];
            Time2 = [Time2 time];
        end
    
        if atstopf(ord(3)) == 1 || bus_loc(ord(3)) ~= loc3(size(loc3,2))
            loc3 = [loc3 bus_loc(ord(3))];
            Time3 = [Time3 time];
        end
        if atstopf(ord(4)) == 1 || bus_loc(ord(4)) ~= loc4(size(loc4,2))
            loc4 = [loc4 bus_loc(ord(4))];
            Time4 = [Time4 time];
        end
    
        %Time = [Time time];
    % 
    %     if time > 1450 && time < 2334
    %         fprintf('state at t = %d \n', gencount)
    %         disp(state)
    %     end
    %     if state(1,1) == state(1,2) && t_nxt(1)>t_nxt(2) && atstop(1) == 0 && atstop(2) == 0
    %         disp('we found')
    %         disp(gencount)
    %         
    %         break
    %     end
      
        if time > simtime*60*60
            break
        end
        %count = count + 1;
    
    %     if gencount >30 && gencount < 40
    %         fprintf('%s action on module %d \n', action, im)
    %         fprintf('module %d spends time %f with extra time %f \n', im, tspent, ex_wt)
    %         fprintf('state at t = %d \n', time)
    %         fprintf('new added time is %d \n',M)
    %         disp(state)
    %         disp(ord)
    %         disp(atstop)
    %         disp('tnxt after')
    %         disp(lapass)
    %         disp(t_nxt)
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
    %now finding the area under the curves
    Pb_cum_int = trapz(Time, Pb_cum);
    Pd_cum_int = trapz(Time, Pd_cum);
    Pa_cum_int = trapz(Time, Pa_cum);
    Pw_cum_int = trapz(Time, Pw_cum);
    Pe_cum_int = trapz(Time, Pe_cum);
    Pa_pre_cum_int = trapz(1, Pa_pre_cum);
    
    avg_inveh = (Pb_cum_int - Pd_cum_int)/Pb_cum(size(Pb_cum,2));
    avg_wait = (Pa_cum_int - Pb_cum_int)/Pa_cum(size(Pa_cum,2));
    avg_walk = Pw_cum_int/Pd_cum(size(Pd_cum,2)); %this formulation is wrong
    avg_walk1 = (Pd_cum_int - Pe_cum_int)/Pd_cum(size(Pd_cum,2));
    avg_wait1 = (Pa_pre_cum_int - Pb_cum_int)/Pa_pre_cum(size(Pa_pre_cum,2));
%     fprintf('Average in-vehicle time (m) : %f \n', avg_inveh/60)
%     fprintf('Average waiting time (m) : %f \n', avg_wait1/60)
%     fprintf('Average walk time (m) : %f \n', avg_walk1/60)
%     fprintf('Policy cost (m) : %f \n', (w_wait*avg_wait1 + avg_inveh + w_walk*avg_walk1)/60)
%     fprintf('final time : %f hr \n', time/3600)
%     fprintf('join count : %i \n', sjcount)
%     fprintf('skip count : %i \n', skcount)
%     fprintf('stop count : %i \n', stcount)
    
    
    fprintf('Policy cost (m) : %f \n', (w_wait*avg_wait1 + avg_inveh + w_walk*avg_walk1)/60)
    Costsplit(simno) = (w_wait*avg_wait1 + avg_inveh + w_walk*avg_walk1)/60;
    
end
