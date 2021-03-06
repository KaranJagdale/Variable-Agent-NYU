%Code to simulate the variable agent algorithm
rng(1)
%This code is being run
%dbstop if naninf
n_s = 8; %Number of Stations
n_b = 2; %Number of bus, 1 bus = 2 modules
n_st = 3; %Number of states 
n_a = 4; % Number of actions bus can take
a_par = rand(1,n_s); %These are Ps values for the stops
arr_par = rand(1,n_s)/30*3; % Assuming on 90 passengers arrive in 30 minutes
atstop = zeros(1,n_b); %These are the flags which will be 1 if the corresponding bus is at stop and 0 it it is on the road
fixdt = 30; %fixed time lost per stop
dis_stp = 300*(rand(1,n_s) + 1); %Distance between stops distributed between 300 - 600 meters
v_bus = 20*5/18; % Speed of bus 20 Km/h

cap_bus = 50;
unit_cap = cap_bus/2;
v_pas = 5.4*5/18; %Passenger speed in Km/h
t_bo = 5; %boarding time per passenger in seconds
t_al = 2; %Alighting time per passenger in seconds

gam_k = 2; gam_theta = 2; %parameters of Gamma rv

cumtime = zeros(1,n_b);

state0 = zeros(n_st,n_b);
state0(n_st,:) = ones(1,n_b);
state0(1,:) = ones(1,n_b);
state = state0;

%Busses will depart with headay of 90 seconds
%We will work on the basis of time a bus required to reach the next stop.
%It will we stored in the array n_xt. Note that this array will change as
%the code proceeds and its size will also change.
t_nxt = zeros(1,n_b); %this variable stores the time bus will require to reach the upcoming stop
for i=1:n_b  %At start each bus leaves with an interval of 90 seconds
    t_nxt(i) = 500*(i-1) + dis_stp(1)/v_bus;
end
%disp(t_nxt)
l_action = zeros(n_st,n_b); %last action taken by th bus

count = 1; %This will account the number of rounds
%consec_t = dis_stp/v_bus;  %We are assuming that, in the first round the passengers begin to arrive at the stop after the bus leaves the previous stop 
gencount = 0;
skcount = 0;
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
Pd_cum = 0;
Pa_cum = 0;
Pb_cum = 0;
Pw_cum = 0;
pw_cum = 0;
pd_cum = 0; %cumulative deboarding passengers
pb_cum = 0; %cumulative boarding passengers
pa_cum = 0; %cumulative arriving passengers
ord = [1 2 3 4]; %will keep track of the order of the buses. This is basically 
%the mapping between the state index and the bus number
while true
    %fprintf('n_b = %f \n',n_b)
    %fprintf('size of t_nxt = %f \n',size(t_nxt,2))
%     if gencount < 100
%         disp(state(2,:))
%     end
    gencount = gencount + 1;
    %disp(gencount)
%     if gencount == 2
% 
%         disp(t_nxt)
%     end
    hw = headwaycu(state, t_nxt, atstop, dis_stp, v_bus, n_b,count);

    for i=hw
            if i<0
                disp('danger')
                disp(gencount)
                disp(hw)
                disp(state)
                disp(t_nxt)
                disp(atstop)
                break
            end
    end

    [M,im] = min(t_nxt);  %This step finds which bus reaches/leaves the next stop first
    T = [T M];
    
    implus = iplus(im,n_b);
    for i=1:n_b
            if i ~= im
                t_nxt(i) = t_nxt(i) - t_nxt(im);  %nxt for im will be calculated as the stop time at the stop
            end       
    end
    

    if atstop(im) == 0
        if state(1,im) == n_s %Because our stops are circular n_s->1
            state(1,im) = 1;
        else
            state(1,im) = state(1,im) + 1; %state(1,im) is storing the last bus-stop number the bus has visited
        end
        
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
        %lapass
        if l_action(im) ==1
            lapass = state(2,im)*a_par(state(1,im));  % this is a deterministic quantity
        else
            lapass = 0;
        end
        
        %acap: accomodating capacity
        acap = cap_bus - state(2,im);
    
        %lpass: for now ignoring the left over passengers
    %     if l_action(im+1) ==1
    %         lpass = arr_par(im);
    %     end
        %Computing rewards
        lpass = 0;

        
        %implus is the index of the bus behind im
        r_st = Reward(im, state, 0, a_par,arr_par, dis_stp, v_pas, ...
        hw,unit_cap,lpass,l_action(im), lapass, hw(implus), t_bo, t_al);
    
        r_sk = Reward(im, state, 1, a_par,arr_par, dis_stp, v_pas, ...
        hw,unit_cap,lpass,l_action(im), lapass, hw(implus), t_bo, t_al);
    
        r_sp = Reward(im, state, 2, a_par,arr_par, dis_stp, v_pas, ...
        hw,unit_cap,lpass,l_action(im), lapass, hw(implus), t_bo, t_al);
        
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
            n_b = n_b +1; %The number of modules increase due to splitting
            %compute th number of passengers getting down at the current stop
            pds_c = a_par(state(1,im))*state(2,im);
            %Assuming that only the passengers getting down at the current stop
            %are in the rear module
            
            load_r = floor(pds_c);
            if state(2,im) - load_r > unit_cap
                load_r = state(2,im) - unit_cap;
            end

            load_f = ceil(state(2,im) - load_r);
%             if gencount == 400
%                 fprintf('pds_c = %f \n',pds_c)
%                 fprintf('load_f = %f \n',load_f)
%             end
%             if gencount == 1
%                 disp('good')
%             end
            
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
            atstop_i(im + 1) = 1; %as rear module has stopped as the stop
            for i = im + 2 : n_b
                atstop_i(i) = atstop(i-1);
            end
            atstop = atstop_i;

            
            
            %At this point we have distributed the load and updated the states
            %due to increase in the number of busses
    
            %Now the rear bus will stop at the stop and the front bus will
            %leave
            
            pdsr = binornd(state(2,im)+state(2,im+1),a_par(state(1,im))); %all the passengers in the joint bus will be moved to the rear bus
            if pdsr > state(2,im+1)  %pdsr is random so can be more than the number of passengers in the rear bus
                pdsr = state(2,im+1);
            end
            pd_cum = pd_cum + pdsr;
            
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
            bpass = arr_par(state(1,im+1))*hw(im);
            pa = poissrnd(bpass);
            pa_cum = pa_cum + pa;
            
            pbsr = min(pa, unit_cap - state(2,im+1));
            pb_cum = pb_cum + pbsr;
            
%             if gencount == 498
%                 disp(state(2,im+1))
%                 disp(pbsr)
%                 fprintf('bpass  = %f \n',bpass)
%                 fprintf('headway = %f \n',hw(im))
%             end
            state(2,im+1) = state(2,im+1) + pbsr;
            tspent = t_al*pdsr + t_bo*pbsr + fixdt + ex_wt;
           
            
            %Updating t_nxt as the number of buses has changed
            t_nxt_i = zeros(1,n_b);
            for i = 1:im-1
                t_nxt_i(i) = t_nxt(i);
            end

            t_nxt_i(im) = dis_stp(state(1,im))/v_bus + ex_wt; %As the leading module does not stop.
            t_nxt_i(im + 1) = tspent; %time to leave the stop
            
            for i = im+2:n_b
                t_nxt_i(i) = t_nxt(i-1);
            end
            t_nxt = t_nxt_i;

            %Updating headway. Here doing manually instead of calling the
            %function as it is more accurate
            hw_i = zeros(1,n_b);
            for i = 1:im
                hw_i(i) = hw(i);
            end
            hw_i(im+1) = tspent;
            for i = im+2:n_b
                hw_i(i) = hw(i-1);
            end
            hw = hw_i;
        
    
        else  %if the bus is splitted already then the competition will be in stop and skip
            if r_st > r_sk
                action = 'stop';
                pdsr = binornd(state(2,im),a_par(state(1,im)));
                pd_cum = pd_cum + pdsr;
                

                state(2,im) = state(2,im) - pdsr;
                pa = poissrnd(arr_par(state(1,im))*hw(im));
                pa_cum = pa_cum + pa;
                
                pbsr = min(pa, unit_cap - state(2,im));
                pb_cum = pb_cum + pbsr;
                
                state(2,im) = state(2,im) + pbsr;
                tspent = t_al*pdsr + t_bo*pbsr + fixdt + ex_wt;
                %Updating t_nxt, atstop, 
                atstop(im) = 1;
                t_nxt(im) = tspent;
                %Not updating the headway as we are neglecting the time
                %spent at the stop in the headway calculation

                 
            else  %If skip action is taken then the time spent at the stop will be zero
                pw = binornd(state(2,im),a_par(state(1,im))); %Number of passengers have to walk
                pw_cum = pw_cum + pw;
                
                skcount = skcount + 1;
                action = 'skip';
                tspent = 0 + ex_wt;
                t_nxt(im) = dis_stp(state(1,im))/v_bus + ex_wt;
            end
        end 
     
    
    else %If the bus is not at the stop we will evalueate join possibility
        
%         if gencount == 2
%             disp(r_jn)
%             disp(state(3,im))
%             disp(state(3,implus))
%         end
        if state(3,im) == 0 && state(3,implus) == 0 %Join have positive reward an it is indeed possible to join
            r_jn = Reward(im, state, 3, a_par,arr_par, dis_stp, v_pas, ...
                   hw,unit_cap,lpass,l_action(im), lapass, hw(implus), t_bo, t_al);
        else
            r_jn = -1; %This is dummy value if the statement in 'if' does not executte then join is anyway not feasible
        end
            if r_jn > 0 && state(3,im) == 0 && state(3,implus) == 0

                if gencount == 11
                    disp('passing join')
                    disp(n_b)
                end
                action = 'join';
                tspent = hw(implus);  %Front bus has to wait for this time
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
                    t_nxt_i(n_b) = dis_stp(state(1,n_b))/v_bus + tspent;
                    t_nxt = t_nxt_i;
                    %updating hw

                    hw_i = zeros(1,n_b);

                    for i=2:n_b
                        hw_i(i-1) = hw(i);
                    end

                    hw_i(n_b) = hw(n_b+1) + tspent;

                    hw = hw_i;

                    %updating atstop

                    atstop_i = zeros(1,n_b);

                    for i=2:n_b
                        atstop_i(i-1) = atstop(i);
                    end
                    atstop_i(n_b) = 0;
                    atstop = atstop_i;


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
                    
                    t_nxt_i(im) = dis_stp(state(1,im))/v_bus + tspent;
    
                    for i =im+1:n_b
                        t_nxt_i(i) = t_nxt(i+1);
                    end
                    t_nxt = t_nxt_i;
    
                    %Updating headway
                    hw_i = zeros(1,n_b);
                    for i=1:im-1
                        hw_i(i) = hw(i);
                    end
                    
                    hw_i(im) = hw(im) + tspent;
                    
                    for i = im+1:n_b
                        hw_i(i) = hw(i+1);
                    end
                    hw = hw_i;
    
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
                end

        %If join is not feasible  
            else
                action = 'nextbs';
                t_nxt(im) = dis_stp(state(1,im))/v_bus;
                atstop(im) = 0;
                tspent = dis_stp(state(1,im))/v_bus;  %this tspent can be the t_nxt two lines above

            end
        
           
    end
    
    time  = time + M;
    %Now we construct arrays that consist the bus locations at different
    %time instants
    Pa_cum = [Pa_cum pa_cum];
    Pd_cum = [Pd_cum pd_cum];
    Pb_cum = [Pb_cum pb_cum];
    Pw_cum = [Pw_cum pw_cum];
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
  
    if gencount == 1000
        break
    end
    count = count + 1;

%     if gencount == 713 || gencount == 714        
%         fprintf('%s action on module %d \n', action, im)
%         fprintf('module %d spends time %f with extra time %f \n', im, tspent, ex_wt)
%         fprintf('state at t = %d \n', time)
%         fprintf('new added time is %d \n',M)
%         disp(state)
%         disp(ord)
%         disp(atstop)
%         disp('tnxt after')
%         disp(t_nxt)
%         disp(gencount)
%     end


    
     %disp(gencount)   
        %Currently we are formulating the problem as if the stop, skip and
        %split will happen just at the stop but the join will happen after the
        %passengers have boarded and deboarded at the stop.
        
        %For now we are going forward with formulation such that the join
        %action will be considered only if the stop action has taken by the bus
        %in the front.    
end 
%now fining the area under the curves



figure(1)
plot(Time1, loc1)
xlabel('time (seconds)')
ylabel('Bus Stop')
hold on 
%plot(Time, locf1)
figure(2)
plot(Time2,loc2)
xlabel('time (seconds)')
ylabel('Bus Stop')
figure(3)
plot(Time3,loc3)
xlabel('time (seconds)')
ylabel('Bus Stop')
figure(4)
plot(Time4,loc4)
xlabel('time (seconds)')
ylabel('Bus Stop')
pltsz = 100:175;
figure(5)
plot(Time1(pltsz), loc1(pltsz))
hold on
plot(Time2(pltsz),loc2(pltsz))
plot(Time3(pltsz),loc3(pltsz))
plot(Time4(pltsz),loc4(pltsz))
legend('m1','m2','m3','m4')
xlabel('time (seconds)')
ylabel('Bus Stop')
figure(6)
plot(Time, Pb_cum)
hold on
plot(Time, Pa_cum)
plot(Time, Pd_cum)
plot(Time, Pw_cum)
legend('boarding', 'awaiting','de-boarding','walking')
% hold on
% plot(Time, State2(1,:))
%plot()