%Code to simulate the variable agent algorithm
rng(1)
n_s = 8; %Number of Stations
n_b = 2; %Number of bus, 1 bus = 2 modules
n_st = 3; %Number of states 
n_a = 4; % Number of actions bus can take
a_par = rand(1,n_s); %These are Ps values for the stops
arr_par = rand(1,n_s)/30; % Assuming on 30 passengers arrive in 30 minutes
atstop = zeros(1,n_b); %These are the flags which will be 1 if the corresponding bus is at stop and 0 it it is on the road

dis_stp = 300*(rand(1,n_s) + 1); %Distance between stops distributed between 300 - 600 meters
v_bus = 20*5/18; % Speed of bus 40 Km/h

cap_bus = 50;
unit_cap = cap_bus/2;
v_pas = 5.4*5/18; %Passenger speed in Km/h
t_bo = 5; %boarding time per passenger in seconds
t_al = 2; %Alighting time per passenger in seconds

gam_k = 2; gam_theta = 2; %parameters of Gamma rv

cumtime = zeros(1,n_b);

state0 = zeros(n_st,n_b);
state0(n_st,:) = ones(1,n_b);
state = state0;

%Busses will depart with headay of 90 seconds
%We will work on the basis of time a bus required to reach the next stop.
%It will we stored in the array n_xt. Note that this array will change as
%the code proceeds and its size will also change.
t_nxt = zeros(1,n_b); %this variable stores the time bus will require to reach the upcoming stop
for i=1:n_b  %At start each bus leaves with an interval of 90 seconds
    t_nxt(i) = 90*(i-1) + dis_stp(1)/v_bus;
end

l_action = zeros(n_st,n_bus); %last action taken by th bus

count = 1; %This will account the number of rounds
%consec_t = dis_stp/v_bus;  %We are assuming that, in the first round the passengers begin to arrive at the stop after the bus leaves the previous stop 
while true
    [M,im] = min(t_nxt);  %This step finds which bus reaches/leaves the next stop first
    
    if atstop(im) == 0
        if state(1,im) == n_s %Because our stops are circular n_s->1
            state(1,im) = 1;
        else
            state(1,im) = state(1,im) + 1; %state(1,im) is storing the last bus-stop number the bus has visited
        end
        
        for i=1:n_b
            if i == im
                t_nxt(i) = t_nxt(i);  %nxt for im will be calculated as the stop time at the stop
            else
                t_nxt(i) = t_nxt(i) - t_nxt(im);
            end       
        end
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
        hw = headway(state, t_nxt, atstop, dis_stp, v_bus, n_b,count);

        %The current headway computation does not consider the time spent at
        %the bus-stop. Will improve this in future
    
        %Now computing rewards
        %Need to calculate  lapass, lpass
        %lapass
        if l_action(im) ==1
            lapass = state(2,im)*a_par(im);  % this is a deterministic quantity
        end
    
        %acap: accomodating capacity
        acap = cap_bus - state(2,im);
    
        %lpass  for now ignoring the left over passengers
    %     if l_action(im+1) ==1
    %         lpass = arr_par(im);
    %     end
        %Computing rewards
        lpass = 0;

        r_st = Reward(im, state(:,im), 0, a_par,arr_par, dis_stp, v_pas, ...
        hw(im),cap_bus,lpass,l_action(im), lapass, hw(im-1), t_bo, t_al);
    
        r_sk = Reward(im, state(:,im), 1, a_par,arr_par, dis_stp, v_pas, ...
        hw(im),cap_bus,lpass,l_action(im), lapass, hw(im-1), t_bo, t_al);
    
        r_sp = Reward(im, state(:,im), 2, a_par,arr_par, dis_stp, v_pas, ...
        hw(im),cap_bus,lpass,l_action(im), lapass, hw(im-1), t_bo, t_al);
    
        if state(3,im) == 1 %If the bus is joined then the split action will always have highest reward as we are using the AVERAGE in rewards
            n_b = n_b +1; %The number of modules increase due to splitting
            %compute th number of passengers getting down at the current stop
            pds_c = a_par(state(1,im))*state(2,im);
            %Assuming that only the passengers getting down at the current stop
            %are in the rear module
            load_r = pds_c;
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
            %In the code we have distributed the load already so summing
            %loads of both rear and the front module in the above line

            state(2,im) = state(2,im) - pdsr; %This will be zero because the way we have distributed the passengers

            %In Zaid's work he has modelled the split such as the rear module
            %stops at the stop and deboards the passengers but don't boards any
            %passenger so that the rear module can catch up with the front one.
            %Right now I am allowing the passengers to board to make the
            %solution more general
            pbsr = max(poissrnd(arr_par(im)*hw(im)), cap_bus - state(2,im));
            state(2,im) = state(2,im) + pbsr;
            tspent = t_al*pdsr + t_bo*pbsr;

            %Updating t_nxt as the number of buses has changed
            t_nxt_i = zeros(n_b);
            for i = 1:im-1
                t_nxt_i(i) = t_nxt(i);
            end

            t_nxt_i(im) = dis_stp(state(1,im))/v_bus; %As the leading module does not stop.
            t_nxt_i(im + 1) = tspent; %time to leave the stop
            
            for i = im+2:n_b
                t_nxt_i(i) = t_nxt(i);
            end
            t_nxt = t_nxt_i;

            %Updating headway. Here doing manually instead of calling the
            %function as it is more accurate
            hw_i = zeros(n_b);
            for i = 1:im
                hw_i(i) = hw(i);
            end
            hw_i(im+1) = tspent;
            for i = im+1:n_b
                hw_i(i) = hw(i);
            end
            hw = hw_i;
        end
    
        if state(3,im) == 0  %if the bus is splitted already then the competition will be in stop and skip
            if r_st > r_sk
                pdsr = binornd(state(2,im),a_par(state(1,im)));   
                state(2,im) = state(2,im) - pdsr;
                pbsr = max(poissrnd(arr_par(im)*hw(im)), bus_cap - state(2,im));
                state(2,im) = state(2,im) + pbsr;
                tspent = t_al*pdsr + t_bo*pbsr;
                %t_nxt, hw, atstop, 
               
                 
            else  %If skip action is taken then the time spent at the stop will be zero
                tspent = 0;
                t_nxt(im) = dis_stp(state(1,im))/v_bus;
            end
        end 
    end 
    
    if atstop(im) == 1
        r_jn = Reward(im, state(:,im), 3, a_par,arr_par, dis_stp, v_pas, ...
        hw(im),cap_bus,lpass,l_action(im), lapass, hw(im-1), t_bo, t_al);
        
        if r_jn > 0
            n_b = n_b -1;

            %Updating state
            state_i = zeros(n_st,n_b);
            for i=1:im-1
                state_i(:,i) = state(:,i);
            end
            state_i(im) = [state(1,im);state(2,im)+state(2,im+1);1];
            for i=im+1:n_b
                state_i(:,i) = state(:,i+1);
            end
            
            %Updating t_nxt
            t_nxt_i = zeros(1,n_b);
            for i = 1:im-1
                t_nxt_i(i) = t_nxt(i)

    end
        
    
        %Currently we are formulating the problem as if the stop, skip and
        %split will happen just at the stop but the join will happen after the
        %passengers have boarded and deboarded at the stop.
        
        %For now we are going forward with formulation such that the join
        %action will be considered only if the stop action has taken by the bus
        %in the front.
    


         
end    