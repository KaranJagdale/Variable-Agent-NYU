%Code to simulate the variable agent algorithm
rng(1)
close all
%This code is being run
%dbstop if naninf
n_s = 8; %Number of Stations
n_b = 4; %all module separated
n_st = 3; %Number of states 
n_a = 4; % Number of actions bus can take
a_par = rand(1,n_s); %These are Ps values for the stops
arr_par = rand(1,n_s)/30*3; % Assuming on 90 passengers arrive in 30 minutes
atstop = ones(1,n_b); %These are the flags which will be 1 if the corresponding bus is at stop and 0 it it is on the road
fixdt = 30; %fixed time lost per stop
dis_stp = 300*(rand(1,n_s) + 1); %Distance between stops distributed between 300 - 600 meters
v_bus = 20*5/18; % Speed of bus 20 Km/h
w_wait = 2.1; w_walk = 2.2;% weights of walk time and wait time in final cost
gamma = 1.5; %threshold headway multiplier
headwayt = 120; %threshold headway for 120 the graph looks good!
cap_bus = 50;
unit_cap = cap_bus/2;
v_pas = 5.4*5/18; %Passenger speed in Km/h
t_bo = 5; %boarding time per passenger in seconds
t_al = 2; %Alighting time per passenger in seconds
stsk = zeros(1,n_s);
bsk = zeros(1,n_b);

gam_k = 2; gam_theta = 2; %parameters of Gamma rv

cumtime = zeros(1,n_b);

state0 = zeros(n_st,n_b);
state0(1,:) = ones(1,n_b);
state = state0;

%Busses will depart with headay of 90 seconds
%We will work on the basis of time a bus required to reach the next stop.
%It will we stored in the array n_xt. Note that this array will change as
%the code proceeds and its size will also change.
t_nxt = zeros(1,n_b); %this variable stores the time bus will require to reach the upcoming stop
for i=1:n_b  %At start each bus leaves with an interval of 363 seconds
    t_nxt(i) = 363*(i-1) + dis_stp(1)/v_bus;
end
%disp(t_nxt)
l_action = zeros(1,n_b); %last action taken by the bus
lapass = zeros(1,n_b);
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
    if state(1,1) == 8 %count contains whether the leading module has covered a circle or not
        count = 0;
    end
    gencount = gencount + 1;

    hw = headwaycuu(state, t_nxt, atstop, dis_stp, v_bus, n_b,count);

%     if gencount > 420 && gencount < 450
%             
%                 disp(gencount) 
%                 disp(hw)
%                 disp(state)
%                 disp(atstop)
%                 disp(t_nxt)               
%                 %break
%             
%     end
    for i = hw
        if i< 0
            disp('danger')
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
        
        
%         if l_action(im) ==1
%             lapass = state(2,im)*a_par(state(1,im));  % this is a deterministic quantity
%         else
%             lapass = 0;
%         end
        
        %acap: accomodating ,kjyjtcapacity
        acap = cap_bus - state(2,im);
    
        %lpass: for now ignoring the left over passengers
    %     if l_action(im+1) ==1
    %         lpass = arr_par(im);
    %     end
        %Computing rewards
        lpass = 0;

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

        if hw(im) < gamma*headwayt || hw(im) == gamma*headwayt
                stsk(state(1,im)) = 0;  %limit on a stop being skipped
                bsk(im) = 0;
                action = 'stop';
                l_action(im) = 0;
                pdsr = lapass(im) +  binornd((state(2,im)-lapass(im)),a_par(state(1,im)));
                lapass(im) = 0;
                pd_cum = pd_cum + pdsr;                
                state(2,im) = state(2,im) - pdsr;
                pa = poissrnd(arr_par(state(1,im))*hw(im));
                %disp(pa)
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

                 
        elseif hw(im) > gamma*headwayt && stsk(state(1,im))==0 && bsk(im) == 0  %If skip action is taken then the time spent at the stop will be zero
                pw = binornd(state(2,im),a_par(state(1,im))); %Number of passengers have to walk
                pw_cum = pw_cum + pw;
                lapass(im) = pw;
                stsk(state(1,im)) = 1;
                bsk(im) = 1;
                skcount = skcount + 1;
                action = 'skip';
                tspent = 0 + ex_wt;
                t_nxt(im) = dis_stp(state(1,im))/v_bus + ex_wt;
        end
         
     
    
    else %If the bus is not at the stop we will evalueate join possibility
                action = 'nextbs';
                t_nxt(im) = dis_stp(state(1,im))/v_bus;
                atstop(im) = 0;
                tspent = dis_stp(state(1,im))/v_bus;  %this tspent can be the t_nxt two lines above
    end           

    
    time  = time + M;
    %Now we construct arrays that consist the bus locations at different
    %time instants
%     if gencount > 3 && gencount < 6
%         disp(arr_par(state(1,im)))
%         disp(hw)
%     end
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
    %count = count + 1;

%     if gencount < 6        
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
%now finding the area under the curves
Pb_cum_int = trapz(Time, Pb_cum);
Pd_cum_int = trapz(Time, Pd_cum);
Pa_cum_int = trapz(Time, Pa_cum);
Pw_cum_int = trapz(Time, Pw_cum);
% for i=(Pb_cum - Pd_cum)
%     if i<0
%         disp('scam')
%     end
% end
%avg_inveh = (Pb_cum_int - Pd_cum_int)/Pb_cum(size(Pb_cum,2));
avg_inveh = (Pb_cum_int - Pd_cum_int)/Pb_cum(size(Pb_cum,2));

%avg_wait = (Pa_cum_int - Pb_cum_int)/Pa_cum(size(Pa_cum,2));
avg_wait = (Pa_cum_int - Pb_cum_int)/Pa_cum(size(Pa_cum,2));

avg_walk = Pw_cum_int/Pd_cum(size(Pd_cum,2));
fprintf('Average in-vehicle time: %f \n', avg_inveh)
fprintf('Average waiting time: %f \n', avg_wait)
fprintf('Average walk time: %f \n', avg_walk)
fprintf('Policy cost: %f \n', w_wait*avg_wait + avg_inveh + w_walk*avg_walk)

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
%legend('boarding', 'de-boarding')
% hold on
% plot(Time, State2(1,:))
%plot()