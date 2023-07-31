 function reward = Rewarduhwi(bus, state, action, a_par,arr_par, dis_stp, v_pas, ...
    hw,cap_bus,lpass,lsact, lapass, phead, t_bo, t_al, pphead, hwt)
    %In main code, need to calculate acap, lpass, lapass
    % 0 = "Stop", 1 = "Skip", 2 = "Split", 3 = "Join"

    n_b = size(state,2);
    ipl = iplus(bus,n_b); %bus behind the current bus
    stateb = state(:,ipl); %State is state of all the modules
    state = state(:,bus);
    n_s = size(dis_stp,2);
    fixdt  = 30;
    if state(3) == 1
        cap_bus = cap_bus*2;
    end
    
    %hw = hw(bus);
    %hwt = 150; %target headway
    hwg = 10*0; %headway gain
    
    %nxt_stp is the next stop to the stop to which the bus is approaching.
    nxt_stp = iplus(state(1), n_s);
    
    %nnext_stp is next to next stop to the stop which bus is approaching
    nnxt_stp = iplus(nxt_stp,n_s);
    
    pds = a_par(state(1))*(state(2)-lapass) + lapass;
    pbst = arr_par(state(1))*hw + lpass;
    pbs = min(pbst,cap_bus-(state(2)-pds)); %this is approximation
    pdsplus = a_par(nxt_stp)*(state(2)-pds + pbs);
    pbsplus = min(arr_par(nxt_stp)*hw + lpass,cap_bus - (state(2)+pbs - pds - pdsplus));  %assuming deboarding people = boarding people
     
    w_wait = 2.1; w_walk = 2.2;
    %fixJr = 1500; %this is fixed positive reward for join action since the logical computation is not possible
    %here the fix positive reward is used as joining the bus allows it to
    %split on some next station
    wst = dis_stp(state(1))/v_pas;
    wstp  = dis_stp(nxt_stp)/v_pas;
    rskip = -pds*wst*w_walk - pbs*phead*w_wait + (state(2)-pds)*(pds*t_al + pbs*t_bo + fixdt) + hwg*(hw - hwt);% + pbsplus*(pds*t_al + pbs*t_bo);
    rstop = -rskip;
    rsplit = (state(2)-pds)*(pds*t_al + pbs*t_bo + fixdt) + pds*wst*w_walk + pbs*phead*w_wait + hwg*(hw - hwt); % pbsplus*(pds*t_al + pbs*t_bo);
%     fprintf('rskip : %f \n', rskip)
%     fprintf('rsplit : %f \n', rsplit)
    %fprintf('rskip : %f \n', rskip)
    %Now the reward for join
%     fprintf('(state(2)-pds)*(pds*t_al + pbs*t_bo) = %f \n',(state(2)-pds)*(pds*t_al + pbs*t_bo))
%     fprintf('pds*wst = %f \n',pds*wst)
%     fprintf('pbs*phead = %f \n',pbs*phead)
%     fprintf('10*(hw - hwt) = %f \n', 10*(hw - hwt))
        
    if action == 3 || action == 4   
        stbkbus = [];
        if stateb(1) < state(1) || stateb(1) == state(1)
            for i=stateb(1):state(1)
                stbkbus = [stbkbus i];
            end
        else
            stbkbus = 1:size(a_par,2);
            for i=state(1):stateb(1)
                stbkbus(stbkbus == i) = [];
            end
        end
        rjskipc = 0;
        
        for i=stbkbus
            p_ds = a_par(i)*stateb(2);
            l_i= stateb(2) - p_ds;  %intermediate load
            p_bs = min(phead*arr_par(i),cap_bus - l_i);
            %stateb(2) = stateb(2) + p_bs;
            nxt_stp = iplus(i,size(a_par,2));
            wst = dis_stp(nxt_stp)/v_pas;
            rjskip = -p_ds*wst*w_walk - p_bs*pphead*w_wait + (stateb(2)-p_ds)*(p_ds*t_al + p_bs*t_bo + fixdt);% + 10*(phead - hwt);
            rjskipc = rjskipc + rjskip;
            %disp(rjskip)
        end
        %disp(rjskipc)
        %above is the cumulative skip reward as the back module has to skip some stops to catch up with the front module 


        rskipp = -pdsplus*wstp*w_walk - pbsplus*phead*w_wait + (state(2)-pdsplus)*(pdsplus*t_al + pbsplus*t_bo + fixdt) + hwg*(hw - hwt);
        rsplitp = (state(2)-pdsplus)*(pdsplus*t_al + pbsplus*t_bo + fixdt) +pdsplus*wst*w_walk +pbsplus*phead*w_wait + hwg*(hw - hwt);
        R1 = rsplitp;
        R2 = abs(rskipp)*0;
        rjoin = R1- R2 - phead*state(2) + rjskipc - hwg*(hw - hwt);
    end
    
    switch action
        case 0
            reward = rstop;
        case 1
            reward = rskip;
        case 2
            reward = rsplit;
        case 3
            reward = rjoin;
        case 4
            if state(3) == 0
                reward = -rjoin;
            else
                reward = 0; %In this case only nextbs is available so assigning it zero reward.
            end
    end
 end
%Reward for action 3(join) is probably impossible to compute as it depends
%on the future load which depends on the current action which depends on
%the current reward computation. So, writting the reward for action 3(join)
%is a circular problem. Currenct reward for action 3(join) is written
%assuming that current action taken is stop. So, may be we can replace the
%poaitive part of reward with some constant instead of the current
%formulation. Need to think on this front.


