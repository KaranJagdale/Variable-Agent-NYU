 function reward = Reward(bus, state, action, a_par,arr_par, dis_stp, v_pas, ...
    hw,cap_bus,lpass,lsact, lapass, phead, t_bo, t_al)
    %In main code, need to calculate acap, lpass, lapass
    % 0 = "Stop", 1 = "Skip", 2 = "Split", 3 = "Join"
    %acap - accomodating capacity of the bus, i.e., capacity - load
    %Depending on if the bus is joint or a separate module, the value of
    %capacity will have to vary
    %lpass - left over passengers due to capacity constraint or skip actiom
    %of las bus
    %hw - headway of the bus
    %lapass - left over alighting passengers this has to be calculated in
    %the main code
%     switch lact
%         case 0
%             pbs = min(arr_par(state(1)+1)*hw + lpass,acap);
%         case 1
%             pbs = min(arr_par(state(1)+1)*hw + lpass, acap);
%         case 2
%             pbs = min(arr_par(state(1)+1)*hw + lpass,acap);
%         case 3
%             pbs = min(arr_par(state(1)+1)*hw + lpass,acap);
%     end    
    n_s = 8;
    if state(3) == 1
        cap_bus = cap_bus*2;
    end
    hw = hw(bus);
    hwt = 200; %target headway
    pbs = min(arr_par(state(1))*hw + lpass,cap_bus-state(2)); %this is approximation
    %nxt_stp is the next stop to the stop to which the bus is approaching
    if state(1) == n_s
        nxt_stp = 1;
    else 
        nxt_stp = state(1) + 1;
    end

    
    
    %nnext_stp is next to next stop to the stop which bus is approaching
    if state(1) == n_s-1
        nnxt_stp = 1;
    elseif state(1) == n_s
        nnxt_stp = 2;
    else        
        nnxt_stp = state(1) + 2;
    end
         
%     switch lsact
%         case 0
%             pds = a_par(bus)*state(2);
%         case 1
%             pds = a_par(bus)*(state(2)-lapass) + lapass;  
%         case 2
%             pds = a_par(bus)*state(2);
%         case 3
%             pds = a_par(bus)*state(2);
%     end
    if lsact == 1
        pds = a_par(bus)*(state(2)-lapass) + lapass;
    else
        pds = a_par(bus)*state(2);
    end

    pbsplus = min(arr_par(nxt_stp)*hw + lpass,cap_bus - (state(2)+pbs - pds));  %assuming deboarding people = boarding people
    pdsplus = a_par(nxt_stp)*(state(2)-pds + pbs); 
  
    fixJr = 1500; %this is fixed positive reward for join action since the logical computation is not possible
    %here the fix positive reward is used as joining the bus allows it to
    %split on some next station
    wst = dis_stp(nxt_stp)/v_pas;
    wstp  = dis_stp(nnxt_stp)/v_pas;
    rskip = -pds*wst - pbs*phead + (state(2)-pds)*(pds*t_al + pbs*t_bo) + 10*(hw - hwt);% + pbsplus*(pds*t_al + pbs*t_bo);
    rstop = -rskip;
    rsplit = (state(2)-pds)*(pds*t_al + pbs*t_bo) +pds*wst +pbs*phead + 10*(hw - hwt); % pbsplus*(pds*t_al + pbs*t_bo);
    %disp(pbsplus)
    rskipp = -pdsplus*wstp - pbsplus*phead + (state(2)-pdsplus)*(pdsplus*t_al + pbsplus*t_bo) + 10*(hw - hwt);
    rsplitp = (state(2)-pdsplus)*(pdsplus*t_al + pbsplus*t_bo) +pdsplus*wst +pbsplus*phead + 10*(hw - hwt);
    R1 = rsplitp;
    R2 = abs(rskipp);
    rjoin = R1- R2  - phead*state(2);
    
    switch action
        case 0
            reward = rstop;
        case 1
            reward = rskip;
        case 2
            reward = rsplit;
        case 3
            reward = rjoin;
    end
%Reward for action 3(join) is probably impossible to compute as it depends
%on the future load which depends on the current action which depends on
%the current reward computation. So, writting the reward for action 3(join)
%is a circular problem. Currenct reward for action 3(join) is written
%assuming that current action taken is stop. So, may be we can replace the
%poaitive part of reward with some constant instead of the current
%formulation. Need to think on this front.


