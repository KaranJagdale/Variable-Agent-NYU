 function reward = Rewarduhwi_i(bus, state, action, a_par,arr_par, dis_stp, v_pas, ...
    hw,cap_bus,lpass,lsact, lapass, phead, t_bo, t_al, pphead, hwt, count, v_bus, gencount, pos_code, actn)
    %In main code, need to calculate acap, lpass, lapass
    % 0 = "Stop", 1 = "Skip", 2 = "Split", 3 = "Join"

    n_b = size(state,2);
    ipl = iplus(bus,n_b); %bus behind the current bus
    imin = iminus(bus, n_b);
    iminm = iminus(imin, n_b);
    stateff = state(:,iminm);
    statef = state(:,imin);
    stateb = state(:,ipl); %State is state of all the modules
    state = state(:,bus);
    n_s = 8;
    fixdt  = 30;
    if state(3) == 1
        cap_bus = cap_bus*2;
    end
    
    %hw = hw(bus);
    %hwt = 150; %target headway
    hwg = 10*0; %headway gain
    fg = 0.5;
    
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
    if count == 1
        bstf = state(1);
        st_sep = floor(n_s/n_b);
        for i = 1:st_sep
            bstf = iminus(bstf,n_b);
        end
        lambs = lambsum(state(1), bstf, arr_par, n_s);
    else
        lambs =lambsum(state(1), statef(1), arr_par, n_s); %lambs accountss the number of passengers arriving at the stops between two successive buses
    end
    stop_time = (pds*t_al + pbs*t_bo + fixdt);
    wst = dis_stp(state(1))/v_pas;
    wstp  = dis_stp(nxt_stp)/v_pas;
    %fprintf('hw : %f, st : %f, phead : %f \n', hw, stop_time,  phead)
    rskip = -pds*wst*w_walk - pbs*phead*w_wait + (state(2)-pds)*stop_time + fg*stop_time*lambs*(hw - stop_time - phead) + hwg*(hw - hwt);% + pbsplus*(pds*t_al + pbs*t_bo);
    
    rstop = -rskip;
    
    rsplit = (state(2)-pds)*stop_time + pds*wst*w_walk + pbs*phead*w_wait + fg*stop_time*lambs*(hw - stop_time + phead) + hwg*(hw - hwt); % pbsplus*(pds*t_al + pbs*t_bo);
%     if gencount == 31 && pos_code(1) == '0' && pos_code(2) == '2' && actn == 2
%         fprintf('reward stop : %f \n', rstop)
%         fprintf('reward split : %f \n', rsplit)
%         disp((state(2)-pds)*stop_time)
%         disp(fg*stop_time*lambs*(hw - stop_time) )
%         disp(fg*stop_time*lambs*phead)
%         fprintf('phead : %f, hw : %f, s_t : %f \n', phead, hw, stop_time)
% 
%     end
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
        phead_pred = phead;
%         if gencount == 26 && (pos_code(1) == '3' && pos_code(2) == '4') && actn == 2
%             pphead = 93;
%         end
%         if gencount == 26 && ((pos_code(1) == '4' && pos_code(2) == '4') || (pos_code(1) == '3' && pos_code(2) == '4')) && actn ==2 %&& actn ==1
%             fprintf('phead : %f, pphead : %f \n', phead, pphead)
%         end
        
        for i=stbkbus
            p_ds = a_par(i)*stateb(2);
            l_i= stateb(2) - p_ds;  %intermediate load
            p_bs = min(phead*arr_par(i),cap_bus - l_i);
            %stateb(2) = stateb(2) + p_bs;
            nxt_stp = iplus(i,size(a_par,2));
            wst = dis_stp(nxt_stp)/v_pas;
            st_time = (p_ds*t_al + p_bs*t_bo + fixdt);
            lambs =lambsum(i, state(1), arr_par, n_s);
            
            rjskip = -p_ds*wst*w_walk - p_bs*pphead*w_wait + (stateb(2)-p_ds)*st_time + fg*(phead_pred - st_time - pphead)*st_time*lambs;% + 10*(phead - hwt);
            pphead = pphead + st_time;
            if  phead_pred - dis_stp(i)/v_bus- st_time  > 0
                phead_pred = phead_pred - dis_stp(i)/v_bus -st_time;
            else
                phead_pred = 0;
            end
            if gencount == 26 && pos_code(1) == '4' && pos_code(2) == '4' && actn ==2
                fprintf('rjskip : %f \n', rjskip)
            end
            rjskipc = rjskipc + rjskip;
            %disp(rjskip)
        end
        %disp(rjskipc)
        %above is the cumulative skip reward as the back module has to skip some stops to catch up with the front module 
        %lambs =lambsum(statf(1), stateff(1), arr_par, n_s);
        stop_time_p = (pdsplus*t_al + pbsplus*t_bo + fixdt);
        rskipp = -pdsplus*wstp*w_walk - pbsplus*phead*w_wait + (state(2)-pdsplus)*stop_time_p + stop_time_p*lambs*(hw - stop_time - phead) +  hwg*(hw - hwt);
        rsplitp = (state(2)-pdsplus)*stop_time_p +pdsplus*wst*w_walk + pbsplus*phead*w_wait + stop_time_p*lambs*(hw - stop_time + phead) + hwg*(hw - hwt);
        R1 = rsplitp*0;
        R2 = abs(rskipp)*0;
        rjoin = R1- R2 - phead*state(2) + rjskipc - hwg*(hw - hwt);
    end
%     if gencount == 26 && ((pos_code(1) == '4' && pos_code(2) == '4') || (pos_code(1) == '3' && pos_code(2) == '4')) && actn ==2 %&& actn ==1
%         fprintf('pos_code : %s \n', pos_code)
%         fprintf('skip con : %f, pass con : %f, phead :  %f \n', rjskipc, phead*state(2), phead)
%         disp(stbkbus)
%         fprintf('stateb(1) : %i, state(1) : %i \n', stateb(1), state(1))
%         fprintf('im : %i, imp : %i \n', bus, ipl)
%         %fprintf('bus : %i \n', bus)
%     end
    
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


