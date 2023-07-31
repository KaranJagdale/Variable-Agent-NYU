 function timeGain = Reward_t(bus, state, action, a_par,arr_par, dis_stp, v_pas, ...
    hw,cap_bus,lpass,lsact, lapass, phead, t_bo, t_al, pphead, hwt, count, v_bus, gencount, pos_code, actn, fPas1, fPas2, jP1, jP2, ex_wt)
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
    n_s = size(dis_stp,2);
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
%         bstf = state(1);
%         st_sep = floor(n_s/n_b);
%         for i = 1:st_sep
%             bstf = iminus(bstf,n_b);
%         e
%         fprintf('bstf : %i, state(1) : %i \n', bstf, state(1))
%         lambs = lambsum(state(1), bstf, arr_par, n_s);
        lambs = 0;
    else
        lambs =lambsum(state(1), statef(1), arr_par, n_s); %lambs accountss the number of passengers arriving at the stops between two successive buses
    end
    %stop_time = (pds*t_al + pbs*t_bo + fixdt);
    stop_time = max(pds*t_al , pbs*t_bo) + fixdt + ex_wt;
    wst = dis_stp(state(1))/v_pas;
    wstp  = dis_stp(nxt_stp)/v_pas;
%     fPas1 = 1; %1 is nomina value
%     fPas2 = 1; %0.2 is nominal value
    %fprintf('hw : %f, st : %f, phead : %f \n', hw, stop_time,  phead)
    iPas1 = 1; iPas2 = 1;
    %fPas1= 0; fPas2 = 1.5;

    rskip = (phead*pbs*w_wait + pds*wst*w_walk)*iPas1 + stop_time*lambs/2*phead*w_wait*fPas1;

    rstop = (state(2)-pds)*stop_time*w_wait*iPas2 + (lambs*hw*w_wait*stop_time/2)*fPas2;

    rsplit = 0;

    rnb = 0;
%     if gencount > 3000 && gencount < 3050 && (pos_code== '1' || pos_code == '0')
%         fprintf('pos_code, :%s, stop_time :%f, lambs*hw : %f, state(2)-pds :%f \n', pos_code(actn), stop_time, lambs*hw, state(2) - pds)
%         fprintf('state(1) : %i, statef(1) : %i, lsum : %f, hw : %f \n', state(1), statef(1), lambs, hw)
%         fprintf('term st 1 : %f, term st 2 : %f \n', (state(2)-pds)*stop_time, (lambs*hw*w_wait*stop_time/2))
%         fprintf('term sk 1 : %f, term sk 2 : %f \n', phead*pbs*w_wait + pds*wst*w_walk, stop_time*lambs/2*phead*w_wait)
%     end

    if action == 3 || action == 4
        stbkbus = [];
        if stateb(1) < state(1) || stateb(1) == state(1)
            for i=stateb(1):state(1)
                stbkbus = [stbkbus i];
            end
        else
            stbkbus = 1:size(a_par,2);
            for i=state(1)+1:stateb(1)
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
%         if gencount > 10 && gencount < 20
%             disp('stbkbus')
%             disp(stbkbus)
%             fprintf('bus: %i, buap : %i \n', state(1), stateb(1))
%         end
        for i=stbkbus
            p_ds = a_par(i)*stateb(2);
            l_i= stateb(2) - p_ds;  %intermediate load
            p_bs = min(phead*arr_par(i),cap_bus - l_i);
            %stateb(2) = stateb(2) + p_bs;
            nxt_stp = iplus(i,size(a_par,2));
            wst = dis_stp(nxt_stp)/v_pas;
            st_time = (p_ds*t_al + p_bs*t_bo + fixdt);
            %lambs =lambsum(i, state(1), arr_par, n_s);
            
            rjskip = pphead*p_bs*w_wait + p_ds*wst*w_walk;% + 10*(phead - hwt);
            pphead = pphead + st_time;
            if  phead_pred - dis_stp(i)/v_bus- st_time  > 0
                phead_pred = phead_pred - dis_stp(i)/v_bus -st_time;
            else
                phead_pred = 0;
            end
%             if gencount == 26 && pos_code(1) == '4' && pos_code(2) == '4' && actn ==2
%                 fprintf('rjskip : %f \n', rjskip)
%             end
            rjskipc = rjskipc + rjskip;
            %disp(rjskip)
        end    
%         if gencount == 3 && pos_code(1) == '3' && pos_code(2) == '0' && pos_code(3) == '4' && actn == 1
%             fprintf('pos_code : %s, hw : %f, phead :%f, rjskipc : %f \n', pos_code(actn), hw, phead, rjskipc)
%         end
        rjoin = state(2)*phead + jP1*rjskipc + jP2*hw*lambs*w_wait*phead;
    end
    
    switch action
        case 0
            timeGain = rstop;
        case 1
            timeGain = rskip;
        case 2
            timeGain = rsplit;
        case 3
            timeGain = rjoin;
        case 4
            timeGain = rnb;
    end
 end


