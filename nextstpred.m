function [n_state, t_nxt, atstop, lapass, lpass, tR, tL, Hway] = nextstpred(im, state, act, t_nxt, atstop, hway, arr_par, a_par, unit_cap, t_al, t_bo, dis_stp, lapass, lpass, v_bus,ord,tR, tL, tSim, pos_code, actn, gencount, hw, Hway)
    splittime = 1; fixdt = 30;
    n_b = size(state,2);
%     disp('Hway in nextstpred')
%     disp(Hway)
    env_conhor;
    n_s = size(dis_stp,2);
    %fprintf('n_b before hw %i \n', n_b)
    %hw = headwaycuuu(state, t_nxt, atstop, dis_stp, v_bus, n_b,10,ord);     
    %fprintf('n_b after hw %i \n', n_b)

        
    st_reach = state(1,im);
    for i=1:n_b
            if i ~= im
                t_nxt(i) = t_nxt(i) - t_nxt(im);  %nxt for im will be calculated as the stop time at the stop
            end       
    end
    ex_wt = 0; %extra wait the bus has to do because of already reached buses
    for i=1:n_b
        if state(1,i) == st_reach && atstop(i) == 1
            %b_st = [b_st i];
            if t_nxt(i) > ex_wt 
                ex_wt = t_nxt(i);
            end
        end
    end

    if act == 0
        
        hway = tSim - tR(state(1,im));
        pds = round((state(2,im)-lapass(im))*a_par(state(1,im)) + lapass(im));
        pa = round(hway*arr_par(state(1,im)));
        pbs =  min(pa + lpass(state(1,im)), unit_cap - state(2,im) + pds);
%         if gencount ==  409 && pos_code(1) == '0' && pos_code(2) == '0' && pos_code(3) == '0' && actn == 1
%             fprintf('(from nextstpred) hway : %f, state(1,im) : %i \n', hway, state(1,im))
%             fprintf('(from nextstpred) pds : %i, pbs : %i \n', pds, pbs)
%         end
        lapass(im) = 0;
        lpass(state(1,im)) =  pa + lpass(state(1,im)) - pbs;
        tspent = pds*t_al + pbs*t_bo + fixdt + ex_wt;
        n_state = state;
        %n_state(1,im) = iplus(state(1,im),n_s);
        n_state(2,im) = round(state(2,im) - pds + pbs);
        t_nxt(im) = tspent;
        atstop(im) = 1;
        tR(1,state(1,im)) = tSim;
    elseif act == 1
        hway = tSim - tR(1,state(1,im));
        tspent = 1 + ex_wt;
        atstop(im) = 1;
        lapass(im) = lapass(im) + round((state(2,im) - lapass(im))*a_par(state(1,im)));
%         if gencount ==  443 && pos_code(1) == '1' && pos_code(2) == '1' && pos_code(3) == '4' && actn == 2
%             fprintf('(from nextstpred) hway : %f, arr_par(state(1,im)) : %f, state(1,im) : %i \n', hway, arr_par(state(1,im)), state(1,im))
%             %fprintf('(from nextstpred) pds : %i, pbs : %i \n', pds, pbs)
%         end
        lpass(state(1,im)) = lpass(state(1,im)) + round(hway*arr_par(state(1,im)));
        n_state = state;
        %n_state(1,im) = iplus(state(1,im),n_s);
        t_nxt(im) = tspent;
        tR(1,state(1,im)) = tSim;
    elseif act == 2
        hway = tSim - tR(1,state(1,im));
        n_b = n_b +1;
        pds = round(a_par(state(1,im))*(state(2,im)-lapass(im)) + lapass(im));
        if pds > state(2,im)/2
            pds = state(2,im)/2;
        end
        pa = round(hway*arr_par(state(1,im)));
        pbs =  min(pa + lpass(state(1,im)), unit_cap/2 - floor(state(2,im)/2) + pds); %unitcap here is twice the capacity of the module
        tspent = pds*t_al + pbs*t_bo + fixdt + splittime + ex_wt;
        if gencount ==785 && pos_code(1) == '0' && pos_code(2) == '2' && pos_code(3) == '4' && actn ==2
            fprintf('(from nextstpred) pds : %i, pbs : %i, ex_wt : %f \n', pds, pbs, ex_wt)
            %fprintf('a_par : %f, state(2,im) : %i, ')
        end
        %Updating Hway
%         Hway_i = zeros(1,n_b);
%         Hway_i(1:im) = Hway(1:im);
%         Hway_i(im+1) = tspent;
%         Hway_i(im+2:end) = Hway(im+1:end);
%         Hway = Hway_i;

        n_state = zeros(3, n_b);
        
        for i=1:im-1
                n_state(:,i) = state(:,i);
        end
        n_state(:,im) = [state(1,im);ceil(state(2,im)/2);0];
        n_state(:,im+1) = [state(1,im);floor(state(2,im)/2 - pds + pbs);0];
        for i=im+2:n_b
            n_state(:,i) = state(:,i-1);
        end
        lpass(state(1,im)) = lpass(state(1,im)) + pa - pbs;
        atstop_i = zeros(1,n_b);
        for i=1:im-1
            atstop_i(i) = atstop(i);
        end
        atstop_i(im) = 1;
        atstop_i(im + 1) = 1; %as rear module has stopped at the stop
        for i = im + 2 : n_b
            atstop_i(i) = atstop(i-1);
        end
        atstop = atstop_i;

        lapass_i = zeros(1,n_b);
        for i=1:im-1
            lapass_i(i) = lapass(i);
        end
        lapass_i(im) = 0;
        lapass_i(im+1) = 0;
        for i = im+2:n_b
            lapass_i(i) = lapass(i-1);
        end
        lapass = lapass_i; 

        t_nxt_i = zeros(1,n_b);
        for i = 1:im-1
            t_nxt_i(i) = t_nxt(i);
        end
        ct = dis_stp(state(1,im))/v_bus ; %+ gamrnd(gam_k, gam_theta) - gam_theta*gam_k;  %this tspent can be the t_nxt two lines above
        
        t_nxt_i(im) = ex_wt + 1; %adding 1 to ensure the module behind leaves later than the module in the front
        t_nxt_i(im + 1) = tspent; %time to leave the stop
        atstop(im) = 1;
        for i = im+2:n_b
            t_nxt_i(i) = t_nxt(i-1);
        end
        t_nxt = t_nxt_i;
        tR(1,state(1,im)) = tSim;

    elseif act == 3

       % fprintf('n_b in join %i \n', n_b)
        implus = iplus(im,n_b);
%         if gencount == 796 && pos_code(1) == '1' && pos_code(2) == '3' && pos_code(3) == '0' && actn == 2
%             fprintf('(from nextstpred) if cond : %i \n', atstop(implus) == 1 && state(1, implus) == state(1,im))
%         end
        if atstop(implus) == 1 && state(1, implus) == state(1,im)
            tspent = t_nxt(implus);  %Front bus has to wait for this time
        else
            tspent = hw(implus);
        end
        %implus = iplus(im, n_b);
        load_rj = state(2,implus); % this value will require to caculate lapass
        stop_rj = state(1,implus);
        lap_ac = 0;
        if stop_rj < state(1,im)
            bskip = stop_rj+1:state(1,im);
        else
            bskip = stop_rj+1:n_s;
            bskip = [bskip 1:state(1,im)];
        end
        for st = bskip
            lap = floor((load_rj - lapass(implus) - lap_ac)*a_par(st));
            lap_ac = lap_ac + lap;
        end
        lap_ac = floor(lap_ac);
        if im == n_b
                    
            n_b = n_b - 1;
            
            state_i = zeros(n_st,n_b);
            for i=2:im-1
                state_i(:,i-1) = state(:,i);
            end
            state_i(:,n_b) = [state(1,im);state(2,im)+state(2,implus);1];
            n_state = state_i;

            %Updating Hway
%             Hway_i = zeros(1,n_b);
%             Hway_i(1:n_b-1) = Hway(2,n_b);
%             Hway_i(n_b) = Hway(n_b+1) + tspent;
%             Hway = Hway_i;

            %updating t_nxt
            t_nxt_i = zeros(1,n_b);
            for i=2:n_b
                t_nxt_i(i-1) = t_nxt(i);
            end
            t_nxt_i(n_b) = tspent;%t_nxt(n_b + 1) + tspent;
            %t_nxt_i(n_b) =  tspent;
            t_nxt = t_nxt_i;
            
            lapass_i = zeros(1,n_b);
            for i =2:n_b
                lapass_i(i-1) = lapass(i);
            end
            lapass_i(n_b) = lapass(1)+ lapass(n_b+1) + lap_ac;% + pw;
            lapass = lapass_i;

            atstop_i = zeros(1,n_b);

            for i=2:n_b
                atstop_i(i-1) = atstop(i);
            end
            atstop_i(n_b) = 1;
            atstop = atstop_i;
        else
            n_b = n_b - 1;
            %As the behind module is going to skip the stops we
            %have to calculate the number of waiting passengers
            %Updating state
            
            state_i = zeros(n_st,n_b);
            for i=1:im-1
                state_i(:,i) = state(:,i);
            end
            state_i(:,im) = [state(1,im);state(2,im)+state(2,implus);1];
            for i=im+1:n_b
                state_i(:,i) = state(:,i+1);
            end
            n_state = state_i;
            
            %Updating Hway

%             Hway_i = zeros(1,n_b);
%             Hway_i(1:im-1) = Hway(1:im-1);
%             Hway_i(im) = Hway(im) + tspent;
%             Hway_i(im+1:end) = Hway(im+2:end);
%             Hway = Hway_i;
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

            lapass_i = zeros(1,n_b);
            for i=1:im-1
                lapass_i(i) = lapass(i);
            end
            lapass_i(im) = lapass(im) + lapass(im+1) + lap_ac;% + pw;
            for i=im+1:n_b
                lapass_i(i) = lapass(i+1);
            end
            lapass = lapass_i;
            
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
    else
        mxtraveltime = 0;
        for bus = 1:n_b %this takes care of overtaking
            if state(1,bus) == state(1,im) && atstop(bus) == 0
                if t_nxt(bus) > mxtraveltime
                    mxtraveltime = t_nxt(bus);
                end
            end
        end
        ct = dis_stp(state(1,im))/v_bus;
        tspent = max(ct, mxtraveltime) + 0.1; %bus leaving late reaches late due to addition of 1
        
        t_nxt(im) = tspent;
        n_state = state;
        atstop(im) = 0;
    end
%     if size(atstop) ==1
%         disp('issue')
%         fprintf('n_b %i \n', n_b) 
%         disp(state)
%     end
%     disp('nstate')
%     disp(n_state)
%     disp('state')
%     disp(state)
%     fprintf('im %i \n', im)
%     fprintf('act %i \n', act)

end

        
     



