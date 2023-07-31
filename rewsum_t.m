function bestcomb = rewsum_t(k, state, t_nxt, atstop, t_al, t_bo, arr_par, a_par, dis_stp, im, hway, v_bus, v_pas, lpass, lapass, cap_bus,hwt,ord,count, gencount, reachtime, leavetime, time, Hway, Pa_all_cum, fPas1, fPas2, jP1, jP2, ex_wt)
    %k is the length of the horizon
    %this function computes the action to be taken in the current state
    %such that the cumulative reward over horizon 'k' is maximized if
    %average formulation is considered
    %total possible states are 5^k
    dis_f = 1; %discount factor
    n_pos = 5^k;
    state0 = state;
    minT = inf;
%     if gencount == 2621
%         srno = 1:size(state,2);
%         disp(im)
%         disp(table(srno', state', atstop', t_nxt' ))
%     end
    n_s = size(dis_stp,2);
    for pos = 0:n_pos -1
        tR = reachtime; tL = leavetime;
        tSim = time;
        simstate = state0; n_b = size(simstate,2);
        t_nxt_sim = t_nxt; atstop_sim = atstop; lapass_sim = lapass; lpass_sim = lpass;
        im_sim = im; 
        pos_code = dec2base(pos, 5);
        
        hw = headwaycuuu(simstate, t_nxt_sim, atstop_sim, dis_stp, v_bus, n_b,10,ord);
        if length(pos_code) < k
            ldiff = k - length(pos_code);
            pos_code = [repmat('0', 1, ldiff), pos_code]; %if the length of the pos_code is less than k then this line just adds zeros
        end
       
            %now we evaluate cumulative reward of all possibilities
        cumrew = 0;
        %disp(size(state0,2))
        for actn = 1:length(pos_code)
            
            n_b = size(simstate,2);

            all_a = allowed_actions(simstate, im_sim, atstop_sim);

            if contains(all_a,pos_code(actn)) %if the action is allowed in the current state
                implus = iplus(im_sim, n_b); implusp = iplus(implus, n_b);
                gone_else = false;
               
                if actn == 1
                    hw(im_sim) = hway; %using actyual value of headway at the first time instant
                end
                if count == 1 %this is to ensure the equlibrium state at the start
                    hw(1) = hwt;
                end
              
                act = str2num(pos_code(actn));

                
                if simstate(1, implus) == n_s && atstop_sim(implus) == 1
                    hw(implus) = hwt;
                end
               
                 

                rew = Reward_t(im_sim, simstate, act, a_par,arr_par, dis_stp, v_pas, ...
                hw(im_sim),cap_bus,lpass_sim(simstate(1,im_sim)),0, lapass_sim(im_sim), hw(implus), t_bo, t_al,hw(implusp), hwt, count,v_bus, gencount, pos_code, actn, fPas1, fPas2, jP1, jP2, ex_wt);
%                 if gencount == 2621 && pos_code(1) == '4' && pos_code(2) == '0' && pos_code(3) == '0' 
%                     disp(simstate)
%                     disp(atstop_sim)
%                     fprintf('pos_code : %s on module : %i, reward : %f \n', pos_code(actn), im_sim, rew)
%                 end
                cumrew = cumrew + rew*dis_f^(actn - 1);
%                 if pos_code(1) == '0' && pos_code(2) == '0' && pos_code(3) == '0' && pos_code(4) == '4' && gencount == 2001           
%                     fprintf('(from rewsum_t) action : %i, im :%i,headway: %f, reward : %f, actn ; %i \n', act, im_sim, hw(im_sim), rew, actn)                    
%                     %disp(table(simstate', atstop_sim', t_nxt_sim', lapass_sim'))
%                 end
                if simstate(3,im_sim) == 1
                    unit_cap = cap_bus*2;
                else
                    unit_cap = cap_bus;
                end
                if atstop_sim(im_sim) == 0
                    reaching = true;
                else
                    leaving = true;
                end
%                 disp('Hway in rewsum')
%                 fprintf('gencount : %i \n', gencount)
%                 disp(Hway)
                [n_state, t_nxt_sim, atstop_sim, lapass_sim, lpass_sim,tR, tL, Hway] = nextstpredImp(im_sim, simstate, act, t_nxt_sim, atstop_sim, hw(im_sim), arr_par, a_par, unit_cap, t_al, t_bo, dis_stp, lapass_sim, lpass_sim, v_bus,ord,tR, tL, tSim, pos_code, actn, gencount,hw,Hway, Pa_all_cum);
                simstate = n_state;
%                 if (gencount ==3) && (pos_code(1) == '4' || pos_code(1) == '3') %&& pos_code(2) == '0'  %&& actn == 1
%                     fprintf('pos_code : %s, gencount : %i, actn : %i \n', pos_code, gencount, actn)
%                     nums = 1:size(simstate,2);
%                     disp(table(nums', simstate', atstop_sim', t_nxt_sim', lapass_sim'))
%                     disp(lpass_sim)
%                 end
                hw = headwaycuuu(simstate, t_nxt_sim, atstop_sim, dis_stp, v_bus, n_b,10,ord);
                [~, im_sim] = min(t_nxt_sim);
                t_nxt_sim = t_nxt_sim - min(t_nxt_sim);
                
                tSim = tSim + min(t_nxt_sim);
                
                if atstop_sim(im_sim) == 0 %updating the state of the bus that have reached at the next stop
                    simstate(1,im_sim) = iplus(simstate(1,im_sim), n_s);
                end
            else
                cumrew = 1e10; %just setting to some very +ve value
                gone_else = true;
                break
            end
        end
%         if gencount == 3179 % && (pos_code(1) == '3' || pos_code(1) == '4')
%             fprintf('pos_code : %s, cumrew : %f \n', pos_code, cumrew)
%         end
        if cumrew < minT           
            minT = cumrew;
            bestcomb = pos_code;
        end
%         if pos_code(1) == '0' && pos_code(2) == '0' && pos_code(3) == '0' && pos_code(4) == '4' && gencount == 501
%            
%             disp('from rewsum_t for hor')
%             disp(cumrew)
%         end
%         if gencount == 78
%             fprintf('pos_code: %s, rew : %f \n', pos_code, cumrew)
%             
%         end
    end
%     if gencount > 2490 && gencount <2498
%         fprintf('gencount : %i, bestcomb : %s \n', gencount, bestcomb)
%     end
    %fprintf('maxrew %f \n', maxrew)
%     if bestcomb(1) == '3'
%         disp('found')
%         disp(gencount)
%     end
%     if gencount == 501
%         disp('from rewsum_t')
%         disp(minT)
%     end
end
                





    
    