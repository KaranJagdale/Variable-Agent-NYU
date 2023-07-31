function bestcomb = rewsum_expt(~, state, t_nxt, atstop, t_al, t_bo, arr_par, a_par, dis_stp, im, hway, v_bus, v_pas, lpass, lapass, cap_bus,hwt,ord,count, gencount)
    %k is the length of the horizon
    %this function computes the action to be taken in the current state
    %such that the cumulative reward over horizon 'k' is maximized if
    %average formulation is considered
    %total possible states are 5^k
    k = 1;
    dis_f = 1; %discount factor
    n_pos = 5^k;
    state0 = state;
    maxR = -inf;
    n_s = size(dis_stp,2);
    gamma = 1.5;
    bet = 1;
    for pos = 0:n_pos -1
        simstate = state0;
        t_nxt_sim = t_nxt; atstop_sim = atstop; lapass_sim = lapass; lpass_sim = lpass;
        im_sim = im;
        pos_code = dec2base(pos, 5);
        
        if length(pos_code) < k
            ldiff = k - length(pos_code);
            pos_code = [repmat('0', 1, ldiff), pos_code]; %if the length of the pos_code is less than k then this line just adds zeros
        end
       
            %now we evaluate cumulative reward of all possibilities
        cumrew = 0;
        for actn = 1:length(pos_code)
            
            n_b = size(simstate,2);

            all_a = allowed_actions(simstate, im_sim, atstop_sim);

            if contains(all_a,pos_code(actn)) %if the action is allowed in the current state
                implus = iplus(im_sim, n_b); implusp = iplus(implus, n_b);
                gone_else = false;

                hw = headwaycuuu(simstate, t_nxt_sim, atstop_sim, dis_stp, v_bus, n_b,10,ord);
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
               
                 
%                 if gencount == 21 && (pos_code(1) == '0' && pos_code(2) == '3')
%                     fprintf('im_sim : %i \n', im_sim)
%                     fprintf('implus : %i, hw(implus) : %f \n',implus, hw(implus))
%                 end
                rew = Reward_expt(im_sim, simstate, act,hw(im_sim),n_s, gamma, bet, hwt);
                cumrew = cumrew + rew*dis_f^(actn - 1);
%                 if gencount == 32 && (pos_code(1) == '0' || pos_code(1) == '1' || pos_code(1) == '2')
%                     fprintf('pos_code : %s, act : %i, rew : %f, module : %i  \n', pos_code, act, rew, im_sim)
%                     
%                 end
%                 if count == 1
%                     fprintf('act : %i, rew : %f \n', act, rew)
%                 end
                %computing next state
                if simstate(3,im_sim) == 1
                    unit_cap = cap_bus*2;
                else
                    unit_cap = cap_bus;
                end
%                 if gencount == 30 && pos_code(1) == '0' && pos_code(2) == '4' && actn == 1
%                     disp('state init')
%                     disp(simstate)
%                 end
                [n_state, t_nxt_sim, atstop_sim, lapass_sim, lpass_sim] = nextstpred(im_sim, simstate, act, t_nxt_sim, atstop_sim, hw(im_sim), arr_par, a_par, unit_cap, t_al, t_bo, dis_stp, lapass_sim, lpass_sim, v_bus,ord);
                simstate = n_state;
%                 if gencount == 30 && pos_code(1) == '0' && pos_code(2) == '4' && actn == 1
%                     disp('simstate')
%                     disp(simstate)
%                     disp('t_nxt_sim')
%                     disp(t_nxt_sim)
%                 end
                [~, im_sim] = min(t_nxt_sim);

            else
                cumrew = -1e8; %just setting to some very -ve value
                gone_else = true;
                break
            end
        end
%         if gencount == 21 && (pos_code(1) == '0' || pos_code(1) == '1')
%             fprintf('pos_code : %s, cumrew : %f \n', pos_code, cumrew)
%         end
        if cumrew > maxR           
            maxR = cumrew;
            bestcomb = pos_code;
        end
%         if gencount < 10
%             fprintf('pos_code: %s, rew : %f \n', pos_code, cumrew)
%             
%         end
    end
%     if gencount == 78
%         fprintf('bestcomb : %s \n', bestcomb)
%     end
    %fprintf('maxrew %f \n', maxrew)
%     if bestcomb(1) == '3'
%         disp('found')
%         disp(gencount)
%     end
end
                





    
    