function [t_nxt, atstop, state, pa_cum, pb_cum, pd_cum, pe_cum, lapass, lpass, spcount, n_b, w_pass, tw_pass] = nextst(n_st, n_s, a_par,arr_par, im,hway, state, ...
    pa_cum, pb_cum,pd_cum, pe_cum, lapass, lpass, dis_stp, unit_cap, reachtime, Pa_all_cum, time, spcount,n_b, w_pass, tw_pass, t_nxt, atstop, ex_wt)
            envmainprop;
            fixdt = 30;
            rng(1)
            action  = 'split';
            spcount = spcount +1;
            n_b = n_b +1; %The number of modules increase due to splitting
           
            pds_c = binornd((state(2,im)-lapass(im)),a_par(state(1,im))) + lapass(im);
            disp(state(2,im))
            disp(a_par)
            fprintf('pdc : %i, lapass(im) : %i a_par : %f \n', pds_c, lapass(im), a_par(state(1,im)))
            stmin = iminus(state(1,im), n_s);
            for iw=1:4
                if w_pass(iw,stmin) == 0
                    break
                end
            end

            w_pass(iw,stmin) = lapass(im); %If the deboarding passengers are more than unit capacity then the left over passengers gets priority for deboarding
            if w_pass(iw,stmin) ~=0

                tw_pass(iw,stmin) = dis_stp(stmin)/v_pas;
            end
            
            load_r = pds_c;  %floor is not necessary but having it will not affect as well
            lapassf = 0;
            if state(2,im) - load_r > unit_cap
                load_r = state(2,im) - unit_cap;
            end
            
            if load_r > unit_cap
                lapassf = load_r - unit_cap;
                load_r = unit_cap;
                pds_c = load_r;    %as all passengers cannot deboard in this case               
            end
            pdsr = pds_c;
            pd_cum = pd_cum + pdsr;
            pe_cum = pe_cum + pdsr - lapass(im);
%             if gencount == 167
%                 disp(load_r)
%             end
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
            atstop_i(im + 1) = 1; %as rear module has stopped at the stop
            for i = im + 2 : n_b
                atstop_i(i) = atstop(i-1);
            end
            atstop = atstop_i;

            lapass_i = zeros(1,n_b);
            for i=1:im-1
                lapass_i(i) = lapass(i);
            end
            lapass_i(im) = lapassf;
            lapass_i(im+1) = 0;
            for i = im+2:n_b
                lapass_i(i) = lapass(i-1);
            end
            lapass = lapass_i; 

            state(2,im+1) = state(2,im+1) - pdsr; %This will be zero because the way we have distributed the passengers

            f_ind = floor(reachtime(state(1,im))); s_ind = floor(time);
            if s_ind >= size(Pa_all_cum,2)
                    s_ind = size(Pa_all_cum,2)-1;
            end
            %padec(state(1,im), f_ind+1:s_ind+1) = 1;
            pa = Pa_all_cum(state(1,im),s_ind+1) - Pa_all_cum(state(1,im),f_ind+1);
            pa_cum = pa_cum + pa;
            
            pbsr = min(pa + lpass(state(1,im)), unit_cap - state(2,im+1));
            pb_cum = pb_cum + pbsr;
            lpass(state(1,im)) = pa + lpass(state(1,im)) - pbsr; %assuming the leftover passengers gets priority boarding and they are never more than unitcap
            %papbcum = papbcum + pa-pbsr;

            state(2,im+1) = state(2,im+1) + pbsr;
            tspent = t_al*pdsr + t_bo*pbsr + fixdt + ex_wt + splittime;
           
            t_nxt_i = zeros(1,n_b);
            for i = 1:im-1
                t_nxt_i(i) = t_nxt(i);
            end

            t_nxt_i(im) = dis_stp(state(1,im))/v_bus + ex_wt + splittime; %As the leading module does not stop.
            t_nxt_i(im + 1) = tspent; %time to leave the stop
            
            for i = im+2:n_b
                t_nxt_i(i) = t_nxt(i-1);
            end
            t_nxt = t_nxt_i;
end