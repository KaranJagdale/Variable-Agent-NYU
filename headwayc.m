function hw = headwayc(state, nxt, atstop, dis_stp, v_bus, n_b,count)
    hw = zeros(1,n_b);
       
        for i=2:n_b
            
            if atstop(i)==0 && atstop(i-1) == 0
                k = 0;
                for j=state(1,i)+1 : state(1,i-1)
                    k = k+dis_stp(j);
                end
                hw(i) = nxt(i) + k/v_bus - nxt(i-1);
            end

            if atstop(i)==1 && atstop(i-1) == 0
                k = 0;
                %disp(atstop)
                for j=state(1,i) : state(1,i-1) 
                    %disp(j)
                    k = k+dis_stp(j);
                end
                hw(i) = nxt(i) + k/v_bus - nxt(i-1);
            end

            if atstop(i)==0 && atstop(i-1) == 1
                k = 0;
                for j=state(1,i)+1 : state(1,i-1) -1
                    k = k+dis_stp(j);
                end
                hw(i) = nxt(i) + k/v_bus;
            end

            if atstop(i)==1 && atstop(i-1) == 1
                k = 0;
                for j=state(1,i)+1 : state(1,i-1) -1
                    k = k+dis_stp(j);
                end
                hw(i) = nxt(i) + k/v_bus;
            end

            
        end
        
    
        %It is assumed that the passengers come at the bus stop for the first
        %time in the day more or less at the same time. The busses are planned
        %such that for the first bus, the passengers start coming to stop i
        %when the bus leaves the stop i-1.
        if count ==1
            hw(1) = dis_stp(state(1))/v_bus;
        else
            if atstop(1)==0 && atstop(n_b) == 0
                k = 0;
                for j=state(1,1)+1 : state(1,n_b) +1
                    k = k+dis_stp(j);
                end
                hw(1) = nxt(1) + k/v_bus - nxt(n_b);
            end

            if atstop(1)==1 && atstop(n_b) == 0
                k = 0;
                for j=state(1,1) : state(1,n_b) +1
                    k = k+dis_stp(j);
                end
                hw(1) = nxt(1) + k/v_bus - nxt(n_b);
            end

            if atstop(1)==0 && atstop(n_b) == 1
                k = 0;
                for j=state(1,1)+1 : state(1,n_b)
                    k = k+dis_stp(j);
                end
                hw(1) = nxt(1) + k/v_bus;
            end

            if atstop(1)==1 && atstop(n_b) == 1
                k = 0;
                for j=state(1,1) : state(1,n_b) 
                    k = k+dis_stp(j);
                end
                hw(1) = nxt(1) + k/v_bus;
            end

        end
end