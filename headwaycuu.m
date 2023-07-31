function hw = headwaycuu(state, nxt, atstop, dis_stp, v_bus, n_b,count)
    hw = zeros(1,n_b);
%         disp('atstop')
%         disp(atstop)
        for i=2:n_b
            
            if atstop(i)==0 && atstop(i-1) == 0
                k = 0;
                if state(1,i) < state(1,i-1) 
                    for j=state(1,i) + 1 : state(1,i-1)
                        k = k+dis_stp(j);
                    end
                    hw(i) = nxt(i) + k/v_bus - nxt(i-1);
                elseif state(1,i) == state(1,i-1)
                    if atstop(i) == 1 && atstop(i-1) == 0
                        for j=state(1,i) + 1 : state(1,i-1)
                            k = k+dis_stp(j);
                        end
                        hw(i) = nxt(i) + k/v_bus - nxt(i-1);
                    elseif atstop(i) == 0 && atstop(i-1) == 1                     
                        for j=state(1,i-1) +1 : state(1,i)
                            k = k+dis_stp(j);
                        end
                        k = sum(dis_stp) - k;
                        hw(i) = nxt(i) + k/v_bus - nxt(i-1);
                    else
                        if nxt(i)>nxt(i-1)
                            for j=state(1,i) + 1 : state(1,i-1)
                                k = k+dis_stp(j);
                            end
                            hw(i) = nxt(i) + k/v_bus - nxt(i-1);
                        else
                            for j=state(1,i-1) +1 : state(1,i)
                                k = k+dis_stp(j);
                            end
                            k = sum(dis_stp) - k;
                            hw(i) = nxt(i) + k/v_bus - nxt(i-1);
                        end
                    end
                else
                    for j=state(1,i-1) +1 : state(1,i)
                            k = k+dis_stp(j);
                    end
                    k = sum(dis_stp) - k;
                    hw(i) = nxt(i) + k/v_bus - nxt(i-1);
                end


            elseif atstop(i)==1 && atstop(i-1) == 0
                k = 0;
                if state(1,i) < state(1,i-1) 
                    for j=state(1,i) : state(1,i-1)
                        k = k+dis_stp(j);
                    end
                    hw(i) = nxt(i) + k/v_bus - nxt(i-1);
                elseif state(1,i) == state(1,i-1)
                    if atstop(i) == 1 && atstop(i-1) == 0
                        for j=state(1,i) : state(1,i-1)
                            k = k+dis_stp(j);
                        end
                        hw(i) = nxt(i) + k/v_bus - nxt(i-1);
                    elseif atstop(i) == 0 && atstop(i-1) == 1                     
                        for j=state(1,i-1) +1 : state(1,i) -1
                            k = k+dis_stp(j);
                        end
                        k = sum(dis_stp) - k;
                        hw(i) = nxt(i) + k/v_bus - nxt(i-1);
                    else
                        if nxt(i)>nxt(i-1)
                            for j=state(1,i) : state(1,i-1)
                                k = k+dis_stp(j);
                            end
                            hw(i) = nxt(i) + k/v_bus - nxt(i-1);
                        else
                            for j=state(1,i-1) +1 : state(1,i) -1
                                k = k+dis_stp(j);
                            end
                            k = sum(dis_stp) - k;
                            hw(i) = nxt(i) + k/v_bus - nxt(i-1);
                        end
                    end
                else
                    for j=state(1,i-1) +1 : state(1,i) -1
                        k = k+dis_stp(j);
                    end
                    k = sum(dis_stp) - k;
                    hw(i) = nxt(i) + k/v_bus - nxt(i-1);
                end
           

            elseif atstop(i)==0 && atstop(i-1) == 1
                k = 0;
                if state(1,i) < state(1,i-1) 
                    for j=state(1,i)+1 : state(1,i-1) -1
                        k = k+dis_stp(j);
                    end
                    hw(i) = nxt(i) + k/v_bus;
                elseif state(1,i) == state(1,i-1)
                    if atstop(i) == 1 && atstop(i-1) == 0
                        for j=state(1,i)+1 : state(1,i-1) -1
                            k = k+dis_stp(j);
                        end
                        hw(i) = nxt(i) + k/v_bus;
                    elseif atstop(i) == 0 && atstop(i-1) == 1                     
                        for j=state(1,i-1) : state(1,i)
                            k = k+dis_stp(j);
                        end
                        k = sum(dis_stp) - k;
                        hw(i) = nxt(i) + k/v_bus;
                    else
                        if nxt(i)>nxt(i-1)
                            for j=state(1,i)+1 : state(1,i-1) -1
                                k = k+dis_stp(j);
                            end
                            hw(i) = nxt(i) + k/v_bus;
                        else
                            for j=state(1,i-1) : state(1,i)
                                k = k+dis_stp(j);
                            end
                            k = sum(dis_stp) - k;
                            hw(i) = nxt(i) + k/v_bus;
                        end
                    end
                else
                    for j=state(1,i-1) : state(1,i)
                        k = k+dis_stp(j);
                    end
                    k = sum(dis_stp) - k;
                    hw(i) = nxt(i) + k/v_bus;
                end

            

            elseif atstop(i)==1 && atstop(i-1) == 1
                k = 0;
                if state(1,i) < state(1,i-1) 
                    for j=state(1,i) : state(1,i-1) -1
                        k = k+dis_stp(j);
                    end
                    hw(i) = nxt(i) + k/v_bus;
                elseif state(1,i) == state(1,i-1)
                    if atstop(i) == 1 && atstop(i-1) == 0
                        for j=state(1,i) : state(1,i-1) -1
                            k = k+dis_stp(j);
                        end
                        hw(i) = nxt(i) + k/v_bus;
                    elseif atstop(i) == 0 && atstop(i-1) == 1                     
                        for j=state(1,i-1) : state(1,i)-1
                            k = k+dis_stp(j);
                        end
                        k = sum(dis_stp) - k;
                        hw(i) = nxt(i) + k/v_bus;
                    else
                        if nxt(i)>nxt(i-1)
                            for j=state(1,i) : state(1,i-1) -1
                                k = k+dis_stp(j);
                            end
                            hw(i) = nxt(i) + k/v_bus;
                        else
                            for j=state(1,i-1) : state(1,i)-1
                                k = k+dis_stp(j);
                            end
                            k = sum(dis_stp) - k;
                            hw(i) = nxt(i) + k/v_bus;
                        end
                    end
                else
                    for j=state(1,i-1) : state(1,i)-1
                        k = k+dis_stp(j);
                    end
                    k = sum(dis_stp) - k;
                    hw(i) = nxt(i) + k/v_bus;
                end
   
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
                if state(1,1) < state(1,n_b) 
                    for j=state(1,1)+1 : state(1,n_b)
                        k = k+dis_stp(j);
                    end
                    hw(1) = nxt(1) + k/v_bus - nxt(n_b);
                elseif state(1,1) == state(1,n_b)
                    if atstop(1) == 1 && atstop(n_b) == 0
                        for j=state(1,1)+1 : state(1,n_b)
                            k = k+dis_stp(j);
                        end
                        hw(1) = nxt(1) + k/v_bus - nxt(n_b);
                    elseif atstop(1) == 0 && atstop(n_b) == 1                   
                        for j = state(1,n_b)+1: state(1,1)
                            k = k + dis_stp(j);
                        end
                        k = sum(dis_stp) - k;
                        hw(1) = nxt(1) + k/v_bus - nxt(n_b);
                    else
                        if nxt(1) > nxt(n_b)
                            for j=state(1,1)+1 : state(1,n_b)
                                k = k+dis_stp(j);
                            end
                        hw(1) = nxt(1) + k/v_bus - nxt(n_b);
                        else
                            for j = state(1,n_b)+1: state(1,1)
                                k = k + dis_stp(j);
                            end
                            k = sum(dis_stp) - k;
                            hw(1) = nxt(1) + k/v_bus - nxt(n_b);
                        end
                    end
                else
                    for j = state(1,n_b)+1: state(1,1)
                        k = k + dis_stp(j);
                    end
                    k = sum(dis_stp) - k;
                    hw(1) = nxt(1) + k/v_bus - nxt(n_b);
                end
            

            elseif atstop(1)==1 && atstop(n_b) == 0
                k = 0;
                if state(1,1) < state(1,n_b) 
                    for j=state(1,1) : state(1,n_b)
                        k = k+dis_stp(j);
                    end
                    hw(1) = nxt(1) + k/v_bus -nxt(n_b);
                elseif state(1,1) == state(1,n_b)
                    if atstop(1) == 1 && atstop(n_b) == 0
                        for j=state(1,1) : state(1,n_b)
                            k = k+dis_stp(j);
                        end
                        hw(1) = nxt(1) + k/v_bus -nxt(n_b);
                    elseif atstop(1) == 0 && atstop(n_b) == 1                   
                        for j = state(1,n_b)+1: state(1,1)-1
                            k = k + dis_stp(j);
                        end
                        k = sum(dis_stp) - k;
                        hw(1) = nxt(1) + k/v_bus - nxt(n_b);
                    else
                        if nxt(1) > nxt(n_b)
                            for j=state(1,1) : state(1,n_b)
                                k = k+dis_stp(j);
                            end
                            hw(1) = nxt(1) + k/v_bus -nxt(n_b);
                        else
                            for j = state(1,n_b)+1: state(1,1)-1
                                k = k + dis_stp(j);
                            end
                            k = sum(dis_stp) - k;
                            hw(1) = nxt(1) + k/v_bus - nxt(n_b);
                        end
                    end
                else
                    for j = state(1,n_b)+1: state(1,1)-1
                        k = k + dis_stp(j);
                    end
                    k = sum(dis_stp) - k;
                    hw(1) = nxt(1) + k/v_bus - nxt(n_b);
                end

            elseif atstop(1)==0 && atstop(n_b) == 1
                k = 0;
                if state(1,1) < state(1,n_b) 
                    for j=state(1,1) +1 : state(1,n_b)-1
                        k = k+dis_stp(j);
                    end
                    hw(1) = nxt(1) + k/v_bus;
                elseif state(1,1) == state(1,n_b)
                    if atstop(1) == 1 && atstop(n_b) == 0
                        for j=state(1,1) +1 : state(1,n_b)-1
                            k = k+dis_stp(j);
                        end
                        hw(1) = nxt(1) + k/v_bus;
                    elseif atstop(1) == 0 && atstop(n_b) == 1                   
                        for j = state(1,n_b): state(1,1)
                            k = k + dis_stp(j);
                        end
                        k = sum(dis_stp) - k;
                        hw(1) = nxt(1) + k/v_bus;
                    else
                        if nxt(1) > nxt(n_b)
                            for j=state(1,1) +1 : state(1,n_b)-1
                                k = k+dis_stp(j);
                            end
                            hw(1) = nxt(1) + k/v_bus;
                        else
                            for j = state(1,n_b): state(1,1)
                                k = k + dis_stp(j);
                            end
                            k = sum(dis_stp) - k;
                            hw(1) = nxt(1) + k/v_bus;
                        end
                    
                    end
                else
                    for j = state(1,n_b): state(1,1)
                        k = k + dis_stp(j);
                    end
                    k = sum(dis_stp) - k;
                    hw(1) = nxt(1) + k/v_bus;
                end
            

            elseif atstop(1)==1 && atstop(n_b) == 1
                k = 0;
                if state(1,1) < state(1,n_b) 
                    for j=state(1,1) : state(1,n_b) -1
                        k = k+dis_stp(j);
                    end
                    hw(1) = nxt(1) + k/v_bus;
                elseif state(1,1) == state(1,n_b)
                    if atstop(1) == 1 && atstop(n_b) == 0
                        for j=state(1,1) : state(1,n_b) -1
                            k = k+dis_stp(j);
                        end
                        hw(1) = nxt(1) + k/v_bus;
                    elseif atstop(1) == 0 && atstop(n_b) == 1                   
                        for j = state(1,n_b): state(1,1)-1
                            k = k + dis_stp(j);
                        end
                        k = sum(dis_stp) - k;
                        hw(1) = nxt(1) + k/v_bus;
                    else
                        if nxt(1) > nxt(n_b)
                            for j=state(1,1) : state(1,n_b) -1
                                k = k+dis_stp(j);
                            end
                            hw(1) = nxt(1) + k/v_bus;
                        else
                            for j = state(1,n_b): state(1,1)-1
                                k = k + dis_stp(j);
                            end
                            k = sum(dis_stp) - k;
                            hw(1) = nxt(1) + k/v_bus;
                        end
                    end
                else
                    for j = state(1,n_b): state(1,1)-1
                        k = k + dis_stp(j);
                    end
                    k = sum(dis_stp) - k;
                    hw(1) = nxt(1) + k/v_bus;
                end
            end

        end
        