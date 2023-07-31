function hw = headwaycuuu(state, nxt, atstop, dis_stp, v_bus, n_b,count, ord)
    hw = zeros(1,n_b);
%         disp('atstop')
%         disp(atstop)
    n_b = size(state,2);
        tot_mod = n_b;
        for i=1:size(state,2)
            if state(3,i) == 1
                tot_mod = tot_mod + 1;
            end
        end

        for i=1:n_b
%             l = i;
%             for bus =1:size(state,2)
%                 if state(3,bus) == 1 && i > bus
%                     l = l + 1;
%                 end
%             end
%             io = ord(l);
%             io_n = iminus(io, tot_mod);
%             i_n = find(ord == io_n);
%             %fprintf('i_n before : %i, state size : %i \n', i_n, size(state,2))
%             %disp(state)
%             l = i_n;
%             fprintf('i_n: %i \n', i_n)
%             for bus = 1: size(state,2)
%                 if (state(3,bus) == 1) && (i_n > bus)
%                     l = l - 1;
%                 end
%             end
%             i_n = l;
            %fprintf('i : %i, i_n : %i \n', i, i_n)
            i_n = iminus(i, n_b);
            if atstop(i)==0 && atstop(i_n) == 0
                k = 0;
                if state(1,i) < state(1,i_n) 
                    for j=state(1,i) + 1 : state(1,i_n)
                        k = k+dis_stp(j);
                    end
                    hw(i) = nxt(i) + k/v_bus - nxt(i_n);
                elseif state(1,i) == state(1,i_n)
                    if atstop(i) == 1 && atstop(i_n) == 0
                        for j=state(1,i) + 1 : state(1,i_n)
                            k = k+dis_stp(j);
                        end
                        hw(i) = nxt(i) + k/v_bus - nxt(i_n);
                    elseif atstop(i) == 0 && atstop(i_n) == 1                     
                        for j=state(1,i_n) +1 : state(1,i)
                            k = k+dis_stp(j);
                        end
                        k = sum(dis_stp) - k;
                        hw(i) = nxt(i) + k/v_bus - nxt(i_n);
                    else
                        if nxt(i)>nxt(i_n)
                            for j=state(1,i) + 1 : state(1,i_n)
                                k = k+dis_stp(j);
                            end
                            hw(i) = nxt(i) + k/v_bus - nxt(i_n);
                        else
                            for j=state(1,i_n) +1 : state(1,i)
                                k = k+dis_stp(j);
                            end
                            k = sum(dis_stp) - k;
                            hw(i) = nxt(i) + k/v_bus - nxt(i_n);
                        end
                    end
                else
                    for j=state(1,i_n) +1 : state(1,i)
                            k = k+dis_stp(j);
                    end
                    k = sum(dis_stp) - k;
                    hw(i) = nxt(i) + k/v_bus - nxt(i_n);
                end


            elseif atstop(i)==1 && atstop(i_n) == 0
                k = 0;
                if state(1,i) < state(1,i_n) 
                    for j=state(1,i) : state(1,i_n)
                        k = k+dis_stp(j);
                    end
                    hw(i) = nxt(i) + k/v_bus - nxt(i_n);
                elseif state(1,i) == state(1,i_n)
                    if atstop(i) == 1 && atstop(i_n) == 0
                        for j=state(1,i) : state(1,i_n)
                            k = k+dis_stp(j);
                        end
                        hw(i) = nxt(i) + k/v_bus - nxt(i_n);
                    elseif atstop(i) == 0 && atstop(i_n) == 1                     
                        for j=state(1,i_n) +1 : state(1,i) -1
                            k = k+dis_stp(j);
                        end
                        k = sum(dis_stp) - k;
                        hw(i) = nxt(i) + k/v_bus - nxt(i_n);
                    else
                        if nxt(i)>nxt(i_n)
                            for j=state(1,i) : state(1,i_n)
                                k = k+dis_stp(j);
                            end
                            hw(i) = nxt(i) + k/v_bus - nxt(i_n);
                        else
                            for j=state(1,i_n) +1 : state(1,i) -1
                                k = k+dis_stp(j);
                            end
                            k = sum(dis_stp) - k;
                            hw(i) = nxt(i) + k/v_bus - nxt(i_n);
                        end
                    end
                else
                    for j=state(1,i_n) +1 : state(1,i) -1
                        k = k+dis_stp(j);
                    end
                    k = sum(dis_stp) - k;
                    hw(i) = nxt(i) + k/v_bus - nxt(i_n);
                end
           

            elseif atstop(i)==0 && atstop(i_n) == 1
                k = 0;
                if state(1,i) < state(1,i_n) 
                    for j=state(1,i)+1 : state(1,i_n) -1
                        k = k+dis_stp(j);
                    end
                    hw(i) = nxt(i) + k/v_bus;
                elseif state(1,i) == state(1,i_n)
                    if atstop(i) == 1 && atstop(i_n) == 0
                        for j=state(1,i)+1 : state(1,i_n) -1
                            k = k+dis_stp(j);
                        end
                        hw(i) = nxt(i) + k/v_bus;
                    elseif atstop(i) == 0 && atstop(i_n) == 1                     
                        for j=state(1,i_n) : state(1,i)
                            k = k+dis_stp(j);
                        end
                        k = sum(dis_stp) - k;
                        hw(i) = nxt(i) + k/v_bus;
                    else
                        if nxt(i)>nxt(i_n)
                            for j=state(1,i)+1 : state(1,i_n) -1
                                k = k+dis_stp(j);
                            end
                            hw(i) = nxt(i) + k/v_bus;
                        else
                            for j=state(1,i_n) : state(1,i)
                                k = k+dis_stp(j);
                            end
                            k = sum(dis_stp) - k;
                            hw(i) = nxt(i) + k/v_bus;
                        end
                    end
                else
                    for j=state(1,i_n) : state(1,i)
                        k = k+dis_stp(j);
                    end
                    k = sum(dis_stp) - k;
                    hw(i) = nxt(i) + k/v_bus;
                end

            

            elseif atstop(i)==1 && atstop(i_n) == 1
                k = 0;
                if state(1,i) < state(1,i_n) 
                    for j=state(1,i) : state(1,i_n) -1
                        k = k+dis_stp(j);
                    end
                    hw(i) = nxt(i) + k/v_bus;
                elseif state(1,i) == state(1,i_n)
                    if atstop(i) == 1 && atstop(i_n) == 0
                        for j=state(1,i) : state(1,i_n) -1
                            k = k+dis_stp(j);
                        end
                        hw(i) = nxt(i) + k/v_bus;
                    elseif atstop(i) == 0 && atstop(i_n) == 1                     
                        for j=state(1,i_n) : state(1,i)-1
                            k = k+dis_stp(j);
                        end
                        k = sum(dis_stp) - k;
                        hw(i) = nxt(i) + k/v_bus;
                    else
                        if nxt(i)>nxt(i_n)
                            for j=state(1,i) : state(1,i_n) -1
                                k = k+dis_stp(j);
                            end
                            hw(i) = nxt(i) + k/v_bus;
                        else
                            for j=state(1,i_n) : state(1,i)-1
                                k = k+dis_stp(j);
                            end
                            k = sum(dis_stp) - k;
                            hw(i) = nxt(i) + k/v_bus;
                        end
                    end
                else
                    for j=state(1,i_n) : state(1,i)-1
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
        %when the bus leaves the stop i_n.
        if count ==1
            hw(1) = dis_stp(state(1))/v_bus;
        end
end
        