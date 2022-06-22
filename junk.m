%         headway = zeros(1,n_b);
%        
%         for i=2:n_b
%             
%             if atstop(i)==0 || atstop(i+1) == 0
%                 k = 0;
%                 for j=state(1,i)+1 : state(1,i+1) +1
%                     k = k+dis_stp(j);
%                 end
%                 headway(i) = nxt(i) + k/v_bus - nxt(i+1);
%             end
% 
%             if atstop(i)==1 || atstop(i+1) == 0
%                 k = 0;
%                 for j=state(1,i) : state(1,i+1) +1
%                     k = k+dis_stp(j);
%                 end
%                 headway(i) = nxt(i) + k/v_bus - nxt(i+1);
%             end
% 
%             if atstop(i)==0 || atstop(i+1) == 1
%                 k = 0;
%                 for j=state(1,i)+1 : state(1,i+1)
%                     k = k+dis_stp(j);
%                 end
%                 headway(i) = nxt(i) + k/v_bus;
%             end
% 
%             if atstop(i)==1 || atstop(i+1) == 1
%                 k = 0;
%                 for j=state(1,i)+1 : state(1,i+1) +1
%                     k = k+dis_stp(j);
%                 end
%                 headway(i) = nxt(i) + k/v_bus;
%             end
% 
%             
%         end
%         
%     
%         %It is assumed that the passengers come at the bus stop for the first
%         %time in the day more or less at the same time. The busses are planned
%         %such that for the first bus, the passengers start coming to stop i
%         %when the bus leaves the stop i-1.
%         if count ==1
%             headway(1) = disp_stp(state(1))/v_bus;
%         else
%             if atstop(1)==0 || atstop(n_b) == 0
%                 k = 0;
%                 for j=state(1,1)+1 : state(1,n_b) +1
%                     k = k+dis_stp(j);
%                 end
%                 headway(1) = nxt(1) + k/v_bus - nxt(n_b);
%             end
% 
%             if atstop(1)==1 || atstop(n_b) == 0
%                 k = 0;
%                 for j=state(1,1) : state(1,n_b) +1
%                     k = k+dis_stp(j);
%                 end
%                 headway(1) = nxt(1) + k/v_bus - nxt(n_b);
%             end
% 
%             if atstop(1)==0 || atstop(n_b) == 1
%                 k = 0;
%                 for j=state(1,1)+1 : state(1,n_b)
%                     k = k+dis_stp(j);
%                 end
%                 headway(1) = nxt(1) + k/v_bus;
%             end
% 
%             if atstop(1)==1 || atstop(n_b) == 1
%                 k = 0;
%                 for j=state(1,1) : state(1,n_b) 
%                     k = k+dis_stp(j);
%                 end
%                 headway(1) = nxt(1) + k/v_bus;
%             end
% 
%         end


        
    
%         switch max(r_st,r_sk,r_sp,r_jn)
%             case r_st
%                 action = 0;
%             case r_sk
%                 action = 1;
%             case r_sp
%                 action = 2;
%             case r_jn
%                 action = 3;
%         end
        
%if r_jn > 0
%                     n_b = n_b -1;
%                     scale_i = zeros(n_s, n_b);
%                     for i=1:im-2
%                         scale_i(:,i) = scale(:,i);
%                     end
%                     for i=im:n_b
%                         scale_i(:,i) = scale(:,i+1);
%                     end
%                     scale_i(:,im-1) = [state(1,im);state(2,im)+state(2,im-1);1];
%                     scale = scale_i;
%                     %Here, we have joined the imth module to the im-1th module
%                     %and that becomes im-1th module in the new configuration
%                 end