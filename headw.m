function hw = headw(stop1, stop2, atstop, dis_stp,v_bus,nxt,nb)
    
    k = 0;
    %If the rear bus is at stop then the iterant starts from stop1 else it starts at stop+1
    %This is not accurate as the time spent at the stop is not accounted
    %but we are starting with this and then encorporate it later
    if atstop ==1
        l = 0;
    else
        l =1;
    end
   
    for j=stop1+l : stop2
        k = k+dis_stp(j);
    end
    hw = nxt1 + k/v_bus - nxt2;
end