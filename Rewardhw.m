function reward = Rewardhw(bus, state, action, hw, n_st, gamma, beta, hwt)
    
    n_b = size(state,2);
    ipl = iplus(bus,n_b); %bus behind the current bus
    stateb = state(:,ipl); %State is state of all the modules
    state = state(:,bus);
    hw = hw(bus);
    rskip = hw - gamma*hwt;
    rstop = -rskip;
    rsplit = hw - gamma*hwt;
    bnst = iplus(stateb(1),n_st);
    curst = state(1);
    if bnst == curst || stateb(1) == curst
        rjoin = beta*hwt - hw;
    else
        rjoin = -1;
    end
    
    switch action
        case 0
            reward = rstop;
        case 1
            reward = rskip;
        case 2
            reward = rsplit;
        case 3
            reward = rjoin;
    end
 end


