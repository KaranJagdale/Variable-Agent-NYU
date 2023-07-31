function all_a = allowed_actions(state, im, atstop)
    %0 - stop, 1 - skip, 2 - split, 3 - join, 4 - nextbs
    %disp(im)
%     disp('im')
%     disp(im)
    %disp(state)
    n_b= size(state,2);
    if atstop(im) == 1 
        implus = iplus(im,n_b);
        if state(3,im) == 0 && state(3,implus)==0
            all_a = '34';
        else
            all_a = '4';
        end
    else
        if state(3,im) == 0
            all_a = '01';
        else
            all_a = '012';
        end
    end
end