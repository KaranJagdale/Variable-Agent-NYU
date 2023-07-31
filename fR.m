function f_r = fR(state,im)
%function to decide if the given module is front or the rear modules of
%the bus
    imn = im;
    for i=1:im-1
        if state(3,i) == 1
            imn = imn + 1;
        end
    end
    if mod(imn,2) == 1
        f_r = 1;
    else
        f_r = 0;
    end
end