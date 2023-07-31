function lsum = lambsum(st_num1, st_num2, arr_par, n_s)
    %gives the sum of arrival parameters between two buses whose bus stop
    %numbers are given st_num1 and st_num2
   
    if st_num1 <= st_num2 
        if st_num1 ~= n_s
            lsum = sum(arr_par((st_num1+1):st_num2));
        else
            lsum = 0;
        end

    else
        if st_num1 ~= n_s
            lsum = sum(arr_par(st_num1+1 : n_s));
            lsum = lsum + sum(arr_par(1:st_num2));
        else
            lsum = sum(arr_par(1:st_num2));
        end
    end
end