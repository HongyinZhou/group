function [rho] = quasi_level_set(rho, threshold_sum, threshold_mean, threshold_num)
global x_array
global y_array
global lambda_ele_x
global lambda_ele_y
global pixel_size
%scanning every cell
for ii = 1:x_array
    for jj = 1:y_array
        limit_x_min = floor(lambda_ele_x / pixel_size)*(ii - 1) + 1;
        limit_x_max = floor(lambda_ele_x / pixel_size)*ii; 
        limit_y_min = floor(lambda_ele_y / pixel_size)*(jj - 1) + 1;
        limit_y_max = floor(lambda_ele_y / pixel_size)*jj;
        lable_x_min = 0;
        lable_x_max = 0;
        lable_y_min = 0;
        lable_y_max = 0;
        for k =  limit_x_min : limit_x_max 
            if sum(rho(limit_y_min : limit_y_max, k)) > threshold_sum&&...
            mean(rho(limit_y_min : limit_y_max, k))> threshold_mean&&...
            sum(sum(rho(limit_y_min : limit_y_max, k)>0.7))> threshold_num
            lable_x_min = k;
            break;
            end
        end
        for k = lable_x_min : limit_x_max 
            if sum(rho(limit_y_min : limit_y_max, k)) > threshold_sum&&...
            mean(rho(limit_y_min : limit_y_max, k))> threshold_mean&&...
            sum(sum(rho(limit_y_min : limit_y_max, k)>0.7))> threshold_num
            lable_x_max = k;
            break;
            end
        end
        for k = lable_y_min : limit_y_max 
            if sum(rho(k, limit_y_min : limit_y_max)) > threshold_sum&&...
            mean(rho(k, limit_y_min : limit_y_max))> threshold_mean&&...
            sum(sum(rho(k, limit_y_min : limit_y_max)>0.7))> threshold_num
            lable_y_min = k;
            break;
            end
        end
        for k = lable_y_min : limit_y_max 
            if sum(rho(k, limit_y_min : limit_y_max)) > threshold_sum &&...
            mean(rho(k, limit_y_min : limit_y_max))> threshold_mean&&...
            sum(sum(rho(k, limit_y_min : limit_y_max)>0.7))> threshold_num
            lable_y_max = k;
            break;
            end
        end
        rho(lable_y_min : lable_y_max, lable_x_min : lable_x_max ) = 1;
    end
end
end