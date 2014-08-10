function [tt_moment, tt_integral] = parameter_cost_function_prepare(input_t, mu_order, m_index, d, x_lim)


[n_m, n_x] = size(mu_order);

tt_integral = tt_moment_full( input_t, zeros(n_x, 1), m_index, d, x_lim, 0);

mean = zeros(n_x, 1);

for i = 1 : n_x
    
    ei = zeros(1, n_x); ei(i) = 1;
    
    mean(i) = tt_moment_full(input_t, ei, m_index, d, x_lim, 0);
    
end

tt_moment = cell(n_m,1);

for i = 1 : n_m
    
    tt_moment{i} = tt_moment_full( input_t, mu_order(i,:), m_index, d, x_lim, mean);

end
    