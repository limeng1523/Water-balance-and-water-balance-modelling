function [x_op,f_op] = pso( N, M, c1, c2, w, ntime, P, Qobs, R_glo, T_air, rel_hum, v_wind, p_air )
% N£ºnumber of particle
% c1,c2£ºlearning rate
% w£ºinertia weight
% M£ºmaximum number of iterations
% Initialize the individual of the population

prior = [0.6025 92.91 0.1889 0.0443 0.004 0.00051 5.58 10.52 1.896 23.46];
pmin = [0 50 0.01 0.01 0.0001 0.0001 2 2 2 2];
pmax = [1 500 1 1 1 1 250 250 250 250];
format long;
for i = 1:N
    x(i,:) = abs(randn * prior / 20 + prior);
    v(i,:) = randn * prior / 20;
end

% first calculate the fitness of each particle£¬than initialize pi and pg
for i = 1:N
    p(i) = hbv_nl( ntime, P, Qobs, R_glo, T_air, rel_hum, v_wind, p_air, x(i,:) );
    y(i,:) = x(i,:);
end 
pg = x(N,:);  % pg is global optimal solution
for i = 1:(N-1)
    if hbv_nl( ntime, P, Qobs, R_glo, T_air, rel_hum, v_wind, p_air, x(i,:)) > hbv_nl( ntime, P, Qobs, R_glo, T_air, rel_hum, v_wind, p_air, pg )
       pg = x(i,:);
    end
end

for t = 1:M
    for i = 1:N
        v(i,:) = w * v(i,:) + c1 * rand * (y(i,:) - x(i,:)) + c2 * rand * (pg - x(i,:));
        x(i,:) = min(pmax, max(pmin, (x(i,:) + v(i,:))));
        if hbv_nl( ntime, P, Qobs, R_glo, T_air, rel_hum, v_wind, p_air, x(i,:)) > p(i)
            p(i) = hbv_nl( ntime, P, Qobs, R_glo, T_air, rel_hum, v_wind, p_air, x(i,:));
            y(i,:) = x(i,:);
        end
        if p(i) > hbv_nl( ntime, P, Qobs, R_glo, T_air, rel_hum, v_wind, p_air, pg)
            pg = y(i,:);
        end
    end
    Pbest(t) = hbv_nl( ntime, P, Qobs, R_glo, T_air, rel_hum, v_wind, p_air, pg);
end
plot(1:M, Pbest);
x_op = pg;
f_op = hbv_nl( ntime, P, Qobs, R_glo, T_air, rel_hum, v_wind, p_air, pg);
end

