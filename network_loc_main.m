% network_loc_main.m
%% Problem Initialization
N = 100; R = 0.4; s = 0.05;
[E, pos, K] = network_loc_data(N, R);
pos_anchor = pos(N-K+1:N, :);
L = size(E, 1); 
d = sqrt(sum( (pos(E(:, 1), :) - pos(E(:, 2), :)).^2, 2));
rho = (1 + s*randn(L, 1)) .* d;
pos_free = network_loc(N, E, pos_anchor, rho); 

%% Levenberg-Marquardt Visualization
scatter(pos(1:N-K, 1), pos(1:N-K, 2), "go", "filled"); hold on;
scatter(pos_free(:, 1), pos_free(:, 2), "bo");
scatter(pos_anchor(:, 1), pos_anchor(:, 2), "s", "ko", "filled");

for i=1:N-K
    plot([pos_free(i, 1), pos(i, 1)], [pos_free(i, 2), pos(i, 2)], "k");
end

title("Sensor Coordinate Estimates (N=100, R=0.4, s=0.05)");
xlabel("u"); ylabel("v"); grid on;
xlim([-0.2, 1.2]); ylim([-0.2, 1.2]);
