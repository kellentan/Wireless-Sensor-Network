% network_loc.m
function pos_free = network_loc(N, E, pos_anchor, rho)
%% Parameter Initialization
beta1 = 0.8; beta2 = 2.0; lambda = 10^-5; constraint = 10^-5; 
L = size(E, 1);
disp(L);

M = N - size(pos_anchor, 1); % number of free points
pos_free = rand(M, 2); % random initialization
x_k = [pos_free(:, 1); pos_free(:, 2)]; % unraveled coordinates
loc_total = [pos_free; pos_anchor]; % all node coordinates
grad = 1 / 0; % begin with infinite gradient magnitude

%% Levenberg-Marquardt Method
k = 1; % iteration number
lambdas = []; % stores lambda values per iteration
costs = []; % stores cost function values per iteration
k_vals = []; % stores iteration numbers per iteration

while grad > constraint  
    % evaluate f(x^(k))
    f_vals = zeros(L, 1);
    f_vals(:, 1) = sqrt( sum((loc_total(E(:, 1), :) ...
        - loc_total(E(:, 2), :)).^2, 2) ) - rho(:, 1);
    
    % compute derivative matrix
    df = zeros(L, 2*M);
    for i = 1:L
        % 1) both nodes known, i.e. entire row is zero
        if (E(i, 1)) > M && (E(i, 2)) > M 
            continue
        % 2) node 1 known
        elseif (E(i, 1)) > M
            uv_i = pos_anchor(E(i, 1) - M, :);
            uv_j = pos_free(E(i, 2), :);
            temp = 1 / (sqrt( sum((uv_i - uv_j).^2) ));
            
            % input derivatives w/ respect to node 2
            df(i, E(i, 2)) = temp * -(uv_i(:, 1) - uv_j(:, 1));
            df(i, E(i, 2) + M) = temp * -(uv_i(:, 2) - uv_j(:, 2));
        % 3) node 2 known
        elseif (E(i, 2)) > M % node 2 is known
            uv_i = pos_free(E(i, 1), :);
            uv_j = pos_anchor(E(i, 2) - M, :);
            temp = 1 / sqrt( sum((uv_i - uv_j).^2) );
            
            % input derivatives w/ respect to node 1
            df(i, E(i, 1)) = temp * (uv_i(:, 1) - uv_j(:, 1));
            df(i, E(i, 1) + M) = temp * (uv_i(:, 2) - uv_j(:, 2));   
        % 4) neither node known
        else
            uv_i = pos_free(E(i, 1), :);
            uv_j = pos_free(E(i, 2), :);
            temp = 1 / sqrt( sum((uv_i - uv_j).^2) );
            
            % input derivatives w/ respect to node 1 and node 2
            df(i, E(i, 1)) = temp * (uv_i(:, 1) - uv_j(:, 1));
            df(i, E(i, 1) + M) = temp * (uv_i(:, 2) - uv_j(:, 2));   
            df(i, E(i, 2)) = temp * -(uv_i(:, 1) - uv_j(:, 1));
            df(i, E(i, 2) + M) = temp * -(uv_i(:, 2) - uv_j(:, 2));  
        end
    end
    
    % compuate iteration gradient and iteration gradient magnitude
    temp = 2 * df' * f_vals; 
    grad_temp = norm(temp, 1);
    
    % compute x_hat using backslash operator
    x_hat = x_k - (df' * df + lambda * eye(2 * M)) \ (temp / 2);
    
    % store the potential new coordinates
    pred_E = [x_hat(1:M, 1), x_hat(1+M:size(x_hat, 1))];
    tmp_loc = [pred_E; pos_anchor]; 
    
    % evaluate f(x_hat)
    f_x_vals = zeros(L, 1);
    f_x_vals(:, 1) = sqrt( sum((tmp_loc(E(:, 1), :) ...
        - tmp_loc(E(:, 2), :)).^2, 2) ) - rho(:, 1);
    
    % store cost function, lambda parameter, and iteration number values
    lambdas(end + 1) = lambda;
    costs(end + 1) = norm(f_vals)^2;
    k_vals(end + 1) = k;
    
    % update parameters on next iteration
    if norm(f_x_vals)^2 < norm(f_vals)^2
        x_k = x_hat;
        lambda = lambda * beta1;
        
        pos_free = [x_k(1:M), x_k(1+M:2*M)];
        loc_total = [pos_free; pos_anchor];
    else
        lambda = lambda * beta2;
    end
    
    % update iteration number and iteration gradient magnitude
    k = k + 1;
    grad = grad_temp;
end

%% Cost Function & Lambda v. Iteration
plot(k_vals, lambdas); title("Lambda v. Iteration Number"); grid on;
xlabel("Iteration Number"); ylabel("Lambda"); hold off;
figure;
plot(k_vals, costs); title("Cost Function v. Iteration Number"); grid on;
xlabel("Iteration Number"); ylabel("Cost Function");
figure;

%% Sensor Coordinates (Utilized in Calling Script)
% scatter(pos(1:N-K, 1), pos(1:N-K, 2), "go", "filled"); hold on;
% scatter(pos_free(:, 1), pos_free(:, 2), "bo");
% scatter(pos_anchor(:, 1), pos_anchor(:, 2), "s", "ko", "filled");
% 
% for i=1:N-K
%     plot([pos_free(i, 1), pos(i, 1)], [pos_free(i, 2), pos(i, 2)], "k");
% end
% 
% title("Sensor Coordinate Estimates (N=100, R=0.4, s=0.05)");
% xlabel("u"); ylabel("v"); grid on;
% xlim([-0.2, 1.2]); ylim([-0.2, 1.2]);