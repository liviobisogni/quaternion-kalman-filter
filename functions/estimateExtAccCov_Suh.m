function [Q_hat_a_b, lambda, mu] = estimateExtAccCov_Suh(r_a, lambda, mu, H_a, P__next_prev, R_a)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [Eq. 34 - 35 Suh]
% ESTIMATEEXTACCCOV_SUH Adaptive estimation algorithm by Y. S. Suh for estimating external acceleration
% covariance matrix Q__a_b.
%
% INPUT:
%   * r_a,              Residual in the accelerometer measurement update        (3 x M_1) vector        [m / s^2]
%   * lambda,           Threshold between first and second condition            (3 x (M_2+1)) matrix
%   * mu,               Defined after [Eq. 32 Suh]                              (3 x (M_2+1)) matrix
%   * H_a,              Accelerometer measurement matrix                        (3 x 9) matrix
%   * P__next_prev,     State covariance matrix                                 (9 x 9) matrix
%   * R_a               Covariance measurement matrix of the accelerometer      (3 x 3) matrix
%
% OUTPUT:
%   * Q_hat_a_b,        Estimated external acceleration covariance              (3 x 3) matrix
%   * lambda,           (newly computed) lambda                                 (3 x (M_2+1)) matrix
%   * mu,               (newly computed) mu                                     (3 x (M_2+1)) matrix
%
% Notation: k = 1 is the most recent sample, whereas k = 3 is the oldest one
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    global M_1 M_2;
    
    global estimateExtAccCov_COUNTER    % TEMP (DEBUG ONLY)!!!!!!!!!!!!!!

    % Check number of arguments
    narginchk(6,6);
    
    if (~isequal(size(r_a), [3 M_1]))
        error('r_a must be a (3 x M_1) vector.');
    end
    if (~isequal(size(lambda), [3 (M_2 + 1)]))
        error('lambda must be a (3 x (M_2 + 1)) vector.');
    end
    if (~isequal(size(mu), [3 (M_2 + 1)]))
        error('mu must be a (3 x (M_2 + 1)) vector.');
    end
        
    % [Eq. 30 Suh]
    % U: a symmetric (3 x 3) matrix
    temp = 0;
    for i = 1:M_1
        temp = temp + r_a(:,i) * r_a(:,i)';     % N.B.: U = U(k), i.e., it is the most recent U matrix
    end
    U = temp / M_1;

    lambda = circshift(lambda,1,2); % shift each row of lambda 1 position to the right (i.e., age lambda's of 1 sample)
    mu = circshift(mu,1,2);         % shift each row of mu 1 position to the right (i.e., age mu's of 1 sample)
    
    % Returns diagonal matrix D of eigenvalues and matrix V whose columns
    % are the corresponding right eigenvectors, so that: U * V = V * D.
    [V,D] = eig(U);
    
    % (M_2 + 1) most recent eugenvalues
    lambda(:,1) = diag(D);
    
    % (M_2 + 1) most recent eugenvectors
    u = zeros(M_2 + 1);
    for i = 1:(M_2 + 1)                     % redundant (you could just type u = V)
        u(:,i) = V(:,i);                    % u(:,i,1) 
    end

%     %% VERIFY THIS - ONLY - BY "HAND"     <--- Okie, verified
%     U2 = zeros(length(r_a));
%     for i = 1:3
%         U2 = U2 + lambda(i,1) * u(:,i,1) * u(:,i,1)';
%     end
%     U - U2                    % NB: should be ideally 0

    % [Eq. post-32 Suh]
    for i = 1:3 % NB: k = 1 is the most recent sample, whereas k = (M_2 + 1) is the oldest one
        mu(i,1) = u(:,i)' * (H_a * P__next_prev * H_a' + R_a) * u(:,i);
    end
    
    % Compute acceleration mode (it is able to detect external acceleration)
    mode = computeAccMode(lambda, mu);
    
    % [Eq. 34 Suh]
    if (mode == 1)                      % 'Mode 1' aka 'No external acceleration Mode'
        Q_hat_a_b = zeros(3);
        
    % [Eq. 35 Suh]
    elseif (mode == 2)                  % 'Mode 2' aka 'External acceleration Mode'
        estimateExtAccCov_COUNTER = estimateExtAccCov_COUNTER + 1;    % TEMP (DEBUG ONLY)!!!!!!!!!!!!!!

        Q_hat_a_b = zeros(3);
        for i = 1:3
            Q_hat_a_b = Q_hat_a_b + max(lambda(i,1) - mu(i,1), 0) * u(:,i) * u(:,i)';
        end
        
%         Q_hat_a_b = Q_hat_a_b / 100;        % CANC!!!!
%         Q_hat_a_b_SUH = Q_hat_a_b % DEBUG ONLY
    else
        error('Invalid acceleration mode');
    end

end