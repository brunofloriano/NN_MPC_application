%% Accumulate Laplacian matrices
% This function takes a Laplacian matrix L (current communication graph) and adds 
% it to the cumulative matrix Lc (cumulative graph)

function [Lc] = accumulate_laplacian(L,Lc)
N = length(L);

% Add matrices
Lc = Lc + L;

for i = 1:N
    for j = 1:N
        % Normalize (make non-diagonal = -1)
        if Lc(i,j) < 0
            Lc(i,j) = -1;
        end

        % Routing
        if L(i,j) < 0
            for k = 1:N
                if Lc(j,k) < 0
                    Lc(i,k) = -1;
                end
            end
        end

    end
end

% Build diagonal
a = Lc.*(ones(N) - eye(N));
b = sum(a');
Lc = a - diag(b);
