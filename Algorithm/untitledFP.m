% Define parameters
D = 1; % Diffusion coefficient
V = @(x) -x.^3; % Potential function
T = 1; % Total time
N = 100; % Number of time steps
dt = T/N; % Time step
x_min = -5; % Minimum x
x_max = 5; % Maximum x
Nx = 100; % Number of spatial grid points
dx = (x_max-x_min)/Nx; % Spatial step
x = linspace(x_min+dx,x_max-dx,Nx)'; % Spatial grid

% Construct the discretized Fokker-Planck matrix
A = sparse(Nx,Nx);
A = A + (1/dx^2)*spdiags([-2*ones(Nx-1,1),ones(Nx-1,1),ones(Nx-1,1)],[-1,0,1],Nx-1,Nx-1);
for i=1:Nx-1
    A(i,i) = A(i,i) - V(x(i))/D;
end

% Compute the probability density function over time
P = zeros(Nx,N);
P(:,1) = exp(-V(x)./D); % Initial condition
for n=2:N
    P(:,n) = expm(A*dt).*P(:,n-1);
end

% Plot the probability density function over time
figure;
for n=1:10:N
    plot(x,P(:,n),'LineWidth',2);
    xlim([x_min x_max]);
    ylim([0 max(P(:,1))]);
    xlabel('x');
    ylabel('P(x,t)');
    title(['Fokker-Planck equation, t=' num2str((n-1)*dt)]);
    drawnow;
end