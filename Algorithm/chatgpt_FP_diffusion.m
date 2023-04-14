% Fokker-Planck equation solver
% Solve ∂p(x,t)/∂t = - ∂/∂x (A(x,t)p(x,t)) + (1/2) ∂²/∂x² (B(x,t)p(x,t))
% using a finite difference method

% Set parameters
L = 10;         % length of domain
T = 1;          % final time
N = 1000;       % number of grid points
M = 1000;       % number of time steps
dx = L/N;       % spatial step size
dt = T/M;       % time step size
x = linspace(-L/2,L/2,N);  % spatial grid

% Set initial condition
p0 = exp(-x.^2/2);

% Set drift and diffusion coefficients
A = -x;
B = ones(size(x));

% Construct finite difference matrices
D1 = (diag(-ones(N-1,1),-1) + diag(ones(N-1,1),1))/(2*dx); % central difference operator
D2 = (diag(-2*ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1))/(dx^2); % central difference operator

% Compute matrix for time step
L1 = eye(N) - dt*D1*A;
L2 = eye(N) + dt*(1/2)*D2*B;

% Time stepping loop
p = p0;
for i = 1:M
    p = L2\(L1*p);  % solve linear system using backslash operator
end

% Plot results
plot(x,p);
xlabel('x');
ylabel('p(x,t)');