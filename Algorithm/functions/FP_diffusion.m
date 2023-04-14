%% Fokker-Planck drift and diffusion
% This function uses the Fokker-Planck equation to compute the time
% eovlution of a single variable PDF 

function p_new = FP_diffusion(p,D,u,delta_t,delta_x)

n = length(p);
% Initialize the PDF of the next time step
%p_new = zeros(n,1);

% Initialize the PDF shifted in space
% p_minus = [0;p(1:n-1)];
% p_plus = [p(2:n);0];

% First order derivative term
first_term = u.*p;

% Second order derivative term
second_term = D.*p;

% Shift in space
first_term_plus = [first_term(2:n) 0];
%first_term_minus = [0 first_term(1:n-1)];
second_term_minus = [0 second_term(1:n-1)];
second_term_plus = [second_term(2:n) 0];

% first_term_plus = [first_term(2:n) first_term(end)];
% second_term_minus = [second_term(1) second_term(1:n-1)];
% second_term_plus = [second_term(2:n) second_term(end)];


% Compute derivates
first_derivate = (first_term_plus - first_term)/delta_x;
%first_derivate = (first_term - first_term_minus)/delta_x;
second_derivate = (second_term_plus - 2*second_term + second_term_minus)/delta_x^2;

% Compute time step
p_new = p + delta_t*(-first_derivate + second_derivate);

