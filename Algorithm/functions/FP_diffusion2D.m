%% Fokker-Planck drift and diffusion 2D
% This function uses the Fokker-Planck equation to compute the time
% eovlution of a single variable PDF

function p_new = FP_diffusion2D(p,D,u,delta_t,delta_x,p_initial,boundary)

epsilon = 1e-6;

n = size(p,2);  % size of X dimension
m = size(p,1);  % size of Y dimension

% Initialize the PDF of the next time step
%p_new = zeros(n,1);

% Initialize the PDF shifted in space
% p_minus = [0;p(1:n-1)];
% p_plus = [p(2:n);0];

% First order derivative term
%A = u.*p;
Axv = u(:,:,1).*p;
Ayv = u(:,:,2).*p;

% Second order derivative term
B = D.*p;

%% Shift in space
Axplus = [Axv(:,2:n) zeros(m,1)];
%Ayplus = [zeros(1,n); Ayv(1:m-1,:)];
Ayplus = [Ayv(2:m,:); zeros(1,n)];
%Axminus = [zeros(m,1) Axv(:,1:n-1)];

Bxplus = [B(:,2:n) zeros(m,1)];
%Byplus = [zeros(1,n); B(1:m-1,:)];
Byplus = [B(2:m,:); zeros(1,n)];

Bxminus = [zeros(m,1) B(:,1:n-1)];
%Byminus = [B(2:m,:); zeros(1,n)];
Byminus = [zeros(1,n); B(1:m-1,:)];

%Bxyplus = [Byplus(:,2:n) zeros(m,1)];
Bxy11 = [Byplus(:,2:n) zeros(m,1)];
Bxy01 = [zeros(m,1) Byplus(:,1:n-1)];
Bxy10 = [Byminus(:,2:n) zeros(m,1)];
Bxy00 = [zeros(m,1) Byminus(:,1:n-1)];

%% Shift in space but not with zeros
% Axplus = [Axv(:,2:n) Axv(:,n)];
% Ayplus = [Ayv(1,:); Ayv(1:m-1,:)];
%
% Bxplus = [B(:,2:n) B(:,n)];
% Byplus = [B(1,:); B(1:m-1,:)];
%
% Bxminus = [B(:,1) B(:,1:n-1)];
% Byminus = [B(2:m,:); B(m,:)];
%
% Bxyplus = [Byplus(:,2:n) Byplus(:,n)];

%% Shift in space but not with zeros but smart
% Axplus = [Axv(:,2:n) (2*Axv(:,n) - Axv(:,n-1))];
% Ayplus = [(2*Ayv(1,:)-Ayv(2,:)); Ayv(1:m-1,:)];
% 
% Bxplus = [B(:,2:n) (2*B(:,n)-B(:,n-1))];
% Byplus = [(2*B(1,:) - B(2,:)); B(1:m-1,:)];
% 
% Bxminus = [(2*B(:,1)-B(:,2)) B(:,1:n-1)];
% Byminus = [B(2:m,:); (2*B(m,:)-B(m-1,:))];
% 
% Bxyplus = [Byplus(:,2:n) (2*Byplus(:,n)-Byplus(:,n-1))];


% Compute derivates
Ax = (Axplus - Axv)/delta_x;
Ay = (Ayplus - Ayv)/delta_x;
Bx2 = (Bxplus - 2*B + Bxminus)/delta_x^2;
By2 = (Byplus - 2*B + Byminus)/delta_x^2;
%Bxy = (Bxyplus - Bxplus - Byplus + B)/delta_x^2;
Bxy = (Bxy11 - Bxy01 - Bxy10 + Bxy00)/(4*delta_x^2);

% Compute time step
increment = (-(Ax + Ay) + Bx2 + 2*Bxy + By2);
%increment(:,1) = zeros(n,1);
%increment(1,:) = zeros(1,n);
% for i = 1:n
%     for j = 1:m
%         x = (i - n/2)*delta_x;
%         y = (j - m/2)*delta_x;
%         z = norm([x y]);
%         %if abs(increment(i,j)) > epsilon && (abs(i-n/2) > n/2 - 10 || abs(j-m/2) > m/2 - 10)
%         if true %inside(i,j) == 0
%             %increment(i,j) = p_initial(i,j) - p(i,j);
%             %if abs(p_initial(i,j) - p(i,j)) < epsilon
%             %if norm([i-n/2 j-m/2]) > n/2
%             %    increment(i,j) = 0;
%             %increment(i,j) = increment(i,j)*exp(-(z)^2/(10^2)); %*(1-inside(i,j));
%             %increment(i,j) = increment(i,j)*(1 - 1./(1+exp(-(z-9)))); %*(1-inside(i,j));
%         else
%             increment(i,j) = increment(i,j)*(inside(i,j));
%         end
% 
%     end
% end
%increment = increment.*boundary;
%idx = find(abs(increment) > epsilon;
%increment(idx) = 0;
p_new = p + delta_t*increment;
%p_new = p_new(2:n-1,2:m-1);
