%% Fokker-Planck drift and diffusion 2D
% This function uses the Fokker-Planck equation to compute the time
% eovlution of a single variable PDF

function p_new = FP_diffusion2D2(p_old,D,u,delta_t,delta_x,boundary)


n = size(p_old,1);  % size of X dimension
m = size(p_old,2);  % size of Y dimension

% First order derivative term
Axv = u(:,:,1).*p_old;
Ayv = u(:,:,2).*p_old;

% Second order derivative term
%B = D.*p;
B11 = D(1,1).*p_old;
B12 = D(1,2).*p_old;
B21 = D(2,1).*p_old;
B22 = D(2,2).*p_old;

increment = p_old*0;
for i = 1:n
    for j = 1:m
        %if  %i ~= 1 && j ~= 1 && i ~= n && j ~= m
        if boundary(i,j) == 1
            Ax(i,j) = (Axv(i,j) - Axv(i-1,j))/delta_x;
            Ay(i,j) = (Ayv(i,j) - Ayv(i,j-1))/delta_x;
% 
            %             Bx2 = (B(i+1,j) - 2*B(i,j) + B(i-1,j))/delta_x^2;
            %             By2 = (B(i,j+1) - 2*B(i,j) + B(i,j-1))/delta_x^2;
            %
            %             Bxy = (B(i+1,j+1) - B(i+1,j-1) - B(i-1,j+1) + B(i-1,j-1))/(4*delta_x^2);

            %             increment(i,j) = (-(Ax + Ay) + Bx2 + 2*Bxy + By2);


            Bx2(i,j) = (B11(i+1,j) - 2*B11(i,j) + B11(i-1,j))/delta_x^2;
            By2(i,j) = (B22(i,j+1) - 2*B22(i,j) + B22(i,j-1))/delta_x^2;

            Bxy12 = (B12(i+1,j+1) - B12(i+1,j-1) - B12(i-1,j+1) + B12(i-1,j-1))/(4*delta_x^2);
            Bxy21 = (B21(i+1,j+1) - B21(i+1,j-1) - B21(i-1,j+1) + B21(i-1,j-1))/(4*delta_x^2);

            increment(i,j) = (-(Ax(i,j) + Ay(i,j)) + Bx2(i,j) + Bxy12 + Bxy21 + By2(i,j));

        else
            increment(i,j) = 0;
        end

    end
end


p_new = p_old + delta_t*increment;
