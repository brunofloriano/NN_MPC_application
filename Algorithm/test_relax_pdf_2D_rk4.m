function test_relax_pdf_2D_rk4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% begin user input %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
G.dx=0.1; G.dt=.001; G.xh=G.dx/2; G.d=2; G.dif = 1; 
dist1.mean = [0;0]; dist1.std = [1;1]; dist1.P = [dist1.std(1,1)^2 0; 0 dist1.std(2,1)^2]; dist1.n = round(((6*dist1.std(1,1)/G.dx)+1)*((6*dist1.std(2,1)/G.dx)+1)); 
dist1.Pinv = inv(dist1.P); dist1.Pdet = det(dist1.P); 
dist2.mean = [-3;-3]; dist2.std = [0.2;0.2]; dist2.P = [dist2.std(1,1)^2 0; 0 dist2.std(2,1)^2]; dist2.n = round(((6*dist2.std(1,1)/G.dx)+1)*((6*dist2.std(2,1)/G.dx)+1));
dist2.Pinv = inv(dist2.P); dist2.Pdet = det(dist2.P); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of user input %%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Hole in the middle Distribution
D_t = Initialize_PDF(G,dist1); D_t.P = D_t.P-D_t.P(2);
D = Initialize_PDF(G,dist1); D = Initialize_vuw(D,G); 
D_ag = Initialize_PDF(G,dist2); %D = combine_PDF(D,D_ag,G);

% for i=1:(2*3/G.dx)
%     dist2.mean = dist2.mean + [G.dx; G.dx]; 
%     D_ag = Initialize_PDF(G,dist2); 
%     D = combine_PDF(D,D_ag,G); 
% end
D = boundary_conditions(D,G); plot_initial(D,G);

K = RHS_P(D,G); 
%D_new.P = dictionary(); 
%D_new.P=D.P+G.dt.*K;
%diff = sum(abs(D.P - D_new.P)); 
diff = 1e10;
timestep = 1; clear F_top; clear F_dtv; %clear F_side;  clear F_iso; 
while(diff > 0.02)
  k1=RHS_P(D,G); D.P = D.P + (G.dt/2).*k1;
  k2=RHS_P(D,G); D.P = D.P + (G.dt/2).*k2;
  k3=RHS_P(D,G); D.P = D.P + G.dt.*k3;
  k4=RHS_P(D,G);
  D_new.P=D.P+(G.dt/6).*k1+(G.dt/3).*(k2+k3)+(G.dt/6).*k4;
  [diff, F1, F2] = plot_PDF(D,G,timestep,dist1,D_t,D_new); 
  F_top(timestep) = F1; F_dtv(timestep) = F2;
  timestep=timestep+1; D.P = D_new.P;  
end 

create_video(F_top, 'relaxing_gauss_top_2D_2nd.mp4');
create_video(F_dtv, 'relaxing_gauss_dtv_2D_2nd.mp4');
%create_video(F_side, 'relaxing_gauss_side')
%create_video(F_iso, 'relaxing_gauss_iso')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = Initialize_PDF(G,dist)
    D.P = zeros(dist.n+1,1); D.P(1) = 0; D.key = zeros(dist.n+1,1); D.key(1) = -1;
    D.v = D.P; D.f = D.P; D.u = D.P; D.w = D.P; l = 2;  

    count_i = 0;
    for i=round((dist.mean(1,1)-3*dist.std(1,1))/G.dx):round((dist.mean(1,1)+3*dist.std(1,1))/G.dx)
        count_i = count_i + 1;
        count_j = 0;
        for j=round((dist.mean(2,1)-3*dist.std(2,1))/G.dx):round((dist.mean(2,1)+3*dist.std(2,1))/G.dx)
            count_j = count_j + 1;
            state = [i;j]; D.key(l,1) = state_conversion(state, G);  
            %D.P(l,1) = ((2*pi)^G.d*dist.Pdet)^(-1/2)*exp((-1/2)*(state.*G.dx-dist.mean)'*dist.Pinv*(state.*G.dx-dist.mean));
            %D.P(l,1) = ((2*pi)^G.d*dist.Pdet)^(-1/2)*exp((-1/2)*norm(state.*G.dx-dist.mean)^4);
%             if norm(state.*G.dx-dist.mean) < 2
%                 D.P(l,1) = 1;
%             else
%                 D.P(l,1) = 0;
%             end
            pdf = (2*pi)^(-2/2)*det(dist.std(1,1)*eye(2))^(-1/2)*exp( -( norm(state.*G.dx-dist.mean)^4 )/(2*dist.std(1,1)^4) );
            p = initial_test;
            D.P(l,1) = p(count_i,count_j);
            
            l=l+1;
        end
    end
    prob_sum = (G.dx)^G.d.*sum(D.P); D.P = D.P/prob_sum; D.n = length(D.P); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = boundary_conditions(D,G)
    for l=2:D.n
        key = D.key(l,1); state = key_conversion(key,G); pos = G.dx.*state; 
        if(max(abs(pos))==3)
            D.P(l,1) = 0;
        end
    end     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = Initialize_vuw(D,G)
    for l=2:D.n
        current_key = D.key(l,1); current_state = double(key_conversion(current_key,G));
           
        x=G.dx.*current_state; 
        
        x_state = x; x_state(1) = x_state(1)+G.xh; 
        r = norm(x_state); theta = atan2(x_state(2),x_state(1)); 
        v1 = -G.dif*r^3*cos(theta); % Potential Bowl
        D.v(l,1) = v1;
        D.u(l,1)=min(v1, 0);
        D.w(l,1)=max(v1, 0);

        y_state = x; y_state(2) = y_state(2)+G.xh; 
        r = norm(y_state); theta = atan2(y_state(2),y_state(1)); 
        v2 = -G.dif*r^3*sin(theta);
        D.v(l,2) = v2;
        D.u(l,2)=min(v2, 0);
        D.w(l,2)=max(v2, 0);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_initial(D,G)
    
    figure(1); clf; grid on; hold on; view(30,30);
    xlabel('x', 'Interpreter', 'Latex')
    ylabel('y', 'Interpreter', 'Latex')
    zlabel('Probability', 'Interpreter', 'Latex')
    xlim([-3, 3])
    ylim([-3, 3])
 
    x = zeros(1,D.n-1); y = x; z = x;
    
    for l=2:D.n
        current_key = D.key(l,1); state = key_conversion(current_key,G); 
        x(1,l-1) = G.dx*state(1); y(1,l-1) = G.dx*state(2); z(1,l-1) = D.P(l,1);
    end
    
    xv = linspace(min(x), max(x), sqrt(D.n-1));
    yv = linspace(min(y), max(y), sqrt(D.n-1));
    [X,Y] = meshgrid(xv, yv);
    Z = griddata(x,y,z,X,Y);
    s = surf(X,Y,Z); 
    s.EdgeColor = 'none';
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [diff, F1, F2] = plot_PDF(D,G,timestep,dist,D_t,D_new)
    
    x = zeros(1,dist.n-1); y = x; z = x; z_t = x; diff = 0; 
    
    for l=2:D.n
        current_key = D.key(l,1); state = key_conversion(current_key,G); 
        x(1,l-1) = G.dx*state(1); y(1,l-1) = G.dx*state(2); z(1,l-1) = D.P(l,1);
        z_t(1,l-1) = D_t.P(l,1);
    end
    
    diff = sum(abs(D.P - D_new.P)); 
    xv = linspace(min(x), max(x), sqrt(D.n-1));
    yv = linspace(min(y), max(y), sqrt(D.n-1));
    [X,Y] = meshgrid(xv, yv);
    Z = griddata(x,y,z,X,Y);
    Z_t = griddata(x,y,z_t,X,Y); 
    
    figure(2); clf; grid on; hold on; view(0,90); %Top View
    title(['iter = ', num2str(timestep), ', t = ', num2str(timestep*G.dt), ', \Delta t = ', num2str(G.dt), ', \Delta', char(949), ' = ', num2str(diff)]);
    xlabel('x', 'Interpreter', 'Latex')
    ylabel('y', 'Interpreter', 'Latex')
    zlabel('Probability', 'Interpreter', 'Latex')
    xlim([-3, 3])
    ylim([-3, 3])
    zlim([0, max(z_t*1.5)])
    surf(X,Y,Z, 'EdgeColor','none'); 
    %surf(X,Y,Z_t, 'FaceAlpha', 0.25, 'EdgeColor','none', 'CData',zeros(sqrt(D.n-1),sqrt(D.n-1)));
    F1 = getframe(gcf);
    drawnow


    figure(3); clf; grid on; hold on; view(-45,0); %Down The Valley View
    title(['iter = ', num2str(timestep), ', t = ', num2str(timestep*G.dt), ', \Delta t = ', num2str(G.dt), ', \Delta', char(949), ' = ', num2str(diff)]);
    xlabel('x', 'Interpreter', 'Latex')
    ylabel('y', 'Interpreter', 'Latex')
    zlabel('Probability', 'Interpreter', 'Latex')
    xlim([-3, 3])
    ylim([-3, 3])
    zlim([0, max(z_t*1.5)])
    surf(X,Y,Z, 'EdgeColor','none'); 
    %surf(X,Y,Z_t, 'FaceAlpha', 0.25, 'EdgeColor','none', 'CData',zeros(sqrt(D.n-1),sqrt(D.n-1)));
    F2 = getframe(gcf);
    drawnow

    
    %{
    figure(4); clf; grid on; hold on; view(45,0); %Side View
    title(['iter = ', num2str(timestep), ', t = ', num2str(timestep*G.dt), ', \Delta t = ', num2str(G.dt), ', \Delta', char(949), ' = ', num2str(diff)]);
    xlabel('x', 'Interpreter', 'Latex')
    ylabel('y', 'Interpreter', 'Latex')
    zlabel('Probability', 'Interpreter', 'Latex')
    xlim([-3, 3])
    ylim([-3, 3])
    zlim([0, max(z_t*1.5)])
    surf(X,Y,Z, 'EdgeColor','none'); 
    surf(X,Y,Z_t, 'FaceAlpha', 0.25, 'EdgeColor','none', 'CData',zeros(sqrt(D.n-1),sqrt(D.n-1)));
    if(timestep~=1)
        if (mod(timestep,5)==0)
            F3 = getframe(gcf);
        else
            F3 = 0; 
        end
    else
        F3 = getframe(gcf);
    end
    drawnow

    figure(5); clf; grid on; hold on; view(30,30); %Iso View
    title(['iter = ', num2str(timestep), ', t = ', num2str(timestep*G.dt), ', \Delta t = ', num2str(G.dt), ', \Delta', char(949), ' = ', num2str(diff)]);
    xlabel('x', 'Interpreter', 'Latex')
    ylabel('y', 'Interpreter', 'Latex')
    zlabel('Probability', 'Interpreter', 'Latex')
    xlim([-3, 3])
    ylim([-3, 3])
    zlim([0, max(z_t*1.5)])
    surf(X,Y,Z, 'EdgeColor','none'); 
    surf(X,Y,Z_t, 'FaceAlpha', 0.25, 'EdgeColor','none', 'CData',zeros(sqrt(D.n-1),sqrt(D.n-1)));
    if(timestep~=1)
        if (mod(timestep,5)==0)
            F4 = getframe(gcf);
        else
            F4 = 0; 
        end
    else
        F4 = getframe(gcf);
    end
    drawnow
    %}
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D_p = combine_PDF(D_p, D_ag,G)
    for l=2:D_p.n
        current_key = D_p.key(l,1); ag_pos = find(D_ag.key==current_key); 
        if(~isempty(ag_pos))
            D_p.P(l,1) = max(D_p.P(l,1) - D_ag.P(ag_pos,1),0);
        end
    end 
    prob_sum = (G.dx)^G.d.*sum(D_p.P); D_p.P(3:D_p.n-1) = D_p.P(3:D_p.n-1)/prob_sum;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = RHS_P(D,G)
    diffusion = zeros(D.n,1); K = zeros(D.n,1);

    for l=2:D.n %Flux out
        l_key = D.key(l,1); l_state = key_conversion(l_key,G); 
        for q=1:G.d
            k_state = l_state; k_state(q) = k_state(q)+1; k_key = state_conversion(k_state,G);
            i_state = l_state; i_state(q) = i_state(q)-1; i_key = state_conversion(i_state,G);
            k_pos = find(D.key==k_key); i_pos = find(D.key==i_key); 
            if(isempty(k_pos)),k_pos = 1;end
            if(isempty(i_pos)),i_pos = 1;end
            D.f(l,q) = D.u(l,q)*D.P(k_pos,1) + D.w(l,q)*D.P(l,1);
            diffusion(l,1) = diffusion(l,1) + D.P(k_pos,1)-2*D.P(l,1)+D.P(i_pos,1); 
        end
    end
    
    for d=1:G.d %Flux out w/ CTU
        for l=2:D.n
            l_key = D.key(l,1); l_state = key_conversion(l_key,G); 
            i_state = l_state; i_state(q) = i_state(q)-1; i_key = state_conversion(i_state,G);
            i_pos = find(D.key==i_key); if(isempty(i_pos)),i_pos = 1;end

            F = G.dt*(D.P(l,1)-D.P(i_pos,1))/(2*G.dx);
            for e=1:G.d
                if e~=d
                    j_state = l_state; j_state(e) = j_state(e)-1; j_key = state_conversion(j_state,G);
                    j_pos = find(D.key==j_key); if(isempty(j_pos)),j_pos = 1;end
                    p_state = i_state; p_state(e) = p_state(e)-1; p_key = state_conversion(p_state,G);
                    p_pos = find(D.key==p_key); if(isempty(p_pos)),p_pos = 1;end

                    D.f(l,e) = D.f(l,e)-D.w(l,e)*D.w(i_pos,d)*F;
                    D.f(j_pos,e) = D.f(j_pos,e)-D.u(j_pos,e)*D.w(i_pos,d)*F;
                    D.f(i_pos,e) = D.f(i_pos,e)-D.w(i_pos,e)*D.u(i_pos,d)*F;
                    D.f(p_pos,e) = D.f(p_pos,e)-D.u(p_pos,e)*D.u(i_pos,d)*F;
                end
            end

            if(D.v(i_pos,d) > 0)
                i_i_state = i_state; i_i_state(d) = i_i_state(d)-1; i_i_key = state_conversion(i_i_state,G);
                i_i_pos = find(D.key==i_i_key); if(isempty(i_i_pos)),i_i_pos = 1;end
                th = (D.P(i_pos,1) - D.P(i_i_pos,1))/(D.P(l,1)-D.P(i_pos,1));
            else
                k_state = l_state; k_state(q) = k_state(q)+1; k_key = state_conversion(k_state,G);
                k_pos = find(D.key==k_key); if(isempty(k_pos)),k_pos = 1;end
                th = (D.P(k_pos,1) - D.P(l,1))/(D.P(l,1)-D.P(i_pos,1));
            end

            t = abs(D.v(i_pos,d));
            D.f(i_pos,d) = D.f(i_pos,d) + t*(G.dx/G.dt-t)*F*MC(th);
        end
    end
          
    for l=2:D.n
        l_key = D.key(l,1); l_state = key_conversion(l_key,G); l_pos = abs(G.dx.*l_state);
        if (max(l_pos)~=3)
            K(l,1) = K(l,1) + (1/G.dx^2)*(G.dif*diffusion(l,1));
            for d=1:G.d
                i_state = l_state; i_state(d) = i_state(d)-1; i_key = state_conversion(i_state,G);
                i_pos = find(D.key==i_key); 
                if(isempty(i_pos))
                    K(l,1) = K(l,1)-(1/G.dx)*(D.f(l,d))+(1/G.dx^2)*(G.dif*diffusion(l,1));
                else
                    K(l,1) = K(l,1)-(1/G.dx)*(D.f(l,d)-D.f(i_pos,d));
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phi]=MC(th), phi=max(0,min([(1+th)/2 2 2*th]));   end   % Flux limiters
function [phi]=VL(th), phi=min((th+abs(th))/(1+abs(th)),0); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = state_conversion(state,G)
    shift_state = ShiftState(state,G.d);
    if (G.d==1)
        key = shift_state(1);
    else
        key = CantorPair(shift_state);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = key_conversion(key,G)
    if(G.d==1)
        shift_state = [key]; 
    else
        shift_state = CantorUnpair(key,G.d);
    end
    state = UnshiftState(shift_state,G.d);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function shift_state = ShiftState(state,d)
    shift_state = zeros(d,1);
    for i=1:d
        if(state(i,1)<0)
            shift_state(i,1)=-2*state(i,1)-1;
        else
            shift_state(i,1)=2*state(i,1);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = UnshiftState(shift_state,d)
    state = zeros(d,1);
    for i=1:d
        if(mod(shift_state(i,1),2)==0)
            state(i,1)=shift_state(i,1)/2;
        else
            state(i,1)=(shift_state(i,1)+1)/-2;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = CantorPair(state)
    if(length(state)>2)
        last = state(end,1); state(end,1) = [];
        x = CantorPair(state); y = last;
        key = (1/2)*(x+y)*(x+y+1)+y;
    else
        x=state(1,1); y=state(2,1);
        key = (1/2)*(x+y)*(x+y+1)+y;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = CantorUnpair(key,d)
    state = zeros(d,1);
    for i=2:d
        z=key; w=floor(((8*z+1)^(1/2)-1)/2); t=(w^2+w)/2; 
        y=z-t; x=w-y; state(d-i+1,1)=x; state(d-i+2,1) =y;
        key=x;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function create_video(F,title)
    
    
    writerObj = VideoWriter(title, 'MPEG-4');
    writerObj.FrameRate = 10;
    open(writerObj);
    for i=1:length(F)
        frame = F(i);
        writeVideo(writerObj, frame);
    end

    close(writerObj);
end