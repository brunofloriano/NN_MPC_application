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