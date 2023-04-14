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