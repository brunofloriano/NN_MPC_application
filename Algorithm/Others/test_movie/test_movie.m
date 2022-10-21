%% Define waveform properties
f = 1;          % Frequency (Hz)
t = 0:.001:2;   % Eval time (s)
Y_Fourier = zeros(size(t)); % Preallocate
 
%% Generate the frames
for k = 1:20;
    Y_Fourier = (4/pi)*sin(2*pi*(2*k-1)*f*t)/(2*k-1) + Y_Fourier;
    if k == 1
        % Only create the plot for 1st iteration, update wave for
        % subsequent
        figure();
        hold('on'); grid('on');
        ylim([-1.5 1.5]);
        set(gca, 'xtick', 0:0.5:2);
        h1 = plot(t, square(t*2*pi), 'LineWidth', 2);
        h2 = plot(t,Y_Fourier, 'LineWidth', 2);
        legend('Square Wave (1 Hz)', 'Fourier Approximation');
        xlabel('Time (s)');
        ylabel('Amplitude');
    else
        % Update y data
        set(h2, 'ydata', Y_Fourier)
    end
    title(['k = ' num2str(k)]);
    % Save as png with a resolution of 150 pixels per inch
    print(['Frame ' num2str(k)], '-dpng', '-r150');
end
 
%% Now time to generate the animated gif
GifName = 'SquareWaveAnimation.gif';
delay = 0.5;    % Delay between frames (s)
for ii = 1:20
    [A, ~] = imread(['Frame ' num2str(ii) '.png']);
    [X, map] = rgb2ind(A, 256);
    if ii == 1
        imwrite(X, map, GifName, 'gif', 'LoopCount', inf, 'DelayTime', delay)
    else
        imwrite(X, map, GifName, 'gif', 'WriteMode', 'append', 'DelayTime', delay)
    end
end