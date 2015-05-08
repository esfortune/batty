function mothwinganime(im, data, snd, ttl, lim, mname)
% mothwinganime(im, data, snd, ttl, lim, mname)
% 'im' is the image structure from frontwingtrack
% 'data' includes the measurements from frontwingtrack
% 'snd' is the sound channel structure from Spike2 
% 'lim' is the frame limits - for example, [800 1800]
% 'mname' is optional. This specifies an filename for movie output, 'filename.avi'

%% Preparations
% Only make a movie if the user asks for it
if nargin > 5;
    writerObj = VideoWriter(mname);
    writerObj.FrameRate = 30;
    open(writerObj);
end;

% Some constants    
Fs = 1/snd.interval;
tim = 1/Fs:1/Fs:snd.length/Fs;
env = lpf(abs(real(hilbert(snd.values))), Fs, [500 5]);

% Use figure 1
figure(1); 

%% The main loop
for i = lim(1)+8:lim(2);

    clf;
    
    % Our image
    axes('position', [0.1875 0.1500 0.7750 0.8150]); cla;
    imshow(im(i).orig);

    hold on;
    % Loop for the 7 previous data points for the wing contrails
    for j = 0:7;
        if ~isempty(find(data.r.up == i-j, 1)); % Downstrokes
            plot(data.r.x(i-(j+1):i-j),data.r.y(i-(j+1):i-j), 'r*-', 'LineWidth', 3);
        else % Upstrokes
            plot(data.r.x(i-(j+1):i-j),data.r.y(i-(j+1):i-j), 'r*-', 'LineWidth', 1);
        end;
        
        if ~isempty(find(data.r.up == i-j, 1)); % Downstrokes
            plot(data.l.x(i-(j+1):i-j),data.l.y(i-(j+1):i-j), 'g*-', 'LineWidth', 3);
        else % Upstrokes
            plot(data.l.x(i-(j+1):i-j),data.l.y(i-(j+1):i-j), 'g*-', 'LineWidth', 1);
        end;    
    end;
    hold off;

% Position / Velocity plot    
    axes('position', [0.05 0.7 0.25 0.15]); cla;

    % position
    plot(1:22, -data.r.y(i-7:i+14) -200, 'r', 1:22, -data.l.y(i-7:i+14) -200, 'g'); 
    xlim([1 22]);    
    ylim([-800 600]);
hold on;
    % velocity
    plot(1:22, -data.r.vel(i-7:i+14) +300, 'r', 1:22, -data.l.vel(i-7:i+14) +300, 'g');

    text(2,-500, 'Position');
    text(2,100, 'Velocity');
    set(gca, 'xtick', [], 'ytick', []);
    
    % Vertical line
    
    plot([8 8],[600 -800], 'k-');
hold off;

% Envelope plot
    axes('position', [0.05 0.5 0.25 0.15]); cla;
    tt = find(tim > ttl.times(i-7) & tim < ttl.times(i+14));
    plot(tim(tt), env(tt));
    xlim([tim(tt(1)) tim(tt(end))]);
    ylim([-0.01 0.035]);
    text(tim(tt(1))+0.002, 0, 'AM');
    set(gca, 'xtick', [], 'ytick', []);
    
% Specgram plot
    axes('position', [0.05 0.3 0.25 0.15]); cla;
    specgram(snd.values(tt), 256, Fs, [], 250);
    ylim([30000 60000]);
    colormap('HOT'); caxis([-50 10]);
    set(gca, 'xtick', [], 'ytick', [], 'xlabel', [], 'ylabel', []);

% capture move for movie, if requested 

if nargin > 5;
    out= getframe(gcf);
    writeVideo(writerObj,out);    
end;

    pause(0.0001);
    %close(1);
    
end;

if nargin > 5;
    close(writerObj);
end;
