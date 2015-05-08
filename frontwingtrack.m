function [im, data] = frontwingtrack(lims)
% [im data] = frontwingtrack([startframe endframe]);
% cd into the directory with the jpeg files
% Run this code with the first frame and last frame. Leave empty if you
% want all frames.
% THE KEY IS TO start with wings high or low in the first frame.
% Change some defaults if necessary:\
%    
% 
%% Read the files from the images and set the defaults

    flist = dir('*.jpg');

    thresh = 0.02; % Threshold level for grayscale. 0 is black, 1 is white.
    boxlim = [80 160]; % Size in pixels of our search box [x y]
    cnt = 5; % The number of white pixels to qualify - 5 seems to work. Higher will reject more noise and lower is more sensitive.

    % If the user does not specify the limits, do the whole list of jpegs.
    if nargin == 0; % zero 'arguments'
        lims(1) = 2; % Always skip the first frame...
        lims(2) = length(flist); % Go to the end.
    end;
        
    % Always avoid the first frame of the series (catch user errors)
    if lims(1) == 1; lims(1:end-1) = lims(2:end); end;

    % Pre-allocate the main variables
    im(lims(2)).orig = zeros(600,800); % This is the biggest memory hog - the stack of images.

    
%% This is our main loop to process every frame.

for ff = lims(1):lims(2);

    im(ff).orig = imread(flist(ff).name); 
    im(ff).fname = flist(ff).name;
     
     bw = im2bw(im(ff).orig, thresh); % Make a BW image with threshold 'thresh'
 
       for y=3:-1:1; mybm(:,:,y) = double(bw); end;
       mybm(:,:,3) = double(ones(600,800));                   
       mybm(:,:,2) = double(zeros(600,800));                  


%% We need to treat the first frame differently from the rest - get our clicks
if ff == lims(1);

    figure(3);

    imshow(mybm);

% Step one, get medial limits

    [rlim, rhinge] = ginput(1);
    hold on; 
    plot([rlim, rlim], [1 599], 'w-', 'LineWidth', 2);
    plot([1 rlim], [rhinge rhinge], 'y-', 'LineWidth', 2);

    [llim, lhinge] = ginput(1);
    hold on; 
    plot([llim, llim], [1 599], 'w-', 'LineWidth', 2);
    plot([llim, 799], [lhinge, lhinge], 'y-', 'LineWidth', 2);
    
    vlim = 565; % Ventral limit is to prevent problems with numeric labels

% Step two, get the leg/antenna limits

    [rAnt, ~] = ginput(1);
    plot([rAnt rAnt], [1 599], 'y-', 'LineWidth', 4);
    [lAnt, ~] = ginput(1);
    plot([lAnt lAnt], [1 599], 'y-', 'LineWidth', 4);

% Step three, get the current wing positions    

    [Xrc, Yrc] = ginput(1);
        plot(Xrc, Yrc, 'r*');    
    [Xlc, Ylc] = ginput(1);
        plot(Xlc,Ylc, 'g*');
    
% Step 4 - plot everything

    plot([Xrc-boxlim(1) Xrc-boxlim(1)], [Yrc-boxlim(2) Yrc+boxlim(2)], 'w-', 'LineWidth', 2);
    plot([Xlc+boxlim(1) Xlc+boxlim(1)], [Ylc-boxlim(2) Ylc+boxlim(2)], 'w-', 'LineWidth', 2);
    
    plot([Xrc-boxlim(1) rlim], [Yrc-boxlim(2) Yrc-boxlim(2)], 'w-', 'LineWidth', 2);
    plot([Xrc-boxlim(1) rlim], [Yrc+boxlim(2) Yrc+boxlim(2)], 'w-', 'LineWidth', 2);
    
    plot([llim Xlc+boxlim(1)], [Ylc-boxlim(2) Ylc-boxlim(2)], 'w-', 'LineWidth', 2);
    plot([llim Xlc+boxlim(1)], [Ylc+boxlim(2) Ylc+boxlim(2)], 'w-', 'LineWidth', 2);

% Step 5 - make our first 'search boxes' - rbox on right and lbox on left
    rbox = [Xrc-boxlim(1), Yrc-boxlim(2);...
        rlim, Yrc-boxlim(2);...
        rlim, Yrc+boxlim(2);...
        Xrc-boxlim(1), Yrc+boxlim(2)];
    lbox = [llim, Ylc-boxlim(2);...
        Xlc+boxlim(1), Ylc-boxlim(2);...
        Xlc+boxlim(1), Ylc+boxlim(2);...
        llim, Ylc+boxlim(2)];

    if rbox(1,2) <= 0; rbox(1,2) = 1; rbox(2,2) = 1; end; % Fixes right box that goes below zero
    if lbox(1,2) <= 0; lbox(1,2) = 1; lbox(2,2) = 1; end; % Fixes left box that goes below zero
    if lbox(2,1) > 800; lbox(2,1) = 800; lbox(3,1) = 800; end; % Fixes left box that goes beyond edge of image
    if rbox(3,2) > vlim; rbox(3,2) = vlim; rbox(4,2) = vlim; end; % Fixes right box that goes below bottom
    if lbox(3,2) > vlim; lbox(3,2) = vlim; lbox(4,2) = vlim; end; % Fixes left box that goes below bottom
    rbox = int16(rbox); % convert to a 16 bit integer.
    lbox = int16(lbox); % convert to a 16 bit integer.

    Xrc(2:ff) = Xrc(1);
    Xlc(2:ff) = Xlc(1);
    Yrc(2:ff) = Yrc(1);
    Ylc(2:ff) = Ylc(1);

end;

pause(0.01); hold off;

%% For all other frames, we automagically find the wing locations

% Right side X
clear X; X = zeros(1,rbox(2,1)-rbox(1,1)); % PreAllocate for speed
     for m = rbox(1,1):rbox(2,1);
         X(m) = sum(bw(rbox(1,2):rbox(3,2),m));
     end;
     rxdist = find(X > cnt);
     Xrc(ff) = rxdist(1);

% Left side X   
clear X; X = zeros(1,lbox(2,1)-lbox(1,1)); % PreAllocate for speed
     for m = lbox(1,1):lbox(2,1);
         X(m) = sum(bw(lbox(1,2):lbox(3,2),m));
     end;
     lxdist = find(X > cnt); 
     Xlc(ff) = lxdist(end);

% Right side Y
clear Y; Y = zeros(1,rbox(3,2)-rbox(1,2)); % PreAllocate for speed   
    if rbox(3,2) < rhinge; % Top quadrant
     for m = rbox(1,2):rbox(3,2);
         Y(m) = sum(bw(m,rbox(1,1):rbox(2,1)));
     end;
        rydist = find(Y > cnt);
            if isempty(rydist);
                ff
                figure(1); subplot(121); imshow(im(ff-1).orig); subplot(122); imshow(im(ff).orig);
             return;
            end;
        Yrc(ff) = rydist(1); 
    end;
    if rbox(1,2) > rhinge; % Bottom quadrant
     for m = rbox(1,2):rbox(3,2);
         Y(m) = sum(bw(m,rbox(1,1):rAnt));
     end;
     rydist = find(Y > cnt);
        Yrc(ff) = rydist(end);
    end;
    if length(Yrc) < ff; % Middle quadrant
     for m = rbox(1,2):rbox(3,2);
         Y(m) = sum(bw(m,rbox(1,1):rAnt));
     end;
     rydist = find(Y > cnt);

        if Yrc(ff-2) < Yrc(ff-1); % Down stroke
            Yrc(ff) = rydist(end);
        else
            Yrc(ff) = rydist(1);
        end;
    end;

% Left side Y
clear Y; Y = zeros(1,lbox(3,2)-lbox(1,2)); % PreAllocate for speed

    if lbox(3,2) < lhinge; % Top quadrant
     for m = lbox(1,2):lbox(3,2);
         Y(m) = sum(bw(m,lbox(1,1):lbox(2,1)));
     end;
     lydist = find(Y > cnt);
         if isempty(lydist);
                ff
             figure(1); subplot(121); imshow(im(ff-1).orig); subplot(122); imshow(im(ff).orig);
             return;
         end;

        Ylc(ff) = lydist(1); 
    end;
    if lbox(1,2) > lhinge; % Bottom quadrant
     for m = lbox(1,2):lbox(3,2);
         Y(m) = sum(bw(m,lAnt:lbox(2,1)));
     end;
     lydist = find(Y > cnt);
         if isempty(lydist);
                ff
             figure(1); subplot(121); imshow(im(ff-1).orig); subplot(122); imshow(im(ff).orig);
             return;
         end;
        Ylc(ff) = lydist(end);
    end;
    if length(Ylc) < ff; % Middle quadrant
     for m = lbox(1,2):lbox(3,2);
         Y(m) = sum(bw(m,lAnt:lbox(2,1)));
     end;
     lydist = find(Y > cnt);
          if isempty(lydist);
                ff
             figure(1); subplot(121); imshow(im(ff-1).orig); subplot(122); imshow(im(ff).orig);
             return;
         end;

    if Ylc(ff-2) < Ylc(ff-1); % Down stroke
            Ylc(ff) = lydist(end);
        else
            Ylc(ff) = lydist(1);
     end;
     
    end;

    % Update the box locations
        rbox = [Xrc(ff)-boxlim(1), Yrc(ff)-boxlim(2);...
        rlim, Yrc(ff)-boxlim(2);...
        rlim, Yrc(ff)+boxlim(2);...
        Xrc(ff)-boxlim(1), Yrc(ff)+boxlim(2)];
    
        lbox = [llim, Ylc(ff)-boxlim(2);...
        Xlc(ff)+boxlim(1), Ylc(ff)-boxlim(2);...
        Xlc(ff)+boxlim(1), Ylc(ff)+boxlim(2);...
        llim, Ylc(ff)+boxlim(2)];
    % Same fixes as above... don't want the boxes running outside of the
    % image
    if rbox(1,2) <= 0; rbox(1,2) = 1; rbox(2,2) = 1; end;
    if lbox(1,2) <= 0; lbox(1,2) = 1; lbox(2,2) = 1; end;
    if lbox(2,1) > 800; lbox(2,1) = 800; lbox(3,1) = 800; end;
    if rbox(3,2) > vlim; rbox(3,2) = vlim; rbox(4,2) = vlim; end;
    if lbox(3,2) > vlim; lbox(3,2) = vlim; lbox(4,2) = vlim; end;
    rbox = int16(rbox);
    lbox = int16(lbox);
    
    
end;    

%% Put the data into the output structure

data.l.x = Xlc;
data.l.y = Ylc;
data.r.x = Xrc;
data.r.y = Yrc;

data.l.imnum = lims(1):1:lims(2);
data.r.imnum = lims(1):1:lims(2);

% Determine whether the wing is going up or down on both sides

data.l.up = find(Ylc(2:end) - Ylc(1:end-1) >= 0) + 1;
data.r.up = find(Yrc(2:end) - Yrc(1:end-1) >= 0) + 1;

% Calculate wing velocities
for p = lims(1)+1:lims(2);    
    data.r.vel(p) = -1 * pdist([data.r.x(p-1),data.r.y(p-1); data.r.x(p), data.r.y(p)]);
    data.l.vel(p) = -1 * pdist([data.l.x(p-1),data.l.y(p-1); data.l.x(p), data.l.y(p)]);
end;

data.r.vel(data.r.up) = -data.r.vel(data.r.up);
data.l.vel(data.l.up) = -data.l.vel(data.l.up);

hold on;
