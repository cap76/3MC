%=========================================================================
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%==========================================================================

function createTimer(pHndl, action)

%displays animated gif in the axes specified by pHndl

persistent hImage tm

try
    
    if strcmp(action, 'start')
        img = imread('processing.gif', 'frames', 'all');
        %img is 32-by-32-by-1-by-8
        %there are 8 32*32 frames, The 1 dimension meanit is a 1 bit
        %(grayscale) image, ie would be 3 if RGB image
        %map is 16-by-3, representing 16 grays, each with all 3 rgb values
        %the same
        colormap(pHndl, 'gray');
        
        
        img = img+40;
        %idx goes from 0 (dark grey) to 11 (white)
        %want 11 to be figure colour then convert to colourmap dark grey
        %values [40,51]
        
        %show first frame
        hImage = image(img(:,:,:,1), 'parent', pHndl);
        set(pHndl, 'box', 'off', 'xtick', [], 'ytick', [], 'visible', 'off');
        
        %start timer
        per = 1.5/size(img, 4);
        per = sprintf('%3.2f', per);
        per = str2num(per);
        tm = timer('TimerFcn', {@timer_update, img, hImage}, 'ExecutionMode', 'fixedRate', 'Period', per);
        start(tm);
        
    else
        %finished with image
        stop(tm);
        delete(tm);
        delete(hImage);
    end
catch
    %some error, just do without the image
    return;
end


%===========================================
function timer_update(tm, data, imgdata, hImage)

persistent idx

if isempty(idx) || idx >= size(imgdata, 4)
    idx = 1;
else
    idx = idx + 1;
end

set(hImage, 'CData', imgdata(:,:,:,idx));








%%%%%%%%The old version that contructs the image from circle objects

function oldcreateTimer(fig, action, position)


% draws the animation that rotates as program runs
%
% fig - the figure to draw on
% position - [x y r] a ring of circles of radius r centered around x,y
%

global timer_idx circles tm;

if strcmp(action, 'start')
    
    centre = position(1:2);
    radius = position(3);
    
    u = get(fig, 'Units');
    set(fig, 'Units', 'pixels');
    pos = get(fig, 'position');
    set(fig, 'units', u);
    
    aspect_ratio = pos(3)/pos(4);
    
    %centres of individual circles around the circumference of the main circle
    angle = 0;
    coord = zeros(8,2);
    for c = 1:8
        coord(c,1) = cos(angle) * radius / aspect_ratio;%x
        coord(c,2) = sin(angle) * radius;%y
        angle = angle + pi/4;
    end
    
    for c = 1:8
        circles(c) = annotation(fig, 'ellipse', [(centre(1)+ coord(c,1)) (centre(2)+coord(c,2)) 0.5/100 0.7/100], 'EdgeColor', [0.7 0.7 0.7], 'FaceColor', [0.7 0.7 0.7]);
    end
    
    timer_idx = 1;
    set(circles(timer_idx),'EdgeColor', [0.5 0.5 0.5], 'FaceColor', [0.5 0.5 0.5]); 
    set(circles(timer_idx+1),'EdgeColor', [0.6 0.6 0.6], 'FaceColor', [0.6 0.6 0.6]); 
    
    tm = timer('TimerFcn', @timer_update, 'ExecutionMode', 'fixedRate', 'Period', 0.15);
    
    start(tm);
    
   
else
    stop(tm);
    delete(tm);
    delete(circles);    
end

%===========================================
function oldtimer_update(tm, data)

global timer_idx circles;

set(circles(timer_idx),'EdgeColor', [0.6 0.6 0.6], 'FaceColor', [0.6 0.6 0.6]);
if timer_idx == 8
    set(circles(1),'EdgeColor', [0.7 0.7 0.7], 'FaceColor', [0.7 0.7 0.7]);
else
    set(circles(timer_idx+1),'EdgeColor', [0.7 0.7 0.7], 'FaceColor', [0.7 0.7 0.7]);
end

if timer_idx == 1
    set(circles(timer_idx+1),'EdgeColor', [0.7 0.7 0.7], 'FaceColor', [0.7 0.7 0.7]);
    timer_idx = 8;
else
    timer_idx = timer_idx - 1;
end

set(circles(timer_idx),'EdgeColor', [0.5 0.5 0.5], 'FaceColor', [0.5 0.5 0.5]);
 

