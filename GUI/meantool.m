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


function meantool()

%Generates a gui asking user if they wish to run a new analysis or view an
%existing one.

if ~(isdeployed || ismcc)
    mydir = fileparts(mfilename('fullpath'));
    if ~strcmp(mydir, pwd)
        error('The tool must be launched from its installation directory'); 
    end
% else
%     if ~strcmp(ctfroot, pwd)
%         error('The tool must be launched from its installation directory'); 
%     end
%     %need to include paths somehow
%     
end

ctrlleft = 0.5;fsz=10;

newfig=figure('resize', 'off', 'Menubar', 'none' ,'NumberTitle','off','Visible','off', 'name', '3MC');

set(0,'Units','centimeters')
screen_size = get(0,'ScreenSize');

%centre fig on screen
figwidth = 14;
figheight = 8;
set(newfig, 'Units', 'centimeters', 'Position', [(screen_size(3)-figwidth)/2 (screen_size(4)-figheight)/2 figwidth figheight]);

maincol = get(newfig, 'Color');
panelwidth = figwidth-1;
panelheight = figheight-2;
frmPos=[0.5 1.5 panelwidth  panelheight];


%draw panel with controls in the required position
panel = uipanel('BorderType', 'etchedin', ...
    'BackgroundColor', maincol, ...
    'Units','centimeters', ...
    'Position',frmPos, ...
    'Parent', newfig);


uicontrol('FontWeight', 'bold','HorizontalAlignment', 'center', 'Parent',panel ,'Style', 'text','Units','centimeters','position',[ctrlleft panelheight-1.5 panelwidth-1 0.75],'string','Welcome to 3MC', 'ForegroundColor', 'k', 'BackgroundColor', get(panel, 'backgroundcolor'), 'FontUnits', 'points', 'FontSize', fsz+4);
uicontrol('FontWeight', 'bold','HorizontalAlignment', 'center', 'Parent',panel ,'Style', 'text','Units','centimeters','position',[ctrlleft panelheight-2.25 panelwidth-1 0.5],'string','What would you like to do?', 'ForegroundColor', 'k', 'BackgroundColor', get(panel, 'backgroundcolor'), 'FontUnits', 'points', 'FontSize', fsz+2);

uicontrol( ...
    'Style','pushbutton', ...
    'Units','centimeters', ...
    'position',[(panelwidth-6)/2  panelheight-3.5 6 0.7], ...
    'Parent',panel, ...
    'FontUnits', 'points', 'FontSize', fsz, ...
    'String', 'Generate new trajectories', ...
    'call','chromosomegui');

uicontrol( ...
    'Style','pushbutton', ...
    'Units','centimeters', ...
    'position',[(panelwidth-6)/2  panelheight-4.5 6 0.7], ...
    'Parent',panel, ...
    'FontUnits', 'points', 'FontSize', fsz, ...
    'String', 'Analyse existing trajectorie', ...
    'call','plotchromosome');

% control buttons
uicontrol(...
    'Style','pushbutton', ...
    'Units','centimeters', ...
    'Position',[figwidth-3 0.5 2.5 0.8], ...
    'Parent',newfig, ...
    'string', 'Close', ...
    'FontUnits', 'points', 'FontSize', 10, ...
    'Callback','delete(gcf);');


set(newfig,'Visible','on');
