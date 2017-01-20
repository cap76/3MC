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



function chromosomegui(varargin)



if nargin == 0
    %start program
    action = 'init';
 
else
    action = varargin{3};
end

UserEvent(action);


%=============================================
function UserEvent(varargin)

%called whenever the user does something
%eg starting the program, or clicking  a control

persistent  closeHndl runHndl mainFig browseFileHndl


%inputs
persistent mcHndl nochnHndl cenHndl lenHndl radiusHndl uHndl kappaHndl pHndl fileHndl chrHndl perHndl dirHndl archHndl

args = varargin{1};
if iscell(args)
    action = args{1};
else
    action = args;
end

if strcmp(action,'init')
    
    if ~isdeployed
        addpath(genpath(pwd));
    end
    
    architectures = {'Tethered Tel', ...
                        'Untethered'};
                    
    architectures2 = {'Gaussian Process', ...
                      'Load from kappa.mat'};                    
    
    %create the gui and controls
    mainFig=figure('resize', 'off', 'Menubar', 'none' ,'NumberTitle','off','Visible','off', 'CloseRequestFcn', {@chromosomegui, 'close'}, 'name', '3MC');
    
    set(0,'Units','centimeters')
    screen_size = get(0,'ScreenSize');
    
    %centre fig on screen
    figwidth = 14;
    figheight = 19;
    set(mainFig, 'Units', 'centimeters', 'Position', [(screen_size(3)-figwidth)/2 (screen_size(4)-figheight)/2 figwidth figheight]);
    
    maincol = get(mainFig, 'Color');
    panelwidth = figwidth-0.5;
    panelheight = figheight-1.4;
    frmPos=[0.25 1.15 panelwidth  panelheight];
  
    
    %draw panel with controls in the required position
     panel = uipanel('BorderType', 'etchedin', ...
        'BackgroundColor', get(mainFig, 'Color'), ...
        'Units','centimeters', ...
        'Position',frmPos, ...
        'Parent', mainFig);
    
    ctrlsep = 1;
    ctrlbase = panelheight-ctrlsep;
    ctrlleft = 0.75;
    fsz = 10;

    uicontrol('FontWeight', 'bold','HorizontalAlignment', 'left', 'Parent',panel ,'Style', 'text','Units','centimeters','position',[ctrlleft-0.25 ctrlbase panelwidth-1 0.5],'string','Select the input parameters for Chromosome Analysis', 'ForegroundColor', 'k', 'BackgroundColor', get(panel, 'backgroundcolor'), 'FontUnits', 'points', 'FontSize', fsz);    
    
    %Parameters
    ctrlbase=ctrlbase-ctrlsep;
    chrHndl = createControl('edit', ctrlleft, ctrlbase, 'Chromosome length (base pairs):', '230000', panel);
    
    ctrlbase=ctrlbase-ctrlsep*1.25;    
    perHndl = createControl('popup', ctrlleft, ctrlbase, 'Specify persistence-length function:', architectures2, panel);
          
    %ctrlbase=ctrlbase-ctrlsep*1.25;
    %archHndl = createControl('popup', ctrlleft, ctrlbase, 'Select the tethering condition:', architectures, panel);

    
    
    %[file path]=uigetfile('*.jpg','Select any image');
    ctrlbase=ctrlbase-ctrlsep;   
    pHndl = createControl('edit', ctrlleft, ctrlbase, 'Compaction factor (phi):', '40', panel);
    ctrlbase=ctrlbase-ctrlsep;
    lenHndl = createControl('edit', ctrlleft, ctrlbase, {'Number of beads in the polymer chain:', '(segements in the chromosome)'}, '100', panel); 
    ctrlbase=ctrlbase-ctrlsep*1.25;
    cenHndl = createControl('edit', ctrlleft, ctrlbase, 'Position of centromere (bead number):', '66', panel);
  
    ctrlbase=ctrlbase-ctrlsep*1.25;
    kappaHndl = createControl('edit', ctrlleft, ctrlbase, 'Clustering parameter (kappa):', '10', panel);
    ctrlbase=ctrlbase-ctrlsep;
    dirHndl = createControl('grid', ctrlleft, ctrlbase, {'Direction  of clustering:', 'x', 'y', 'z'}, [0 0 1], panel);
    ctrlbase=ctrlbase-ctrlsep*2.1;
    mcHndl = createControl('edit', ctrlleft, ctrlbase, 'Number of steps in the Markov Chain:', '200000', panel);
    ctrlbase=ctrlbase-ctrlsep;
    nochnHndl = createControl('edit', ctrlleft, ctrlbase, 'Number of Chains:', '2', panel);
    
    ctrlbase=ctrlbase-ctrlsep*1.25;
    radiusHndl = createControl('edit', ctrlleft, ctrlbase, 'Radius of the nuclear periphery (nm):', '1600', panel);
    ctrlbase=ctrlbase-ctrlsep;
    uHndl = createControl('edit', ctrlleft, ctrlbase, 'Confinement potential for beads outside of the NP:', '10000', panel);
    
    ctrlbase=ctrlbase-ctrlsep*1.25;
    archHndl = createControl('popup', ctrlleft, ctrlbase, 'Select the tethering condition:', architectures, panel);
    
    ctrlbase=ctrlbase-ctrlsep;
    uicontrol('FontWeight', 'bold','HorizontalAlignment', 'left', 'Parent',panel ,'Style', 'text','Units','centimeters','position',[ctrlleft-0.25 ctrlbase panelwidth-1 0.5],'string','File name to store results', 'ForegroundColor', 'k', 'BackgroundColor', get(panel, 'backgroundcolor'), 'FontUnits', 'points', 'FontSize', fsz);    
    ctrlbase=ctrlbase-ctrlsep+0.25;
    hndls = createControl('file', ctrlleft, ctrlbase, '', '', panel);

    fileHndl = hndls(1);browseFileHndl=hndls(2);
    
     
    % control buttons
    closeHndl = uicontrol(...
        'Style','pushbutton', ...
        'Units','centimeters', ...
        'Position',[figwidth-2.75 0.25 2.5 0.8], ...
        'Parent',mainFig, ...
        'string', 'Close', ...
        'FontUnits', 'points', 'FontSize', 10, ...
        'Callback',{@chromosomegui, 'close'});
    runHndl = uicontrol(...
        'Style','pushbutton', ...
        'Units','centimeters', ...
        'Position',[figwidth-5.5 0.25 2.5 0.8], ...
        'Parent',mainFig, ...
        'string', 'Run', ...
        'FontUnits', 'points', 'FontSize', 10, ...
        'Callback',{@chromosomegui, 'run'});
        
    set(mainFig,'Visible','on');

elseif strcmp(action, 'run')

    %launches  analysis
    [ok data] = ValidControls(mcHndl, nochnHndl, cenHndl, lenHndl, radiusHndl, uHndl, kappaHndl, pHndl, chrHndl, perHndl, dirHndl, fileHndl, archHndl);

    if ok
       %launch progress window
       RunChrom('go', data, mainFig);
    end

elseif strcmp(action, 'close')
    
    delete(gcf);
    
elseif strcmp(action, 'selectfile')
    
    %browsing for a file
    btn = gcbo;
    fileNameBox = get(btn, 'Userdata');
    
    [fname, pathname] = uiputfile('*.mat', 'Save Results As');
    %what about new file and extension???
    if ~isequal(fname, 0)
        %valid file selected
        fname = fullfile(pathname, fname);
        set(fileNameBox, 'string', fname);
    end
end

%=========================================================================

function hndl = createControl(style, ctrlleft, ctrlbase, caption, init, parent)

ctrlheight=0.6;
captHeight = 0.5;
captionwidth = 9.5;
editwidth=2.25;
fsz = 10; %font size

if strcmp(style, 'edit')
    %caption
    if iscell(caption)
        for i = 1:length(caption)
            cHndl{i} = uicontrol('HorizontalAlignment', 'left', 'Parent',parent ,'Style', 'text','Units','centimeters','position',[ctrlleft ctrlbase-(0.5*(i-1)) captionwidth captHeight],'string', caption{i}, 'ForegroundColor', 'k', 'BackgroundColor', get(parent, 'backgroundcolor'), 'FontUnits', 'points', 'FontSize', fsz);    
        end
    else
        cHndl = uicontrol('HorizontalAlignment', 'left', 'Parent',parent ,'Style', 'text','Units','centimeters','position',[ctrlleft ctrlbase captionwidth captHeight],'string', caption, 'ForegroundColor', 'k', 'BackgroundColor', get(parent, 'backgroundcolor'), 'FontUnits', 'points', 'FontSize', fsz);    
    end
    %control contains a reference to its caption , in case this needs to be
    %set dynamically
    hndl = uicontrol('userdata', cHndl, 'HorizontalAlignment', 'right', 'Parent',parent ,'Style', 'edit','Units','centimeters','position',[captionwidth+ctrlleft ctrlbase-0.1 editwidth ctrlheight],'string',init, 'ForegroundColor', 'k', 'BackgroundColor', 'w','FontUnits', 'points', 'FontSize', fsz);    

elseif strcmp(style, 'file')

     fileHndl =uicontrol( ...
        'Style','edit', ...
        'HorizontalAlignment', 'left', ...
        'Units','centimeters', ...
        'position',[ctrlleft ctrlbase captionwidth ctrlheight], ...
        'Parent',parent, ...
        'FontUnits', 'points', 'FontSize', fsz, ...
        'String', init, ...
        'BackgroundColor', 'w');
    browseHndl =uicontrol( ...
        'Style','pushbutton', ...
        'Units','centimeters', ...
        'position', [captionwidth+ctrlleft ctrlbase editwidth ctrlheight], ...
        'Parent',parent, ...
        'FontUnits', 'points', 'FontSize', fsz, ...
        'String', 'Browse...', ...
        'call',{@chromosomegui, 'selectfile'}, ...
        'Userdata', fileHndl);
    hndl = [fileHndl browseHndl];
    
elseif strcmp(style, 'popup')
    
    %caption
    cHndl = uicontrol('fontweight', 'bold', 'HorizontalAlignment', 'left', 'Parent',parent ,'Style', 'text','Units','centimeters','position',[ctrlleft-0.25 ctrlbase captionwidth captHeight],'string', caption, 'ForegroundColor', 'k', 'BackgroundColor', get(parent, 'backgroundcolor'), 'FontUnits', 'points', 'FontSize', fsz);    
    %control contains a reference to its caption , in case this needs to be
    %set dynamically
    hndl = uicontrol('userdata', cHndl, 'HorizontalAlignment', 'left', 'Parent',parent ,'Style', 'popup','Units','centimeters','position',[captionwidth+ctrlleft-(editwidth/2) ctrlbase-0.1 editwidth*1.5 ctrlheight],'string',init, 'ForegroundColor', 'k', 'BackgroundColor', 'w','FontUnits', 'points', 'FontSize', fsz);    
    
elseif strcmp(style, 'grid')
   
    cHndl = uicontrol('HorizontalAlignment', 'left', 'Parent',parent ,'Style', 'text','Units','centimeters','position',[ctrlleft ctrlbase captionwidth captHeight],'string', caption, 'ForegroundColor', 'k', 'BackgroundColor', get(parent, 'backgroundcolor'), 'FontUnits', 'points', 'FontSize', fsz);    

    numrows = length(init);
    pos = [captionwidth+ctrlleft ctrlbase-((numrows-1)*ctrlheight) editwidth ctrlheight*numrows];
    
    hndl = uitable('userdata', cHndl, 'columneditable', [true], 'columnname', [], 'columnformat', {'numeric'}, 'rowstriping', 'off', 'Parent',parent ,'Units','centimeters','position', pos, 'FontUnits', 'points', 'FontSize', fsz);    
    set(hndl, 'units', 'pixels');
    set(hndl, 'columnwidth', {[40]}); 
    set(hndl, 'units', 'centimeters');
    %data, rowname
    data = {}; rowname = {};
    for i = 1:numrows
       data = [data; init(i)];
       rowname{i} = caption{i+1};
    end
    set(hndl, 'data', data, 'rowname', rowname);
end


%==========================================================================

function [ok data] = ValidControls(numsteps, numchain, cen, numbeads, radius, confpot, kappa, phi, chrlen, perlen, dirtbl, filename, archhndl)


ok = 0;
data = [];

%validate parameters entered

No = str2double(get(numsteps, 'String'));
if isempty(No) || No < 1 || isnan(No)
    ShowError('Please enter a valid number of steps in the Markov chain', '3MC');
    uicontrol(numsteps);
    return;
else
    data.numsteps = No;
end


% ADDED by CP
numchain = str2double(get(numchain, 'String'));
if isempty(numchain) || numchain < 1 || isnan(numchain)
    ShowError('Please enter a valid number of Markov chains', '3MC');
    uicontrol(numchain);
    return;
else
    data.numchain = numchain;
end


beads = str2double(get(numbeads, 'String'));
if isempty(beads) || beads < 10 || isnan(beads)
    ShowError('Please enter a valid number of beads, 10 or greater', '3MC');
    uicontrol(numbeads);
    return;
else
    data.numbeads = beads;
end

p= str2double(get(cen, 'String'));
if isempty(p) || p > beads || isnan(p) || p < 1
    ShowError('Please enter a valid centromere position', '3MC');
    uicontrol(cen);
    return;
else
    data.cen = p;
end

p= str2double(get(phi, 'String'));
if isempty(p) || p >= 500 || isnan(p) || p <= 0
    ShowError('Please enter a valid compaction factor value (0 < phi < 500)', '3MC');
    uicontrol(phi);
    return;
else
    data.phi = p;
end

rad = str2double(get(radius, 'String'));
if isempty(rad) || rad < 1 || isnan(rad)
    ShowError('Please enter a valid radius value', '3MC');
    uicontrol(radius);
    return;
else
    data.radius = rad;
end

u = str2double(get(confpot, 'String'));
if isempty(u) || u < 1 || isnan(u)
    ShowError('Please enter a valid confinement potential value', '3MC');
    uicontrol(confpot);
    return;
else
    data.confpotential = u;
end

kp = str2double(get(kappa, 'String'));
if isempty(kp) || isnan(kp) || kp <= 0 || kp > 20 
    ShowError('Please enter a valid telomere clustering parameter (0 < kappa <= 20)', '3MC');
    uicontrol(kappa);
    return;
else
    data.kappa = kp;
end

chr = str2double(get(chrlen, 'string'));
if isempty(chr) || isnan(chr) || chr <= 0 
    ShowError('Please enter a valid chromosome length, in base pairs', '3MC');
    uicontrol(chrlen);
    return;
else
   data.chrlen = chr; 
end

%per = str2double(get(perlen, 'string'));
%if isempty(per) || isnan(per) || per <= 0 
%    ShowError('Please enter a valid persistence length value', '3MC');
%    uicontrol(perlen);
%    return;
%else
%   data.perlen = per; 
%end

arch2 = get(perlen, 'string')
data.architecture2 = arch2{get(perlen, 'value')}

dircluster = cell2mat(get(dirtbl, 'data'));
dircluster = dircluster./norm(dircluster);
em = 1e-3; %allowed error margin
if sum(dircluster.^2) > (1+em) || sum(dircluster.^2) < (1-em)
    ShowError(['Please enter a valid clustering direction vector. The sum of squares of the x, y and z elements must equal one. The margin of error allowed is ' num2str(em) '.'], '3MC');
    return;
else
   data.dir = dircluster; 
end

arch = get(archhndl, 'string')
data.architecture = arch{get(archhndl, 'value')}

%arch = get(archhndl, 'string');
%data.architecture = arch{get(archhndl, 'value')};

fname = get(filename, 'string');
if isempty(fname) || fname(end) == filesep %entering a path not a file
    ShowError('Please enter a file name', '3MC');
    uicontrol(filename);
    return;
else
    %check the directory exists
    [pathstr, name, ext] = fileparts(fname);
    if exist(pathstr, 'dir') ~= 7
        ShowError('The selected directory does not exist', '3MC');
        uicontrol(filename);
        return;
    else
        %ensure correct extension
        if ~strcmp(ext, '.mat')
            fname = [fname '.mat'];
        end
        
        %check new name is legal
        badchars = '\/:*?"<>|`';
        if sum(ismember(name, badchars))
            ShowError(['The file name cannot contain any of the following characters: ' badchars], '3MC');
            uicontrol(filename);
            return;
        end
        
    end
end

data.filename = fname;


 
a = data.chrlen*0.34/(data.phi*(data.numbeads-1));
%k = data.perlen*data.chrlen*0.34/data.phi;
%data.k = k;
data.a = a;

ok = 1;



