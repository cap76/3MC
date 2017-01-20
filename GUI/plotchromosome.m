function plotchromosome(varargin)

persistent  closeHndl runHndl plotfig browseHndl fileHndl tblHndl filedata  plotboxes widthHndl

global plot_titles plot_funcs

plot_titles= {...
   % 'Radial Distribution Function (Beta fit)'; ...
    %'Locus-Locus Distance (Beta fit)'; ...
    'Intrachromosome Distances'; ...
    'Distance between Homologues'; ...
    'Example Trajectory'...
    %'Sample Chromosome Trajectory'; ...
    };

plot_funcs = {  @plot_vbm;...
                @plot_cdm;...
                @plot_pd;...
                %@plot_polymerpic;...
                };
                    


if nargin == 0
    action = 'init';
elseif ischar(varargin{1})
    %called in code
    action = varargin{1};
else
    %called from a callback like this plotchromosome(src, event, action, ....);
   action = varargin{3}; 
end

if strcmp(action,'init')
    %create the gui and controls
    if ~isdeployed
        addpath(genpath(pwd));
    end
    
    fsz = 10;
    ctrlleft = 0.5;
    
    plotfig=figure('resize', 'off', 'Menubar', 'none' ,'NumberTitle','off','Visible','off', 'name', '3MC');
    
    set(0,'Units','centimeters')
    screen_size = get(0,'ScreenSize');
    
    %centre fig on screen
    figwidth = 14;
    figheight = 17;
    set(plotfig, 'Units', 'centimeters', 'Position', [(screen_size(3)-figwidth)/2 (screen_size(4)-figheight)/2 figwidth figheight]);
    
    maincol = get(plotfig, 'Color');
    panelwidth = figwidth-1;
    panelheight = figheight-2;
    frmPos=[0.5 1.5 panelwidth  panelheight];
  
    
    %draw panel with controls in the required position
     panel = uipanel('BorderType', 'etchedin', ...
        'BackgroundColor', get(plotfig, 'Color'), ...
        'Units','centimeters', ...
        'Position',frmPos, ...
        'Parent', plotfig);
    
   
    uicontrol('FontWeight', 'bold','HorizontalAlignment', 'left', 'Parent',panel ,'Style', 'text','Units','centimeters','position',[ctrlleft panelheight-1.25 panelwidth-1 0.5],'string','Select the results file', 'ForegroundColor', 'k', 'BackgroundColor', get(panel, 'backgroundcolor'), 'FontUnits', 'points', 'FontSize', fsz);    
     %file browser 
    fileHndl =uicontrol( ...
        'Style','edit', ...
        'HorizontalAlignment', 'left', ...
        'Units','centimeters', ...
        'position',[ctrlleft panelheight-2 panelwidth-3 0.6], ...
        'Parent',panel, ...
        'FontUnits', 'points', 'FontSize', 10, ...
        'String', [], ...
        'BackgroundColor', 'w');
    browseHndl =uicontrol( ...
        'Style','pushbutton', ...
        'Units','centimeters', ...
        'position',[panelwidth-2.5 panelheight-2 2 0.6], ...
        'Parent',panel, ...
        'FontUnits', 'points', 'FontSize', 10, ...
        'String', 'Browse ...', ...
        'call','plotchromosome(''selectfile'');');
        
    
    
    %details of selected file
    data = { ...
            'Chromosome length (base pairs)', ' '; ...           %chrlen
            'Persistence length (nm)', ' '; ...                  %perlen
            'Compaction factor (phi)', ' '; ...                  %phi
            'Beads in ploymer chain', ' '; ...                   %numbeads
            'Centromere position', ' '; ...                      %cen            
            'Rigidity constant (k)', ' '; ...                    %k
            'Segment length (a)', ' '; ...                       %a
            'Telomere clustering parameter (kappa)', ' '; ...    %kappa
            'Clustering direction [x y z]', ' '; ...             %dir
            'Steps in Markov chain', ' '; ...                    %numsteps
            'Number of Markov chains', ' '; ...                  %numchain            
            'Radius of nuclear periphery (nm)', ' '; ...         %radius
            'Confinement potential', ' '; ...                    %confpotential
            'Algorithm', ' '; ...                                %architecture
            };

    tblHndl = uitable( ...
                        'Units','centimeters', ...
                        'FontUnits', 'points', 'FontSize', 10, ...
                        'parent', panel, ....
                        'position', [ctrlleft panelheight-10.75 panelwidth-3 8], ...
                        'columnname', [], 'rowname', [], 'data', data, 'userdata', data);
   set(tblHndl, 'units', 'pixels');
   set(tblHndl, 'columnwidth', {[270] [100]});
    set(tblHndl, 'units', 'centimeters');
    
    
      
   %plot options
  uicontrol('FontWeight', 'bold','HorizontalAlignment', 'left', 'Parent',panel ,'Style', 'text','Units','centimeters','position',[ctrlleft panelheight-11.75 panelwidth-1 0.5],'string','Select the required analysis', 'ForegroundColor', 'k', 'BackgroundColor', get(panel, 'backgroundcolor'), 'FontUnits', 'points', 'FontSize', fsz);    
  plotboxes=zeros(1,length(plot_titles));
  width=(panelwidth-1)/2;
  for i = 1:length(plot_titles)
      pos = [ctrlleft+(width*mod(i-1, 2)) (panelheight-12.5)-(floor((i-1)/2))/2 width 0.5];
      plotboxes(i) = uicontrol('max', 1, 'min', 0, 'HorizontalAlignment', 'right', 'Parent',panel ,'Style', 'checkbox','Units','centimeters','position',pos,'string',plot_titles{i}, 'ForegroundColor', 'k', 'BackgroundColor', get(panel, 'backgroundcolor'), 'FontUnits', 'points', 'FontSize', fsz);    
  end      
                    
 uicontrol('HorizontalAlignment', 'left', 'Parent',panel ,'Style', 'text','Units','centimeters','position',[ctrlleft 0.5 panelwidth/2 0.5],'string','Chromosome diameter (nm):', 'ForegroundColor', 'k', 'BackgroundColor', get(panel, 'backgroundcolor'), 'FontUnits', 'points', 'FontSize', fsz);    
 widthHndl = uicontrol('HorizontalAlignment', 'right', 'Parent',panel ,'Style', 'edit','Units','centimeters','position',[5 0.5 2 0.6],'string', '', 'ForegroundColor', 'k', 'BackgroundColor', 'w', 'FontUnits', 'points', 'FontSize', fsz);    
               
    
     
    % control buttons
    closeHndl = uicontrol(...
        'Style','pushbutton', ...
        'Units','centimeters', ...
        'Position',[figwidth-3 0.5 2.5 0.8], ...
        'Parent',plotfig, ...
        'string', 'Close', ...
        'FontUnits', 'points', 'FontSize', 10, ...
        'Callback',{@plotchromosome, 'close'});
    runHndl = uicontrol(...
        'Style','pushbutton', ...
        'Units','centimeters', ...
        'Position',[figwidth-5.75 0.5 2.5 0.8], ...
        'Parent',plotfig, ...
        'string', 'Plot', ...
        'FontUnits', 'points', 'FontSize', 10, ...
        'Callback',{@plotchromosome, 'run'});
        
    filedata = [];
    
    set(plotfig,'Visible','on');

elseif strcmp(action, 'run')

   %do the plots
   width=str2num(get(widthHndl, 'string'));
   if isempty(width) || isnan(width) || width <= 0
       ShowError('Please enter a positive width value', 'Bouquet plot');
       uicontrol(widthHndl);
       return;
   end
   
   %check which plots required
   flags = zeros(1,length(plotboxes));
   for i = 1:length(plotboxes)
       if get(plotboxes(i), 'value')
           flags(i) = 1;
       end
   end
   
   if isempty(filedata)
       ShowError('Please select a valid results file first', 'Plot Chromosome');
   elseif ~any(flags)
       ShowError('Please select at least one plot type', 'Plot Chromosome')
   else
       set(plotfig, 'pointer', 'watch');
       drawplots(width, filedata, flags);
       set(plotfig, 'pointer', 'arrow');
   end

elseif strcmp(action, 'close')
    
    delete(gcf);
    
elseif strcmp(action, 'selectfile')
    
    %browsing for a file

    title = 'Select a results file:';

    [FileName,PathName] = uigetfile('*.mat',title);
    
    
   if ~isequal(FileName, 0)
       %valid file selected
       filedata = [];
       set(fileHndl, 'String', []);
       emptydata = get(tblHndl, 'userdata');
       set(tblHndl, 'data', emptydata);
       file_selected = fullfile(PathName, FileName);
       %preview 
       set(plotfig, 'pointer', 'watch');
       drawnow;
       try
           filedata = load(file_selected);
       catch err
           % file read error
           ShowError('There was an error reading the selected file.', 'Plot Chromosome', err);
           set(plotfig, 'pointer', 'arrow');
           return;
       end
       %read ok, so fill in details, validate contents
       try
           contents = emptydata;    %just titles
           %fill in numbers
            
           contents{1, 2} = filedata.data.chrlen;
           contents{2, 2} = filedata.data.perlen;
           contents{3, 2} = filedata.data.phi;
           contents{4, 2} = filedata.data.numbeads;
           contents{5, 2} = filedata.data.cen;
           contents{6, 2} = filedata.data.k;
           contents{7, 2} = filedata.data.a;
           contents{8, 2} = filedata.data.kappa;
           contents{9, 2} = ['[' num2str(filedata.data.dir(1)) ' ' num2str(filedata.data.dir(2)) ' ' num2str(filedata.data.dir(3)) ']'];
           contents{10, 2} = filedata.data.numsteps;
           contents{11, 2} = filedata.data.numchain;
           contents{12, 2} = filedata.data.radius;
           contents{13, 2} = filedata.data.confpotential;
           contents{14, 2} = filedata.data.architecture;
           
           %cacluate default width
           width = -0.0056*filedata.data.phi^2 + 1.6*filedata.data.phi -3.4;
           

       catch
           ShowError('The selected file has an invalid format', 'Plot Chromosome');
           set(plotfig, 'pointer', 'arrow');
           filedata = [];
           return;
       end
       
       set(tblHndl, 'data', contents);
       set(fileHndl, 'String', file_selected);
       set(widthHndl, 'string', num2str(width));
       
      
       set(plotfig, 'pointer', 'arrow');
    end
end