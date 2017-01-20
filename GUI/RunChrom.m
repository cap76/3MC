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
function RunChrom(varargin)

%Runs the analysis and displays a progress bar

persistent data runfig closeHndl

global txtHndl barHndl barPanelHndl BOUQUET_STOP

if ischar(varargin{1})
    %called in code
    action = varargin{1};
else
    %called from a callback like this RunBouquet(src, event, action, ....);
   action = varargin{3}; 
end


if strcmp(action, 'go')
    %creates controls on the first panel
    data = varargin{2};
    mainfig = varargin{3};
    
    %generate new figure, centered on main figure
    runfig = figure('WindowStyle', 'modal', 'pointer', 'watch', 'resize', 'off', 'Menubar', 'none' ,'NumberTitle','off','Visible','off', 'CloseRequestFcn', @noclose, 'name', data.architecture);

    
    figwidth = 10;
    figheight = 4;
    mainpos = get(mainfig, 'position');
    pos = [mainpos(1)+(mainpos(3)-figwidth)/2 mainpos(2)+(mainpos(4)-figheight)/2 figwidth figheight];
    set(runfig, 'Units', 'centimeters', 'Position', pos);
    
    maincol = get(runfig, 'Color');
    
    %messages
    txtHndl = uicontrol('HorizontalAlignment', 'left', 'Parent',runfig ,'Style', 'text','Units','centimeters','position',[0.5 figheight-1.25 figwidth-1 1],'string','Running Analysis...', 'ForegroundColor', 'k', 'BackgroundColor', maincol, 'FontUnits', 'points', 'FontSize', 10);     
                
    %progress bar

    barPanelHndl=uipanel('Units','centimeters', 'Position', [0.5 figheight-2 figwidth-1 0.5],'Parent',runfig, 'BackgroundColor',  maincol);
    barHndl = uicontrol('visible', 'off', 'style', 'text', 'Units','normalized', 'Position', [0 0 0.001 1],'Parent',barPanelHndl, 'BackgroundColor', 'r', 'ForegroundColor', 'r' );
   
    %cancel button
    closeHndl = uicontrol(...
        'Style','pushbutton', ...
        'Units','centimeters', ...
        'Position',[figwidth/2-1 0.5 2 0.8], ...
        'Parent',runfig, ...
        'string', 'Cancel', ...
        'FontUnits', 'points', 'FontSize', 10, ...
        'Callback',{@RunChrom, 'close'});
    
    set(runfig, 'visible', 'on');
    
    BOUQUET_STOP= false;
 
    interval = 8;   %guess at interval between successive increments of progress bar (secs)
    n = data.numchain * 22;  %number of times bar advanced
    barInc = 1 / n; 

    set(barHndl, 'Userdata', [interval barInc], 'visible', 'on', 'position', [0 0 0.001 1]);
      
   
    UpdateProgress('start');
    [success err] = Run(data);
    UpdateProgress('stop');
   
    
    if ~success       
        if BOUQUET_STOP
            %stopped by user
            set(txtHndl, 'String', 'Analysis cancelled');
        else
            %error
            set(txtHndl, 'String', 'Error running the analysis');
            ShowError('There was an error running the analysis', 'Error', err);
        end
    end
    set(closeHndl, 'string', 'Close');
    set(runfig, 'pointer', 'arrow');
    
elseif strcmp(action, 'close')
    
    lbl = get(gcbo, 'string');
    if strcmp(lbl, 'Cancel')
       %cancel analysis
        set(txtHndl, 'String', 'Stopping chromosome analysis...');
        drawnow;
        BOUQUET_STOP = true;
    elseif strcmp(lbl, 'Close')
        delete(runfig);
    end

end

%==========================================================================
function [success err] = Run(data)

global BOUQUET_STOP

success = 0;
err = [];

try
   % VMF    = [data.kappa,0,0,1;data.kappa,0,0,1];          %von Mises-Fisher parameters.
   % Note: the VMF is a 1 x 4 vector in Rabl, a 2x4 matrix in Bouquet and a
   % 3x4 matrix in Rabl2Bouquet. For the time being we simply replicate the
   % 1x4 matrix, however, the other entries should also be included as
   % inputs
   VMF      = [data.kappa,data.dir(1),data.dir(2),data.dir(3);data.kappa,data.dir(1),data.dir(2),data.dir(3);data.kappa,data.dir(1),data.dir(2),data.dir(3)];
    a       = data.a;
    radius  = data.radius;
    len     = data.numbeads;
    cen     = data.cen;
    
    %NOTE: cen needs to be input as a free parameter (it represents the position of the centromere, hence must be between 1 and len). 
    %Architecture is the choice of scripts to be run and should be input as
    %e.g., checkboxed.
    %cen     =  data.cen;%66;% data.Cen; 
    architecture2 = data.architecture2;
    architecture = data.architecture;
    numchain = data.numchain;%3;
        
    printStr('Initialising...');
    %save params first to validate file name
    save(data.filename, 'data');
    printStr('Done');    
    printStr('Generating random flight...');
    
    
    
    switch upper(architecture)
        
        case 'TETHERED TEL'
            [x0,y0,z0,theta,phi] = tMultivMF(data.dir(1),data.dir(2),data.dir(3),100,radius,1);% last para means ONE LOOP INTERATION
            [x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);
        case 'TETHERED CEN'            
            [x0,y0,z0,theta,phi] = tMultivMF(data.dir(1),data.dir(2),data.dir(3),100,radius,1);% last para means ONE LOOP INTERATION            
            [x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);
            dx = x1(1,cen); dy = y1(1,cen); dz = z1(1,cen);
            x1 = x1 - dx; y1 = y1 - dy; z1 = z1 - dz + radius;
        case 'TETHERED TEL/CEN'
            [x0,y0,z0,theta,phi]    = tMultivMF(0,0,1,100,radius,1);
            [x01,y01,z01,theta,phi] = tMultivMF(0,0,1,100,radius,1);
            [x02,y02,z02,theta,phi] = tMultivMF(0,0,1,100,radius,1);
            [x11,y11,z11] = RandFlightCon([x0,y0,z0],[x01,y01,z01],a,len-cen);
            [x22,y22,z22] = RandFlightCon([x01,y01,z01],[x02,y02,z02],a,cen+1);
            x1 = [x11,x22(1,2:length(x22))];
            y1 = [y11,y22(1,2:length(x22))];
            z1 = [z11,z22(1,2:length(x22))];            
        otherwise
            [x0,y0,z0,theta,phi] = tMultivMF(data.dir(1),data.dir(2),data.dir(3),100,radius,1);% last para means ONE LOOP INTERATION
            [x1,y1,z1] = RandFlightCon([x0,y0,z0],[x0,y0,z0],a,len);         
    end
    

    printStr('Done');
    if BOUQUET_STOP
        return;
    end
    %vary number to be run, in parallel too
    switch upper(architecture2)
        case 'GAUSSIAN PROCESS' %Gen \kappa from GP
    
                try %Try loading in initial GP hyperpramts
                    load GPhypinit.mat
                    H         = GPhypinit.H; %Hyperparameters of covariance function
                    K         = feval(GPhypinit.covfunc{:}, H, linspace(0,1,len-2)'); %Covariance matrix at observation points
                    kappa_val = exp(real(gsamp(zeros(1,len-2)',K,1)) -1);                    
                catch %Otherwise initialise a new GP
                    covfunc   = {'covSum',{'covSEiso','covNoise'}}; 
                    H         = [log(.1); log(0.3); log(1e-5)]; %Hyperparameters of covariance function
                    K         = feval(covfunc{:}, H, linspace(0,1,len-2)'); %Covariance matrix at observation points
                    kappa_val = exp(real(gsamp(zeros(1,len-2)',K,1)) -1);                    
                end                
        otherwise
                load kappa.mat                
    end

    %Nodimensionalise the persistence-length
    %kappa_val = kappa_val*data.chrlen*0.34/data.phi;
    
    parameters = [data.numsteps, len, radius, data.confpotential, a, cen, kappa_val];
    net = cell(1,numchain);
    

    switch upper(architecture)
        case 'TETHERED TEL'
                archfunc = @BouquetTV;
        case 'TETHERED CEN'  
                archfunc = @Rabl;
        case 'TETHERED TEL/CEN'
                archfunc = @Rabl2Bouquet;
        otherwise
                archfunc = @FreeWLCTV;
    end
    
   for chainloop = 1:numchain  
       if UpdateProgress('Running Analysis...', [architecture ' ' num2str(chainloop) ' of ' num2str(numchain)]);
            break;  %stop clicked
       end
       [net{1,chainloop}] = feval(archfunc, parameters,x1,y1,z1,VMF);
   end
    
        
    if BOUQUET_STOP
        return;
    else
        printStr({'Saving results...', [data.architecture ' completed']});
        save(data.filename, 'net', '-append');
        printStr(['Analysis complete. Results saved at ' data.filename]);
    end
   
catch err
    printStr(['There was an error running ' architecture]);
    return;
end
success = 1;

%==========================================================================

function noclose(src,evnt)

%do nothing prevent user closing figure from menubar
return;



