function s = UpdateProgress(varargin)


global BOUQUET_STOP
global barHndl num_div
persistent tm running

if isempty(BOUQUET_STOP)
    %not using GUI
    s = 0;
    return;
end
% increments progress bar
data = get(barHndl, 'Userdata');

%length of time that bar takes to advance one increment, in seconds
%This should be a guess at the interval between successive increments
interval = data(1);
% total number of large increments
barinc = data(2);

barframe =get(barHndl, 'parent');
maxWidth = get(barframe, 'position');
maxWidth = maxWidth(3);
%num_div is number of small steps in one large increment. Scale so each is approx 0.1 mm on
%screen. divide width of full bar by the total number of large increments
num_div = ceil(maxWidth*barinc*10);

if strcmp(varargin{1}, 'start')
    %initialise bar
    %remove precision beyond 2dp to avoid warning
    per = sprintf('%5.2f', interval/num_div);
    per = str2num(per);
    tm = timer('TimerFcn', @timer_update, 'ExecutionMode', 'fixedRate', 'Period', per, 'TasksToExecute', num_div);
    start(tm);
    running = true;
    
elseif strcmp(varargin{1}, 'stop')
    %stop bar
    if running
        stop(tm);
        delete(tm);
        running=false;
    end
    barlen = [0 0 1 1];
    set(barHndl, 'Position', barlen);
    drawnow;
else
    %increment bar
    if strcmp(get(tm, 'Running'), 'on')
        stop(tm);
        num = get(tm, 'TasksExecuted');
        for i = 1:num_div-num
            timer_update(tm, [])
        end
    end
    start(tm);
    printStr(varargin);
end

s = BOUQUET_STOP;


%==========================================================================

function timer_update(tm, data)

global barHndl num_div

barlen = get(barHndl, 'Position');
barinc = get(barHndl, 'Userdata');
barinc = barinc(2);

barlen(3) = min(barlen(3) + barinc/num_div, 1);
set(barHndl, 'Position', barlen);
drawnow;
 
