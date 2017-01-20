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

function t = timeStr(tend)

%converts a number of seconds to a string repesenting hrs, min, sec

if tend < 60
    t = [num2str(round(tend)) ' seconds'];
elseif tend < 3600
    minutes = floor(tend/60);
    seconds = round(tend - minutes*60);
    t = [num2str(minutes) ' minutes ' num2str(seconds) ' seconds'];
else
    hours = floor(tend/3600);
    minutes = floor((tend - hours * 3600)/60);
    seconds = round(tend - (hours * 3600 + minutes*60));
    t = [num2str(hours) ' hours ' num2str(minutes) ' minutes ' num2str(seconds) ' seconds'];
end