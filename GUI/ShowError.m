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


function  ShowError(msg, tl, err)

if nargin > 2
    %show system error message
     se = {err.message , [err.stack(1).file ' line ' num2str(err.stack(1).line)]};
else
    se = '';
end

uiwait(errordlg([msg se],tl,'modal'));
           