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


function drawplots(wd, filedata, flags)

%file is one output by bouquetgui
%flags is which to plot, a vector of length 6

global plot_titles plot_funcs Loci Width

Loci = [3,10,15]';
Width = wd;

%get the data
try
    net = filedata.net; %cell array
    
    %choose positions
    set(0,'Units','centimeters')
    screen_size = get(0,'ScreenSize');
    swidth = screen_size(3);
    sheight = screen_size(4);
    pos = cell(6, 1);
    
    switch sum(flags)
        
        case 1
            %just one so fill screen
            pos{1} = [3 3 swidth-6 sheight-6];
        case 2
            pos{1} = [0 sheight-(swidth/2)-2 swidth/2 swidth/2];
            pos{2} = [swidth/2 sheight-(swidth/2)-2 swidth/2 swidth/2];
        case 3
            pos{1} = [0 sheight-(swidth/3)-2 swidth/3 swidth/3];
            pos{2} = [swidth/3 sheight-(swidth/3)-2 swidth/3 swidth/3];
            pos{3} = [2*swidth/3 sheight-(swidth/3)-2 swidth/3 swidth/3];
        case 4
            pos{1} = [0 sheight/2-1 sheight/2 sheight/2-1];
            pos{2} = [sheight/2-1 sheight/2-1 sheight/2 sheight/2-1];
            pos{3} = [0 0 sheight/2 sheight/2-1];
            pos{4} = [sheight/2-1 0 sheight/2 sheight/2-1];
        otherwise
            pos{1} = [0             sheight-swidth/3-1     swidth/3   swidth/3-1];
            pos{2} = [swidth/3-1   sheight-swidth/3-1     swidth/3   swidth/3-1];
            pos{3} = [2*swidth/3-1 sheight-swidth/3-1     swidth/3   swidth/3-1];
            pos{4} = [0             sheight-2*swidth/3-1     swidth/3   swidth/3-1];
            pos{5} = [swidth/3-1   sheight-2*swidth/3-1     swidth/3   swidth/3-1];
            pos{6} = [2*swidth/3-1 sheight-2*swidth/3-1     swidth/3   swidth/3-1]; 
       
    end
    
    num = 1;
    for i = 1:length(flags)
        if flags(i)
            figure('numbertitle', 'off', 'name', plot_titles{i}, 'units', 'centimeters', 'position', pos{num});
            feval(plot_funcs{i}, net); %width needed by 1, loci by 2 and 3
            title(plot_titles{i});
            num = num + 1;
        end
    end

catch err
    ShowError('The selected file has an invalid format', 'Meiosis Analysis', err);
end
    
    
    
%     if flags(2);
%         %Fit a beta distribution to a set of loci of choice, in this case loci 3,
%         %10 and 15
%         [val] = RDF(net1,Loci);
%         x     = linspace(0,1,100);
%         figure('numbertitle', 'off', 'name', plot_titles{2}, 'units', 'centimeters', 'position', pos{num})
%         for j = 1:size(val,2)
%             Y = betapdf(x,val{j}(1,1),val{j}(1,2));
%             hold on
%             plot(x,Y,'-')
%         end
%         title(plot_titles{2});
%         num = num +1;
%     end
%     
%     if flags(3)
%         %Calclate the physical distances between two loci (or two sets of loci) and
%         %fit a Beta distribution (parameters a and b). In this case Loci 3, 10 and
%         %15 on the first chromosome and 3, 10 and 15 on the second.
%         x     = linspace(0,1,100);
%         [val] = LocusLocusDistance(net1,net2,Loci,Loci);
%         figure('numbertitle', 'off', 'name', plot_titles{3}, 'units', 'centimeters', 'position', pos{num});
%         for j = 1:size(val,2)
%             Y = betapdf(x,val{j}(1,1),val{j}(1,2));
%             hold on
%             plot(x,Y,'-');
%         end
%         title(plot_titles{3});
%         num = num +1;
%     end
    
 




