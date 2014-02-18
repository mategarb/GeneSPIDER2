function varargout = boxplot(Data,varargin)
% <[options ,h, lgh]> = illustrate.boxplot(X,Data,<name>,<value>)
% function for plotting structured results from the Network Inference pipeline
%
% Data: data matrix, if dim=3 then a group for each 3rd dimension.
%
% name value pairs:
% ================
% 'xpos'      : x-coordinates for boxes in Data. If X = [], then
%               group data equidistant with groups (default behaviour)
% 'title'     : title  [str]
% 'x'         : xlabel [str]
% 'y'         : ylabel [str]
% 'legend'    : [str cell-array]
% 'c'         : colors for groups [str array] (or maybe it works with rgb values)
% 'highlight' : [range] mark x-range in plot to highlight
% 'width'     : [value] width of boxes.
% 
%
% output variables:
% ================
% options : all options that can be used as input
% h : figure handle (optional)
% lgh : legend handle (optional)
% 
% Examples:
% ========
% 
% Data(:,:,1) = rand(20,3);
% Data(:,:,2) = 0.5;
% Data(:,:,3) = randn(20,3);
% Data(:,:,4) = 1;
% 
% % Grouping the data
% illustrate.boxplot(Data,'legend',{'rand','0.5','randn','1'},)
% 
% % plotting all on specified positions X
% X = 20*randn(1,12);
% illustrate.boxplot(Data,'xpos',X)

options = struct('title','','x','x','y','y','c','rgbkcmy','width',1,'scale','lin','ylim',[]);
options(1).legend = {};
options(1).highlight = [];
options(1).xpos = [];
options(1).xtick = {' '};
options(1).xticklabel = {' '};
if nargin == 0
    if nargout == 0
        help illustrate.boxplot
    else
        varargout{1} = options;
    end
    warning(' No input data given')
    return
end
options = tools.optionParser(options,varargin{:});
X = options.xpos;
[n,m,d] = size(Data);
Data = reshape(Data,n,m*d);

g = length(X);

h = figure('color',[1,1,1]);
if g == 0
    g1 = repmat(1:m,1,d)';
    g2 = repmat(1:d,m,1);
    g2 = g2(:);
    if ischar(options.c)
        options.c = options.c(1:d);
    end
    boxplot(Data, {g1 g2}, 'factorgap', 10, 'symbol', '+', 'color', options.c)
    set(gca,'xtick',linspace(min(xlim)+d/2,max(xlim)-d/2,m))
    set(gca,'xticklabel',options.xticklabel)
    grid on
    set(gca,'XGrid','off')
else
    [X,D] = sort(X);
    boxplot(Data(:,D), 'positions', X, 'symbol', '+','width',options.width)

    if strcmp(' ',options.xticklabel)
        if strcmp('log',options.scale)
            set(gca,'xtick',logspace(min(log10(xlim)),max(log10(xlim)),8))
            set(gca,'XScale','log');
        else
            set(gca,'xtick',linspace(min(xlim),max(xlim),8))
        end
        xticklabel = strtrim(cellstr(num2str(get(gca,'xtick')','% 0.2g')));
        set(gca,'xticklabel',xticklabel)
    else
        set(gca,'xtick',options.xtick)
        set(gca,'xticklabel',options.xticklabel)
    end
    grid on
end

if ~isempty(options.ylim)
    ylim(options.ylim)
end

if ~isempty(options.highlight)
    hold on
    tmp = ones(1,2);
    plot(tmp*options.highlight(1),ylim,'k')
    plot(tmp*options.highlight(2),ylim,'k')
end


xlabel(options.x)
ylabel(options.y)
title(options.title)
names = fliplr(options.legend);
if ~isempty(options.legend)
    warning('off','MATLAB:legend:IgnoringExtraEntries')
    lgh = legend(findobj(gca,'Tag','Box'),names,'location','Best'); % The legend can be fixed. Check the handle and it's children
end


if nargout >= 1
    varargout{1} = options;
end

if nargout >= 2
    varargout{2} = h;
end

if nargout >= 3
    varargout{3} = lgh;
end

return
