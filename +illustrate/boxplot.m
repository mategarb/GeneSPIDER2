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

options = struct('title','this is a title','x','x','y','y','c','rgbkcmy','width',1);
options(1).legend = {};
options(1).highlight = [];
options(1).xpos = [];
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

g1 = repmat(1:m,1,d)';
g2 = repmat(1:d,m,1);
g2 = g2(:);
g = length(X);

options.c = options.c(1:d);

h = figure('color',[1,1,1]);
if g == 0
    boxplot(Data, {g1 g2}, 'factorgap', 10, 'symbol', '+', 'color', options.c)
    % set(gca,'xticklabel',{' '})
    set(gca,'xtick',linspace(min(xlim)+d/2,max(xlim)-d/2,m))
    % set(gca,'xticklabel',linspace(min(xlim)+d/2,max(xlim)-d/2,m))
    set(gca,'xticklabel',options.xticklabel)
    grid on
    set(gca,'XGrid','off')
else
    boxplot(Data, X, 'positions', X, 'symbol', '+','width',options.width)
    if strcmp(' ',options.xticklabel)
        set(gca,'xtick',linspace(min(xlim),max(xlim),8))
        set(gca,'xticklabel',get(gca,'xtick'))
    else
        set(gca,'xtick',options.xticklabel)
        set(gca,'xticklabel',options.xticklabel)
    end
    grid on
end

if ~isempty(options.highlight)
    hold on
    tmp = ones(1,2);
    plot(tmp*options.highlight(1),ylim,'k')
    plot(tmp*options.highlight(2),ylim,'k')
end


% if ~isempty(options.labels) % rearange xticks to fit to each group
%     set(gca,'XTickLabel',{' '})
% end

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
% , 'positions', X


% stuff = get(findobj(gca,'Type','text'))
% t = get(stuff(1).Parent,'Children') % Gets all
% now look for the Type
% get(t(n))
% tmp = get(t(40))
% tmp.Type should be text

% text objects
% nanan = get(findobj(gca,'Type','text'))
% This could probably be done in an easier way
% tmp = get(findobj(gca,'Type','text'));
% t = get(tmp(1).Parent,'Children');
% tmp = get(t(40));
% tt = get(findobj(get(get(gca,'Children'),'Children'),'Type','text'))
textobj = findobj(get(get(gca,'Children'),'Children'),'Type','text');
for i=1:length(textobj)
    xtext(i,:) = get(textobj(i),'Position');
end

% h = xtext;
h = textobj;
