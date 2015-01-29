function L = regressorLength(v,w)
% helper function for calculating,
% \[
%   L = \frac{\boldsymbol{v}^T \boldsymbol{w}}{||\boldsymbol{v}|| \cdot ||\boldsymbol{w}||}
% \]
% L = regressorLength(v,w)
%
%

if ~isrow(v)
    v = v';
end

if isrow(w)
    w = w';
end


L = v*w/(norm(v)*norm(w));