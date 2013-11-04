function L = regressorLength(v,w)

if ~isrow(v)
    v = v';
end

if isrow(w)
    w = w';
end


L = v*w/(norm(v)*norm(w));