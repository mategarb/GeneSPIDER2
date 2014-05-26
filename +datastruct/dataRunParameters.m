classdef dataRunParameters
% dataRunParameters
% parameters = struct('alpha',0.05, 'cls','fcls', 'variable','SNRm', 'range',logspace(-1,3,5), 'zetavec',logspace(-6,0,100), 'looco','false', 'nit',1);
%

    properties (SetAccess = public)
          alpha = 0.05;
            cls = 'fcls';
       variable = 'SNRm';
          range = logspace(-1,3,5);
        zetavec = logspace(-6,0,100);
          looco = false;
            nit = 1;
         method = '';
        dataset = '';
        network = '';
    end
end
