classdef Exchange < hgsetget
    
    methods
        function obj = Exchange()
        end
        
        function obj = populate(obj,input)
            inputnames = fieldnames(input);
            names = fieldnames(obj);
            for name = inputnames'
                if any(strmatch(name,names))
                    obj.(name{1}) = input.(name{1});
                end
            end
        end
        
        function save(obj,savepath,fending,varargin)
        % saves a datastruct obj to file, either as a .mat, .json, .ubj or .xml file.
        % The name of the file will be the name of the obj set.
        %
        %   Input Arguments: save(obj[,path,<ext>,<other_args>])
        %   ================
        %   (none) :        Save obj file to current directory as a .mat file.
        %   path :          Path to directory.
        %   ext :           File extension: .xml,  eXtensible Markup Language
        %                                   .json, JavaScript Object Notation
        %                                   .ubj,  Universal Binary JSON
        %                               or  .mat   matlab native format (default)
        %   other_args:     These arguments are sent to the specific save function
        %
        %
        % Support for writing json or ubj formats can be found here:
        %     https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files-in-matlab-octave
        %
        % Support for writing xml formats can be found here:
        %     https://www.mathworks.com/matlabcentral/fileexchange/6268-xml4mat-v2-0
        % a slightly optimized version is contained here:
        %    <my-repo>

            warning('off','MATLAB:structOnObject')
            if ~exist('fending','var')
                fending = '.mat';
            end
            if ~exist('savepath','var')
                savepath = './';
                if ~exist('savepath','dir')
                    error('save path does not exist')
                end
            end

            obj_data = struct(obj);

            if nargin == 3
                if ~isa(fending,'char')
                    error('Input arguments must be a string')
                end
                if strcmp(fending,'.xml')
                    if exist('mat2xml') ~= 2
                        error('Save method for xml files does not seem to exist')
                    end
                elseif strcmp(fending,'.json')
                    if exist('savejson') ~= 2
                        error('Save method for json files does not seem to exist')
                    end
                elseif strncmp(fending,'.ubj',4)
                    if exist('saveubjson') ~= 2
                        error('Save method for universal binary json files does not seem to exist')
                    end
                end
            end

            if isa(obj,'datastruct.Dataset')
                name = obj_data.dataset;
                savevar = 'obj_data';
            elseif isa(obj,'datastruct.Network')
                name = obj_data.network;
                savevar = 'obj_data';
            else 
                warning('Unknown object type')
                name = 'unknown_datatype';
                savevar = 'obj_data';
            end

            if strcmp(fending,'.xml')
                xmlString = simplify_mbml( spcharout( regexprep(mat2xml(obj_data,savevar),'\n',' ')) ); % spcharout can not handle newlines
                xmlwrite(fullfile(savepath,[name,fending]), str2DOMnode(xmlString));
            elseif strcmp(fending,'.json')
                savejson(savevar,obj_data,fullfile(savepath,[name,fending]),varargin{:});
            elseif strncmp(fending,'.ubj',4)
                saveubjson(savevar,obj_data,fullfile(savepath,[name,fending]),varargin{:});
            elseif strcmp(fending,'.mat')
                save(fullfile(savepath,name),savevar,varargin{:});
            else
                error('unknown file extension')
            end
        end
    end

    methods (Static = true)
        function varargout = load(varargin)
        % Load a dataset/network file back in to a datastruct.<Object>
        % obj_data = datastruct.<Dataset/Network>.load(['path/file'] or [path,file]);
        %
        %   Input Arguments: datastruct.<Dataset/Network>.loaddata([path,file])
        %   ================
        %   (none) :        Outputs a list of data files availible in the current directory.
        %   path :          Path to directory or full path with filename. If no file is specified
        %                   it will output a list of availible data sets in that directory.
        %   file :          Filename or number of its place in the list of files.
        %                   Have to include the path input variable.
        %
        %   Output Arguments: obj / list
        %   ================
        %   obj :           Populate the datastruct.<Object> with the loaded file.
        %   list :          If no file is specified a list of availible datasets
        %                   is returned for the given directory.
        %

            lpath = pwd;
            lfile = [];
            if nargin == 1
                if isa(varargin{1},'double')
                    lfile = varargin{1};
                else
                    if exist(varargin{1}) == 2
                        [p,f,e] = fileparts(varargin{1});
                        lpath = p;
                        lfile = [f,e];
                    elseif exist(varargin{1}) == 7
                        lpath = varargin{1};
                    else
                        error('Unknown path or file')
                    end
                end
            elseif nargin == 2
                lpath = varargin{1};
                if exist(lpath) ~= 7
                    error('Unknown path: %s',lpath)
                end
                if isa(varargin{2},'double')
                    lfile = varargin{2};
                else
                    if exist(fullfile(lpath,varargin{2})) == 2
                        lfile = varargin{2};
                    else
                        error('Unknown file')
                    end
                end
            elseif nargin > 2
                error('wrong number of input arguments.')
            end

            if isa(lfile,'double')
                obj_datas = dir(lpath);
                j = 0;
                for i = 1:length(obj_datas)
                    if ~obj_datas(i).isdir
                        [pa,fi,fext] = fileparts(obj_datas(i).name);
                        if strcmp(fext,'.mat') || strcmp(fext,'.xml') || strcmp(fext,'.json') || strncmp(fext,'.ubj',4)
                            j = j+1;
                            output{j} = obj_datas(i).name;
                        end
                    end
                end
                if ~exist('output','var')
                    warning('no data files in directory:\n %s',lpath)
                    return
                end
                if nargout == 1 && isempty(lfile)
                    varargout{1} = output;
                    return
                elseif isempty(lfile)
                    if nargout == 2
                        error('To many output arguments defined. Only file list can be returned')
                    end
                    for j=1:length(output)
                        fprintf('%d %s\n',j,output{j});
                    end
                    return
                else
                    lfile = output{lfile};
                end
            end
           
            fetchfile = fullfile(lpath,lfile);
            [p,f,e] = fileparts(fetchfile);
            if strcmp(e,'.mat')
                obj_data = load(fetchfile);
                % For backwards compatability.
                if isfield(obj_data,'dataset')
                    obj_data = obj_data.dataset;
                elseif isfield(obj_data,'network')
                    obj_data = obj_data.network;
                end
            elseif strcmp(e,'.json')
                obj_data = loadjson(fetchfile);
                name = fieldnames(obj_data);
                obj_data = obj_data.(name{1});
            elseif strncmp(e,'.ubj',4)
                obj_data = loadubjson(fetchfile);
                name = fieldnames(obj_data);
                obj_data = obj_data.(name);
            elseif strcmp(e,'.xml')
                [MAT,obj_data] = xml2mat(fetchfile);
                eval([obj_data,'=MAT;']);
            end
            
            
            % regexprep(name,'(\<[a-z])','${upper($1)}')
            % Ugly hack if statement
            if isfield(obj_data,'Y') % is a dataset
                obj = datastruct.Dataset;
            elseif isfield(obj_data,'A')
                obj = datastruct.Network; % is a network
            else
                error('no compatible class for the data file')
            end
            
            populate(obj,obj_data);
            if nargout == 1
                varargout{1} = obj;
            elseif nargout == 2
                varargout{1} = obj;
                varargout{2} = obj_data;
            end        
        end
    end
end
