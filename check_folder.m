function [] = check_folder(varargin)
    for i = 1:nargin
        folder = varargin{i};
        if exist(folder) == 0
            mkdir(folder);
        end
    end
end