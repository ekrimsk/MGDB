classdef mgdb < matlab.mixin.Copyable 
% MGDB (Motor-Gearbox DataBase)   
%
%  MGDB Methods:
%       get_combinations - gets a list of motor keys and gearbox keys 
%               corresponding to compatible motor/gearbox combinations 
%       update_filters - updates selection criteria for determining which 
%                       combinations should be returned by get_combinations
%       clear_filters - sets the selection filters back to defaults 
%       update_settings - updates mgdb print settings 
%
%  MGDB Properties:
%       motors -  a containers.Map of motors 
%       gearboxes - a containers.Map of gearboxes
%       compatibility - a map of maps given the 'score' of each motor/gearbox
%                       combination (initialized to inf)
%       filters - see update_filters
%       settings - struct, settings.verbose sets the print level 
%       dd_key - the key used to refer to a direct drive dummy gearbox
%
%   For a list of motor and gearbox properties see the readme at 
%   https://github.com/ekrimsk/MGDB
%
%   TODO -- copy it here as well?
%
%   Author: 
%       Erez Krimsky, ekrimsky@stanford.edu, 9/9/20
%       Stanford University, Biomechatronics Lab 
%
% See also CORDS
properties (GetAccess = public, SetAccess = private)
    motors          % containers.Map
    gearboxes       % containers.Map
    compatibility 
    settings 
    filters 
    dd_key   
end 


methods (Access = public)


    function obj = mgdb(varargin) %
        % TODO -- add onputs to change database path 

        init_tic = tic; 

        if isempty(varargin)
            % TODO -- robustify 
            database_paths = find_dirs('mgdb'); % look for a directory called mgdb
        else 
            database_paths = find_dirs(varargin{1});
        end 

        obj.dd_key = 'DD'; 
        settings.verbose = 1;   % this is the only settings write now 
        obj.settings = settings; 
        % TODO --- accouint for input redirecting us to a new search path 
        % If no args look in a a default location
        
        % Could look in other places too i guess
        % cell array of paths 

        % motor files called "*_motors.csv"
        % gearbox files called "*_gearboxes.csv"
        % compatibility filed called "*_compatibility.csv"
        % see somewhere else (TODO) for formatting 
        for i = 1:length(database_paths)
            pth = database_paths{i};
            % NOTE -- will this work on windows 
            % TODO -- double checj and maybe rewrite with "fullfile" for 
            % windows compatibility 
            motor_file_tmp = dir(sprintf('%s/**/*_motors.csv', pth));
            gearbox_file_tmp = dir(sprintf('%s/**/*_gearboxes.csv', pth));
            compat_file_tmp = dir(sprintf('%s/**/*_compatibility.csv', pth));

            m_tmp = cellfun(@(x) fullfile(x.folder, x.name),...
                            num2cell(motor_file_tmp), 'UniformOutput', false);
            g_tmp = cellfun(@(x) fullfile(x.folder, x.name),...
                         num2cell(gearbox_file_tmp), 'UniformOutput', false);
            c_tmp = cellfun(@(x) fullfile(x.folder, x.name),...
                           num2cell(compat_file_tmp), 'UniformOutput', false);
            if i == 1
                motor_files = m_tmp;
                gearbox_files = g_tmp;
                compat_files = c_tmp; 
            else 
                motor_files = [motor_files, m_tmp];
                gearbox_files = [gearbox_files, g_tmp];
                compat_files = [compat_files, c_tmp];
            end
        end 

        % Read in motor csvs, NOTE move to own function takinfg in list 
        [motor_keys, motor_values] = add_motors(motor_files); 
        % Read in Gearbox csvs --NOTE: may move to its own function taking in list 
        [gb_keys, gb_values] = add_gearboxes(gearbox_files); 

        % NOTE: A motor might appear mulitpe times in multiple compatibility 
        % files (could be compatible with multiple manufactureres for example)
        obj.motors = containers.Map(motor_keys, motor_values);
        obj.gearboxes = containers.Map(gb_keys, gb_values);
        % add the 'dummy' gearbox for a direct drive option (compatible with all motors)
        

        obj.gearboxes('DD') = direct_drive(); % add direct drive option to the map 

        obj.compatibility = containers.Map('KeyType', 'char', 'ValueType', 'any'); % value type will be other maps 
        % Read in compatibility csvs 
        add_compatibility(obj, compat_files); % TODO -- all files 

        %
        %
        %      Initialize Filters 
        %
        %  filters is a struct where each field is different criteria
        %  initialized to not filter out any possible choices 
        % 
        obj.filters = default_filters;
        %
        %           Save Object Out So Can Be Loaded by motor selection 
        %   




        init_time = toc(init_tic); 
        fprintf('Added %d motors and %d gearboxes to database in %0.2f seconds\n',...
                    obj.motors.Count, obj.gearboxes.Count, init_time); 

    end 



    function [motor_combo_keys, gearbox_combo_keys] = get_combinations(obj)
    %
    %
    %   Return structs instead of keys. This is less space efficient but it 
    %   it makes it easier for the caller to use the outputs 
    %
    %

        % get motor, gearbox key tuples 
        motor_keys = obj.motors.keys; 

        motor_combo_keys = {}; % will be a cell array of motor keys 
        gearbox_combo_keys = {}; % will be a cell array of gb keys 

        % NOTE: may want to allocate pretty big then trim 
        cost_list = {}; % will be used to reindex and sort outputs 

        init_combos = 0;
        % For each motor, get all compatible gearboxes 
        for m_idx = 1:length(motor_keys)
            motor_key = motor_keys{m_idx}; 
            tmp_map = obj.compatibility(motor_key);
            gb_keys = tmp_map.keys; 
            for gb_idx = 1:length(gb_keys) 
                gb_key = gb_keys{gb_idx};
                % Check if combo is valid given the filters
                if obj.valid_combo(motor_key, gb_key)
                    motor_combo_keys{end + 1} = motor_key;
                    gearbox_combo_keys{end + 1} = gb_key; 
                    cost_list{end + 1} = tmp_map(gb_key);
                end 
                init_combos = init_combos + 1; % to keep track of how many removed 
            end 
        end 

        num_combos = length(motor_combo_keys);
        num_removed = init_combos - num_combos;
        % TODO - verb
        fprintf('Found %d initial combinations, removed %d\n', init_combos, num_removed);

        % TODO -- some indicator of how many combos filtering removed
        % will want verbosity settings --> pass in from caller
        costs = cell2mat(cost_list); 
        [~, reindex] = sort(costs);

        % Sort by rankings -- reindex keys much faster than reindexing structs 
        motor_combo_keys = motor_combo_keys(reindex);
        gearbox_combo_keys = gearbox_combo_keys(reindex); 
    end 

    function update_rankings(obj, motor_keys, gearbox_keys, new_costs)

        % check that motor_keys, gearbox_keys, new_costs are same length 

        % support string arrays seperate from cell arrays of chars?? 
        for i = 1:length(motor_keys)
            motor_key = motor_keys{i};
            if obj.compatibility.isKey(motor_key)
                tmp_map = obj.compatibility(motor_key);  % by reference 
                tmp_map(gearbox_keys{i}) = new_costs(i);
            else 
                error('Invalid motor key');
            end 
        end 
    end 

    function clear_filter(obj, varargin) 
    %
    %
    %
        % with no extra args, clears all filters 

        % with args -- only clears an individual filter 

        % TODO -- come up with filter list 
        if isempty(varargin{:})  % clear all 

        else 

        end 


    end 


    function update_filters(obj, new_filters)  % TODO switch to struct base 
    % update_filters
    %   new_filters - struct (TODO )
    %
    %
        fields = fieldnames(new_filters); 
        for ii = 1:numel(fields)
            criteria = fields{ii};
            value = new_filters.(criteria);
            if ~iscell(value)
                assert(value >= 0, 'update_filter: values must be nonnegative'); 
            end 
            if isfield(obj.filters, criteria)
                obj.filters.(fields{ii}) = value;
            else 
                error('update_filter: invalid filter criteria');
            end 
        end 
    end 

    function update_settings(obj, new_settings)
        obj.settings.verbose = new_settings.verbose; 
    end 

    function filters = get_filters(obj)
        filters = obj.filters;
    end 
end 


methods (Access = private)

    function [valid] = valid_combo(obj, motor_key, gearbox_key)
    %
    %
    %
    %
    %
        % Check this combination of motor and gearbox against the filters 
        motor = obj.motors(motor_key);
        gearbox = obj.gearboxes(gearbox_key);

        valid_manufs = obj.filters.manufacturer;

        % Maximum speed 
        valid = true;
        if (obj.filters.omega_max*gearbox.ratio) > motor.max_int_speed
            valid = false; 
        elseif (obj.filters.omega_rms*gearbox.ratio) > motor.max_cont_speed
            valid = false; 
        elseif obj.filters.tau_max > gearbox.max_int_torque
            valid = false;
        elseif obj.filters.tau_rms > gearbox.max_cont_torque
            valid = false;
        elseif motor.mass + gearbox.mass > obj.filters.total_mass 
            valid = false; 
        elseif (motor.inertia*gearbox.inertia)*(gearbox.ratio^2) ...
                                    > obj.filters.effective_inertia
            valid = false;
        elseif ~isempty(valid_manufs)  % if empty than any is fine 
            try 
                % will throw an error if manuf string is now good 
                validatestring(lower(motor.manufacturer), lower(valid_manufs));
                validatestring(lower(gearbox.manufacturer), lower(valid_manufs));
            catch ME
                if strcmp(ME.identifier,'MATLAB:unrecognizedStringChoice')
                    valid = false;
                else 
                    rethrow(ME);
                end 
            end 
        end     
    end 


    function add_compatibility(obj, compat_files)
    %
    %
    %
    %
    %
        % could combine all compat files together to make simpler 
        % assume this is done somewhere above this line of code 

        gearbox_keys = obj.gearboxes.keys; % WILL RETURN THE KEYS SORTED ALPHABETICALLY 

        % TODO -- check that compat files is cell array even if only 1 file 

        for kk = 1:length(compat_files)
            compat_file = compat_files{kk}; 
            compat_fid = fopen(compat_file);

            tline = fgetl(compat_fid);
            while ischar(tline)     % go through whole file
                line_cells = strsplit(tline, ',');  % might be no compat 
                motor_key = line_cells{1};

                % First check that the motor key is actually present 

                % For each motor add direct drive compatibility 

                % Now get gearboxes 

                % Loop through gearbox keys in the row
                % As a O(nlogn) compromise heres the plan
                % get the full list of gearbox keys O(n) 
                % Sort the list of gearbox keys lexocigraphically O(nlogn)
                % If test does NOT contain a wildcard, index directly into map
                % If test key DOES contain a wildcard, find the corresponding keys 
                % Since sorted, just need to find the first (lexicographicaly) wildcard match 
                % and the last wildcard match and then everything in between is a match 
                for i = 2:length(line_cells)
                    gb_key = line_cells{i}; % may be partial key 
                    if obj.gearboxes.isKey(gb_key)
                        new_gb_match_keys = {gb_key};  % assumed cell array 
                    elseif endsWith(gb_key, '*')        % partial key
                        new_gb_match_keys = find_matches(gearbox_keys, gb_key(1:end-1)); % TODO -- can rewrite a faster subroutine 
                    else % error 
                        warning('Gearbox %s not found in database', gb_key);
                    end 

                    % TODO -- if gb_match_keys is isempty? 
                    if i == 2
                        gb_match_keys = new_gb_match_keys;
                    else % append 
                        % could be more efficient to grow cell
                        gb_match_keys = [gb_match_keys, new_gb_match_keys];
                    end 

                    if isempty(gb_match_keys)
                        warning('No matching gearboxes found for %s', motor_key);
                    end 
                end

                if obj.compatibility.isKey(motor_key)
                    % already have things added, append 
                    for j = 1:length(gb_match_keys)
                        gb_key = gb_match_keys{j};
                        tmp_map = obj.compatibility(motor_key);     % should be copy by ref
                        tmp_map(gb_key) = inf; % initialize all combination to have infinite cost 
                        % Could check if already key and issue warning if redundant
                    end 
                else 
                    % Add all the new ones together AND a direct drive option
                    % It is much faster to add lots of things to the map
                    % at once so it does not need to get resized and copied 
                    inf_cell = num2cell(inf(1, length(gb_match_keys) + 1)); 
                    obj.compatibility(motor_key) = containers.Map([obj.dd_key, gb_match_keys], inf_cell); 
                end 


                tline = fgetl(compat_fid);
            end % finish reading through lines of this file 
            fclose(compat_fid); 
        end % end looping through files 

        % Now add direct drive for any motors that were skipped 
        motor_keys = obj.motors.keys;
        for m_idx = 1:length(motor_keys)
            motor_key = motor_keys{m_idx};
            if ~obj.compatibility.isKey(motor_key)
                obj.compatibility(motor_key) = containers.Map(obj.dd_key, inf); 
            end 
        end 
    end 


    function vprintf(obj, priority, varargin)
    % vprintf - verbose printf 
    %   
    %   if 'priority' >= obj.settings.verbose, 
    %   varargin{:} gets passed to printf
    %
    %    
        if priority >= obj.settings.verbose
            fprintf(varargin{:}); 
        end 
    end 


end % end private methods 

%=================== END CLASS DEF =========================
end 


function dirs = find_dirs(tofind)
% copied verbatim from https://www.mathworks.com/matlabcentral/answers/347892-get-full-path-of-directory-that-is-on-matlab-search-path
% shoutout to Walter Robinson
    esctofind = regexptranslate('escape', tofind);   %in case it has special characters
    dirs = regexp(path, pathsep,'split');          %cell of all individual paths
    temp = unique(cellfun(@(P) strjoin(P(1:find(strcmp(esctofind, P),1,...
                'last')),filesep), regexp(dirs,filesep,'split'), 'uniform', 0));    
    dirs = temp(~cellfun(@isempty,temp));     %non-empty results only
end 

function dd_struct = direct_drive()
    dd_struct = struct('key', 'DD', 'manufacturer', 'none',...
                         'ID', 'DD', 'type', 'none',...
                        'stages', 0, 'ratio', 1, 'mass', 0, 'inertia', 0,...
                        'efficiency', 1, 'direction', 1,...
                        'max_int_torque', inf, 'max_cont_torque', inf);  
end 

function def_filters = default_filters()
    def_filters.omega_max = 0;      % max speed 
    def_filters.omega_rms = 0; 
    def_filters.tau_rms = 0;    % rms output torque 
    def_filters.tau_max = 0;    % max output torque 
    def_filters.total_mass = inf; 
    def_filters.effective_inertia = inf; 
    def_filters.manufacturer = {}; % Empty means ANY is ok 
end 

function [motor_keys, motor_values] = add_motors(motor_files)
    for i = 1:length(motor_files)
        [motor_keys_tmp, motor_cells_tmp] = read_motor_file(motor_files{i}); 
        if i == 1
            motor_keys = motor_keys_tmp;
            motor_values = motor_cells_tmp;
        else % append 
            % NOTE may want a faster append method 
            motor_keys = [motor_keys, motor_keys_tmp];
            motor_values = [motor_values; motor_cells_tmp]; 
        end 
    end 
end 

function [gb_keys, gb_values] = add_gearboxes(gearbox_files)
    for i = 1:length(gearbox_files)
        [gb_keys_tmp, gb_cells_tmp] = read_gearbox_file(gearbox_files{i}); 
        if i == 1
            gb_keys = gb_keys_tmp;
            gb_values = gb_cells_tmp;
        else % append 
            % NOTE may want a faster append method 
            gb_keys = [gb_keys, gb_keys_tmp];
            gb_values = [gb_values; gb_cells_tmp]; 
        end 
    end 
end 


function [motor_keys, motor_cells] = read_motor_file(motor_file)
    T = readtable(motor_file);    % NOTE: read table can be slow - may be able to speed up by specifing delims and other things
    motor_struct = table2struct(T);
    motor_keys = {motor_struct(:).key};
    motor_cells = num2cell(motor_struct); 
end 

function [gearbox_keys, gearbox_cells] = read_gearbox_file(gearbox_file)
    T = readtable(gearbox_file);
    gearbox_struct = table2struct(T);
    gearbox_keys = {gearbox_struct(:).key};
    gearbox_cells = num2cell(gearbox_struct); 
end 


function matches = find_matches(key_list, key)
    match_idxs = startsWith(key_list, key); 
    matches = key_list(match_idxs);
end 


%{
% TODO -- move to utility 
function idx = find_in_sorted()

    % TODO -- actually write this, for now, lazy :) 


end 
%} 