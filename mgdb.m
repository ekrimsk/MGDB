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
%       num_combos - number of possible compatible combinations in db
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
    motors          % cells 
    motor_map       % key to index key to index in cell array 

    gearboxes       % cells 
    gearbox_map     % key to index in cell array 

    compatibility   % sparse 
    dd_key   
    num_combos

    gearbox_fields 
    motor_fields

    verbosity 
end 


methods (Access = public)


    function obj = mgdb(varargin) %
        init_tic = tic; 

        if isempty(varargin)
            % TODO -- robustify 
            database_paths = find_dirs('MGDB'); % look for a directory called mgdb
        else 
            database_paths = varargin(1);
        end 

        obj.dd_key = 'DD'; 
        obj.verbosity = 1; 

        obj.motor_fields = get_motor_fields();
        obj.gearbox_fields = get_gearbox_fields();

        %obj.motor_filter_signs = get_motor_filter_signs();
        %obj.gearbox_filter_signs = get_gearbox_filter_signs();


        % motor files called "*_motors.csv"
        % gearbox files called "*_gearboxes.csv"
        % compatibility filed called "*_compatibility.csv"
        % see github readme for formatting 
        for i = 1:length(database_paths)
            pth = database_paths{i};
            % TODO windows compatibility check 
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
        obj.add_motors(motor_files); 
        % Read in Gearbox csvs --NOTE: may move to its own function taking in list 
        obj.add_gearboxes(gearbox_files); 
        % add the 'dummy' gearbox for a direct drive option (compatible with all motors)
        

        % obj.gearboxes('DD') = direct_drive(); % add direct drive option to the map 

        %obj.compatibility = containers.Map('KeyType', 'char', 'ValueType', 'any'); % value type will be other maps 
        
        num_motors = numel(obj.motors);
        num_gearboxes = numel(obj.gearboxes); 

        % Read in compatibility csvs 
        add_compatibility(obj, compat_files); 
        % Now add direct drive option to all gearboxes 
        obj.compatibility(:, 1) = 1; % all motors compatible with direct drive 


        init_time = toc(init_tic); 
        fprintf('Added %d motors and %d gearboxes to database in %0.2f seconds\n',...
                   num_motors, num_gearboxes, init_time); 
    end 

    % TODO -- if verbosity is the only setting remove the struct 
    function update_settings(obj, new_settings)
        obj.settings.verbose = new_settings.verbose; 
    end 

    function filters = get_filters(obj)
        filters = obj.filters;
    end 
end 


methods (Access = private)

    


    


    function vprintf(obj, priority, varargin)
    % vprintf - verbose printf 
    %   
    %   if 'priority' >= obj.settings.verbose, 
    %   varargin{:} gets passed to printf
    %
    %    
        if priority >= obj.verbosity
            fprintf(varargin{:}); 
        end 
    end 


    function [motor_keys, motor_values] = add_motors(obj, motor_files)
    %
    %
    %
    %
    %
        for i = 1:length(motor_files)
            obj.vprintf(1, 'Reading motors from: %s\n', motor_files{i});
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

        keep = true(numel(motor_keys), 1);

        % Remove any motor this is not valid 
        for i = 1:numel(motor_keys)
            [valid, reason_str] = isvalid_motor(motor_values{i});
            if ~valid
                keep(i) = false; % dont keep it
                obj.vprintf(1, reason_str); 
            end 
        end 

        motor_keys = motor_keys(keep);
        motor_values = motor_values(keep); 

        if length(motor_keys) ~= length(unique(motor_keys))
            error('Multiple motors with same key specified');
        end 

    
        obj.motors = motor_values;  % cell array 
        obj.motor_map = containers.Map(motor_keys, 1:length(motor_values)); 
    end 

    function [gb_keys, gb_values] = add_gearboxes(obj, gearbox_files)

        gb_keys = {obj.dd_key};
        gb_values = {direct_drive()}; 
        % NOTE: direct drive added in instantioation 
        for i = 1:length(gearbox_files)
            obj.vprintf(1, 'Reading gearboxes from: %s\n', gearbox_files{i});
            [gb_keys_tmp, gb_cells_tmp] = read_gearbox_file(gearbox_files{i}); 
                            % NOTE may want a faster append method 
     
            gb_keys = [gb_keys, gb_keys_tmp];
            gb_values = [gb_values; gb_cells_tmp]; 
            
        end 

        if length(gb_keys) ~= length(unique(gb_keys))
            error('Multiple gearboxes with same key specified');
        end 

        % gearbox number 1 is always direct drive 
        obj.gearboxes = gb_values;  % cell array 
        obj.gearbox_map = containers.Map(gb_keys, 1:length(gb_values)); 
    end 

    function add_compatibility(obj, compat_files)
    %
    %
    %
    %
    %
        % could combine all compat files together to make simpler 
        % assume this is done somewhere above this line of code 
        num_motors = numel(obj.motors);
        num_gearboxes = numel(obj.gearboxes); 
        compatibility = sparse([], [], [], num_motors, num_gearboxes,...
                                 min(100*num_motors, num_motors*num_gearboxes)); 

        gearbox_keys = obj.gearbox_map.keys; % WILL RETURN THE KEYS SORTED ALPHABETICALLY 

        % matlab containers isKey is very slow when the key is not a key 

        % TODO -- check that compat files is cell array even if only 1 file 
        num_combos = 0; 
        for kk = 1:length(compat_files)
            compat_file = compat_files{kk}; 
            compat_fid = fopen(compat_file);

            tline = fgetl(compat_fid);
            line_count = 1; % for debug 
            while ischar(tline)     % go through whole file
                line_cells = strsplit(tline, ',');  % might be no compat 
                motor_key = line_cells{1};

                try 
                    motor_idx = obj.motor_map(motor_key);   % the index   
                catch ME 
                    error('Motor in compatibility file not in database')
                end 

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

                    if endsWith(gb_key, '*')        % partial key -- any key that ends with star is a partial key
                        new_gb_match_keys = find_matches(gearbox_keys, gb_key(1:end-1)); % TODO -- can rewrite a faster subroutine 
                    % there needs to be a faster way 
                    elseif obj.gearbox_map.isKey(gb_key) 
                        new_gb_match_keys = {gb_key};  % assumed cell array 
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

                % get the indices for all the match keys 
                for i = 1:length(gb_match_keys)
                    gb_idx = obj.gearbox_map(gb_match_keys{i});
                    compatibility(motor_idx, gb_idx) = 1; 
                end 

                %{
                if obj.compatibility.isKey(motor_key)
                    % already have things added, append 
                    for j = 1:length(gb_match_keys)
                        gb_key = gb_match_keys{j};
                        tmp_map = obj.compatibility(motor_key);     % should be copy by ref
                        if tmp_map.isKey(gb_key)
                            warning('redundant compatibilities specified')
                        else 
                            num_combos = num_combos + 1; 
                            tmp_map(gb_key) = inf; % initialize all combination to have infinite cost 
                        end 
                    end 
                else 
                    % Add all the new ones together AND a direct drive option
                    % It is much faster to add lots of things to the map
                    % at once so it does not need to get resized and copied
                    % NOTE: some might already be present though  
                    inf_cell = num2cell(inf(1, length(gb_match_keys) + 1)); 
                    obj.compatibility(motor_key) = containers.Map([obj.dd_key, gb_match_keys], inf_cell); 
                    num_combos = num_combos + length(gb_match_keys) + 1; % + 1 for the dd as well 
                end 
                %}


                tline = fgetl(compat_fid);
                line_count = line_count + 1; 
            end % finish reading through lines of this file 
            fclose(compat_fid); 
        end % end looping through files 

        % Now add direct drive for any motors that were skipped 
        %{
        motor_keys = obj.motors.keys;
        for m_idx = 1:length(motor_keys)
            motor_key = motor_keys{m_idx};
            if ~obj.compatibility.isKey(motor_key)
                obj.compatibility(motor_key) = containers.Map(obj.dd_key, inf); 
                num_combos = num_combos + 1; 
            end 
        end 
        %} 
        num_combos = nnz(compatibility); 
        obj.num_combos = num_combos; 
        obj.compatibility = compatibility; 
    end 





end % end private methods 

%=================== END CLASS DEF =========================
end 


%   NOTE: maybe move specification of fields and limits to another file
%   for better software design (kind of nice to be self contained though)
function fields = get_motor_fields()
    fields = {'key', 'manufacturer', 'ID', 'type', 'V', 'k_t', 'R',...
                     'L', 'mass', 'inertia', 'omega_nl', 'I_nl', 'I_nom'...
                    'max_int_torque', 'max_int_speed', 'max_cont_speed',...
                    'max_cont_power','coulomb_friction','viscous_friction',...
                    'Rth1', 'Rth2'};
end 

function fields = get_gearbox_fields()
    fields = {'key', 'manufacturer', 'ID', 'type', 'stages', 'ratio',...
                     'mass', 'inertia', 'efficiency', 'direction',...
                     'max_int_torque', 'max_cont_torque'};

end 





function [valid, str] = isvalid_motor(motor)
    valid = true; % unless prove otherwise
    str = [];     % default for valid 

    % validate individaul fields 
    % each field has a required type 
    % some fields will have bounds (e.g. nonnegativity constraints)
    % limits on types 
    % motor class and gearbox class to hold some values not a bad idea 

    pos_vals = [motor.V, motor.R, motor.L, motor.mass, motor.inertia, motor.I_nom];

    if ~valid_motor_friction(motor)
        valid = false;
        reason = ['invalid friction parameters. '...
                        'Implied no-load current may be negative'];
                        
    elseif any(isnan(pos_vals)) || any(pos_vals < 0)
        valid = false; 
        reason = ['V, R, L, mass, inertia, I_nom, must all be non-negative reals'];
    end 

    if ~valid 
        str = sprintf('Motor %s skipped: %s\n', motor.key, reason);
    end 
end 

function valid = valid_motor_friction(motor)

    % check if no-load current is nan 
    valid = true; 
    if isnan(motor.I_nl)
        % maybe coulomb or viscous is spefied 
        if ~isnan(motor.coulomb_friction) || ~isnan(motor.viscous_friction)
            if (motor.coulomb_friction < 0) || (motor.viscous_friction < 0)
                valid = false;
            end 
        else
            % do nothing,
            I_nl = (motor.V - motor.k_t*motor.omega_nl)/motor.R; 
            if I_nl < 0 
                valid = false; 
            end 
        end 
    end 
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




function [motor_keys, motor_cells] = read_motor_file(motor_file)
    % NOTE: read table can be slow - may be able to speed up by specifing delims and other things
    T = readtable(motor_file);   
    motor_struct = table2struct(T);
    motor_keys = {motor_struct(:).key};
    motor_cells = num2cell(motor_struct); 
    % validate all the field names types 

    col_names = T.Properties.VariableNames;

    expected_names = get_motor_fields(); 

    if numel(col_names) ~= numel(expected_names)
        error('Motor file %s missing columns', motor_file)
    end 

    for i = 1:length(col_names)
        if ~strcmp(col_names{i}, expected_names{i})
            error('Invalid motor column name %s', col_names{i});
        end 
    end 


end 


function [gearbox_keys, gearbox_cells] = read_gearbox_file(gearbox_file)
    T = readtable(gearbox_file);
    gearbox_struct = table2struct(T);
    gearbox_keys = {gearbox_struct(:).key};
    gearbox_cells = num2cell(gearbox_struct); 
    % validate all the field names and data types 

    col_names = T.Properties.VariableNames;

    expected_names = get_gearbox_fields(); 

    if numel(col_names) ~= numel(expected_names)
        error('Gearbox file %s missing columns', motor_file)
    end 

    for i = 1:length(col_names)
        if ~strcmp(col_names{i}, expected_names{i})
            error('Invalid gearox column name %s', col_names{i});
        end 
    end 


end 


function matches = find_matches(key_list, key)
    % TODO -- because this is sorted we can redo this in Olog(N) with 
    % binary search 
    match_idxs = startsWith(key_list, key); 
    matches = key_list(match_idxs);
end 

function val = cells_contain(key_list, key)
    val = false; 
    % key list is lexicographically sorted key list 
    % key is the FULL key to search for 

    % hardcoding for key_list of length 1 
    if length(key_list) == 1
        if strcmp(key_list{1}, key);
            val = true;
            return;
        end 
    end 

    lb_idx = 1;
    ub_idx = length(key_list);
    cur_idx = round((lb_idx + ub_idx)/2);

    while (ub_idx - lb_idx) >= 1
        cur_idx = round((lb_idx + ub_idx)/2);

        cur_str = key_list{cur_idx}; 

        cmp_tmp = cstrcmp(cur_str, key);
        if cmp_tmp > 0
            ub_idx = cur_idx;
        elseif cmp_tmp < 0 
            lb_idx = cur_idx; 
        else % strings equal 
            ub_idx = cur_idx;  % set same to cause break 
            lb_idx = cur_idx
        end 
    end 

    if strcmp(key_list{lb_idx}, key)
        val = true; 
    end 
end 

%https://www.mathworks.com/matlabcentral/answers/39374-compare-two-strings-based-on-ascii-dictionary-order
function cmp = cstrcmp( a, b )
    % Force the strings to equal length
    x = char({a;b});
    % Subtract one from the other
    d = x(1,:) - x(2,:);
    % Remove zero entries
    d(~d) = [];
    if isempty(d)
        cmp = 0;
    else
        cmp = d(1);
    end
end


