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
    motors          % containers.Map
    gearboxes       % containers.Map
    compatibility 
    settings 
    filters 
    dd_key   
    num_combos

    motor_filter_signs
    gearbox_filter_signs
    motor_fields
    gearbox_fields 
end 


methods (Access = public)


    function obj = mgdb(varargin) %
        % TODO -- add onputs to change database path 

        init_tic = tic; 

        if isempty(varargin)
            % TODO -- robustify 
            database_paths = find_dirs('MGDB'); % look for a directory called mgdb
        else 
            database_paths = varargin(1);
        end 

        obj.dd_key = 'DD'; 
        settings.verbose = 1;   % this is the only settings write now 
        obj.settings = settings; 

        obj.motor_fields = get_motor_fields();
        obj.gearbox_fields = get_gearbox_fields();
        obj.motor_filter_signs = get_motor_filter_signs();
        obj.gearbox_filter_signs = get_gearbox_filter_signs();

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
        obj.add_motors(motor_files); 
        % Read in Gearbox csvs --NOTE: may move to its own function taking in list 
        obj.add_gearboxes(gearbox_files); 

        % NOTE: A motor might appear mulitpe times in multiple compatibility 
        % files (could be compatible with multiple manufactureres for example)
        % TODO -- move into add 
  
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
        obj.filters = obj.default_filters();
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

        obj.vprintf(1, 'Found %d initial combinations\n', obj.num_combos);
        obj.filter_stats(); 

        num_init_combos = 0; % sanittyc checeeasgs

        % For each motor, get all compatible gearboxes 
        for m_idx = 1:length(motor_keys)
            motor_key = motor_keys{m_idx}; 
            tmp_map = obj.compatibility(motor_key);
            gb_keys = tmp_map.keys; 

            num_init_combos = num_init_combos + length(gb_keys);

            for gb_idx = 1:length(gb_keys) 
                gb_key = gb_keys{gb_idx};
                % Check if combo is valid given the filters
                if obj.valid_combo(motor_key, gb_key)
                    motor_combo_keys{end + 1} = motor_key;
                    gearbox_combo_keys{end + 1} = gb_key; 
                    cost_list{end + 1} = tmp_map(gb_key);
                end               
            end 
        end 

        num_valid_combos = length(motor_combo_keys);
        num_removed = obj.num_combos - num_valid_combos;
        obj.vprintf(1, 'Filtering removed %d combinations\n', num_removed);


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

    function clear_filters(obj, varargin) 
    %
    %
    %
        def_filters = default_filters(); 
        if isempty(varargin)  % clear all 
            obj.filters = def_filters(); % overwrite with defaults 
        else 
            % TODO -- make robust (check valid criteria)
            criteria = varargin{1};
            obj.filter.(criteria) = def_filters.(criteria);
        end 
    end 


    function update_filters(obj, new_filters)  % TODO switch to struct base 
    % update_filters
    %   new_filters - struct with (any) of the following fields:
    %       'omega_max' - lower cutoff for gearbox max_int_speed (rad/s) or 
    %                   motor speed (when mulitplied by gear ratio)
    %       'max_cont_speed' - lower cutoff for gearbox max_cont_speed (rad/s)
    %                   or motor speed when (when multiplied by gear ratio)
    %       'tau_max' -  lower bound cutoff for gearbox max_int
    %       'total_mass' - upper bound cutoff of gearbox + motor mass (kg)
    %       'effective_inerita' - upper bound cutoff for combined effective
    %                           inertia in kgm^2
    %       'manufacturer' - cell array of manufacturer names. If cell array
    %                is empty, than ALL motor/gearbox manufacturers are accepted
    %       'motors' - for filtering on motor specific properties
    %                       TODO more info 
    %
    %       'gearboxes' - for filtering on geabox specific properties 
    %                       % TODO -- more info 

        fields = fieldnames(new_filters); 
        for ii = 1:numel(fields)
            criteria = fields{ii};
            value = new_filters.(criteria);
            if strcmp(criteria, 'motors')
                obj.update_motor_filters(values)
            elseif strcmp(criteria, 'gearboxes')
                obj.update_gearbox_filters(values)
            else 
                % TODO -- add char check 
                if ~iscell(value)
                    assert(value >= 0, 'update_filter: values must be nonnegative'); 
                end 
                if isfield(obj.filters, criteria)
                    obj.filters.(criteria) = value;
                else 
                    warning('update_filter: invalid filter criteria: %s', criteria);
                end 
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

    function filter_stats(obj)
        fn = fieldnames(obj.filters); 

        for ii = 1:length(fn)
            criteria = fn{ii}; 
            val = obj.filters.(criteria);
            if ~isempty(val)
                if iscell(val)
                    for kk = 1:length(val)
                        obj.vprintf(1, 'Filtering combination with %s: %s\n',...
                        criteria, val{kk});
                    end 
                else 
                    if all(~isinf(val)) && all(val~=0) 
                        obj.vprintf(1, 'Filtering combination with %s: %0.2f\n',...
                                        criteria, val);
                    end 
                end 
            end 
        end 
    end 


    function update_motor_filters(obj, new_filters)
        motor_filter = obj.filters.motors;
        signs = obj.motor_filter_signs;
        fields = fieldnames(new_filters); 
        for ii = 1:numel(fields)
            field = fields{ii}; 
            if signs ~= '~'   % i.e. this field is setable as a filter 
                error('Cannot set filter for motor field %s', field);
            else 
                motor_filter.(field) = new_filters.(field);
            end 
        end 
        obj.filters.motors = motor_filter; 
    end 


    function update_gearbox_filters(obj, new_filters)
        gearbox_filter = obj.filters.motors;
        signs = obj.gearbox_filter_signs;
        fields = fieldnames(new_filters); 
        for ii = 1:numel(fields)
            field = fields{ii}; 
            if signs{ii} == '~'   % i.e. this field is setable as a filter 
                error('Cannot set filter for gearbox field %s', field); 
            else 
                gearbox_filter.(field) = gearbox_filters.(field);
            end 
        end 
        obj.filters.gearboxes = gearbox_filter; 
    end 


    function [valid] = valid_combo(obj, motor_key, gearbox_key)
    %
    %
    %
    %
    %
        % Check this combination of motor and gearbox against the filters 
        valid = true;

        motor = obj.motors(motor_key);
        gearbox = obj.gearboxes(gearbox_key);

        valid_manufs = obj.filters.manufacturer;

        % filter motor and gearbox individually 
        if ~obj.filter_motor(motor)
            valid = false; 
        elseif ~obj.filter_gearbox(gearbox)
            valid = false;
        end 


        %% Filtering on combination 

        % Maximum speed 
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
                validatestring(lower(motor.manufacturer),lower(valid_manufs));
                validatestring(lower(gearbox.manufacturer),lower(valid_manufs));
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
                    %if obj.gearboxes.isKey(gb_key)  % slow 
                    if endsWith(gb_key, '*')        % partial key -- any key that ends with star is a partial key
                        new_gb_match_keys = find_matches(gearbox_keys, gb_key(1:end-1)); % TODO -- can rewrite a faster subroutine 
                    elseif obj.gearboxes.isKey(gb_key)  % there needs to be a faster way 
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


                tline = fgetl(compat_fid);
                line_count = line_count + 1; 
            end % finish reading through lines of this file 
            fclose(compat_fid); 
        end % end looping through files 

        % Now add direct drive for any motors that were skipped 
        motor_keys = obj.motors.keys;
        for m_idx = 1:length(motor_keys)
            motor_key = motor_keys{m_idx};
            if ~obj.compatibility.isKey(motor_key)
                obj.compatibility(motor_key) = containers.Map(obj.dd_key, inf); 
                num_combos = num_combos + 1; 
            end 
        end 

        obj.num_combos = num_combos; 
    end 



    function valid = filter_motor(motor)
        valid = true; 
        filters = obj.filters.motor_filters; 

        fields = obj.motor_fields; 
        signs = obj.motor_filter_signs; 

        for i = 1:length(fields)
            field = fields{i};
            sgn = signs.(field);
            motor_val = motor.(field); 
            compare_val = filters.(field); 

            if sgn == '~'
                % skip
            elseif sgn == '<'
                if motor_val > compare_val 
                    valid = false;
                    break;
                end 
            elseif sgn == '>'
                if motor_val < compare_val 
                    valid = false;
                    break;
                end 
            elseif sgn == '='
                if ischar(compare_val)
                    compare_val = {compare_val}; % should be cell array of strings
                end 

                try 
                    validatestring(motor_val, compare_val); 
                catch ME
                    if strcmp(ME.identifier,'MATLAB:unrecognizedStringChoice')
                        valid = false;
                        break;
                    else 
                        rethrow(ME);
                    end 
                end 
            else 
                error('invalid sign');
            end 
        end 

    end 


    function valid = filter_gearbox(gearbox)
        valid = true; 
        filters = obj.filters.gearbox_filters; 

        fields = obj.gearbox_fields; 
        signs = obj.gearbox_filter_signs; 

        for i = 1:length(fields)
            field = fields{i};
            sgn = signs.(field);
            val = gearbox.(field); 
            compare_val = filters.(field); 

            if sgn == '~'
                % skip
            elseif sgn == '<'
                if val > compare_val 
                    valid = false;
                    break;
                end 
            elseif sgn == '>'
                if val < compare_val 
                    valid = false;
                    break;
                end 
            elseif sgn == '='
                if ischar(compare_val)
                    compare_val = {compare_val}; % should be cell array of strings
                end 

                try 
                    validatestring(val, compare_val); 
                catch ME
                    if strcmp(ME.identifier,'MATLAB:unrecognizedStringChoice')
                        valid = false;
                        break;
                    else 
                        rethrow(ME);
                    end 
                end 
            else 
                error('invalid sign');
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

    function [motor_keys, motor_values] = add_motors(obj, motor_files)
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

        obj.motors = containers.Map(motor_keys, motor_values);
    end 

    function [gb_keys, gb_values] = add_gearboxes(obj, gearbox_files)
        for i = 1:length(gearbox_files)
            obj.vprintf(1, 'Reading gearboxes from: %s\n', gearbox_files{i});
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

        if length(gb_keys) ~= length(unique(gb_keys))
            error('Multiple gearboxes with same key specified');
        end 

        obj.gearboxes = containers.Map(gb_keys, gb_values);
    end 


    function def_filters = default_filters(obj)
        def_filters.omega_max = 0;      % max speed 
        def_filters.omega_rms = 0; 
        def_filters.tau_rms = 0;    % rms output torque 
        def_filters.tau_max = 0;    % max output torque 
        def_filters.total_mass = inf; 
        def_filters.effective_inertia = inf; 
        def_filters.manufacturer = {}; % Empty means ANY is ok 

        def_filters.motors = obj.default_motor_filter(); 
        def_filters.gearboxes = obj.default_gearbox_filter(); 
    end 

    function def_filter = default_motor_filter(obj)
        fields = obj.motor_fields; 
        signs = obj.motor_filter_signs; 
        for i = 1:length(fields)
            field = fields{i};
            sgn = signs.(field);
            if sgn == '~'
                % skip
            elseif sgn == '<'
                def_filter.(field) = inf;  
            elseif sgn == '>'
                def_filter.(field) = -inf;  
            elseif sgn == '='
                def_filter.(field) = {}; % anything 
            else 
                error('invalid sign');
            end 
        end 
    end 


    function def_filter = default_gearbox_filter(obj)
        fields = obj.gearbox_fields; 
        signs = obj.gearbox_filter_signs; 
        for i = 1:length(fields)
            field = fields{i};
            sgn = signs.(field);
            if sgn == '~'
                % skip
            elseif sgn == '<'
                def_filter.(field) = inf;  
            elseif sgn == '>'
                def_filter.(field) = -inf;  
            elseif sgn == '='
                def_filter.(field) = {}; % anything 
            else 
                error('invalid sign');
            end 
        end 
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

function signs = get_motor_filter_signs()
%   TODO -- this needs lots of documenting 
%
%   '='  - direct match (only on strings)
%   '<'  - motor parameter must be less than this 
%   '>'  - motor parameter must be greater than this 
%   '~'  - filtering on this parameter NOT allowed 
    signs.key = '~';
    signs.manufacturer = '=';
    signs.ID = '~';
    signs.type = '=';
    signs.V = '>';
    signs.k_t = '>';
    signs.R = '<';
    signs.L = '~';
    signs.mass = '<';
    signs.inertia = '<';
    signs.omega_nl = '>';
    signs.I_nl = '~';
    signs.I_nom = '~'; 
    signs.max_int_torque = '>';
    signs.max_int_speed = '>';
    signs.max_cont_speed = '>';
    signs.max_cont_power = '>';
    signs.coulomb_friction = '~';
    signs.viscous_friction = '~';
    signs.Rth1 = '~';
    signs.Rth2 = '~';
end 

function signs = get_gearbox_filter_signs()
%   TODO -- this needs lots of documenting 
%
%   '='  - direct match (only on strings)
%   '<'  - motor parameter must be less than this 
%   '>'  - motor parameter must be greater than this 
%   '~'  - filtering on this parameter NOT allowed 
    signs.key = '~';
    signs.manufacturer = '=';
    signs.ID = '~';
    signs.type = '=';
    signs.stages = '<';
    signs.ratio = '<';
    signs.mass = '<';
    signs.inertia = '<';
    signs.efficiency = '>';
    signs.direction = '='; 
    signs.max_int_torque = '>';
    signs.max_cont_torque = '>';
end


% TODO --- funtion to update the individual motor + gearbox filters 




function [valid, str] = isvalid_motor(motor)
    valid = true; % unless prove otherwise
    str = [];     % default for valid 

    % validate individaul fields 
    % each field has a required type 
    % some fields will have bounds (e.g. nonnegativity constraints)
    % limits on types 
    % motor class and gearbox class to hold some values not a bad idea 

    pos_vals = [motor.V, motor.R, motor.L, motor.mass, motor.inertia];


    if ~valid_motor_friction(motor)
        valid = false;
        reason = ['invalid friction parameters. '...
                        'Implied no-load current may be negative'];
    elseif any(isnan(pos_vals)) || any(pos_vals < 0)
        valid = false; 
        reason = ['V, R, L, mass, inertia must all be non-negative reals'];
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


%{
% TODO -- move to utility 
function idx = find_in_sorted()

    % TODO -- actually write this, for now, lazy :) 


end 
%} 