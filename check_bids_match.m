function [is_matched, file_list1, file_list2] = check_bids_match(file_list1, file_list2, do_sort)
% CHECK_BIDS_MATCH  Verify that two dir() file lists share identical
%   BIDS entities (sub, ses, run, task) entry-by-entry.
%
%   [is_matched, sorted1, sorted2] = check_bids_match(list1, list2)
%   [is_matched, sorted1, sorted2] = check_bids_match(list1, list2, do_sort)
%
%   Inputs
%       file_list1, file_list2 : struct arrays returned by dir()
%       do_sort (optional)     : if true (default), sort both lists by
%                                filename before comparing
%
%   Outputs
%       is_matched  : true if every pair shares the same sub/ses/run/task
%       file_list1  : (possibly sorted) first list
%       file_list2  : (possibly sorted) second list
%
%   Example
%       cifti  = dir(fullfile(sub_dir, '**', 'func', '*bold.dtseries.nii'));
%       conf   = dir(fullfile(sub_dir, '**', 'func', '*confounds*.tsv'));
%       [ok, cifti, conf] = check_bids_match(cifti, conf);
%       if ~ok, error('File mismatch'); end

    if nargin < 3, do_sort = true; end

    % Sort by filename
    if do_sort
        [~, idx1] = sort({file_list1.name});
        file_list1 = file_list1(idx1);
        if any(idx1 ~= 1:length(idx1))
            fprintf('The order of file_list1 is changed. \n');
        end
        [~, idx2] = sort({file_list2.name});
        file_list2 = file_list2(idx2);
        if any(idx2 ~= 1:length(idx2))
            fprintf('The order of file_list2 is changed. \n');
        end
    end

    % Length check
    if length(file_list1) ~= length(file_list2)
        warning('check_bids_match:lengthMismatch', ...
            'File counts differ: %d vs %d', length(file_list1), length(file_list2));
        is_matched = false;
        return;
    end

    is_matched = true;
    for i = 1:length(file_list1)
        e1 = parse_bids_entities(file_list1(i).name);
        e2 = parse_bids_entities(file_list2(i).name);

        fields_to_check = {'sub', 'ses', 'run', 'task'};
        for f = 1:length(fields_to_check)
            key = fields_to_check{f};
            v1 = e1.(key);
            v2 = e2.(key);

            if ~strcmp(v1, v2)
                warning('check_bids_match:entityMismatch', ...
                    'File %d: %s mismatch  ("%s" vs "%s")\n  %s\n  %s', ...
                    i, key, v1, v2, file_list1(i).name, file_list2(i).name);
                is_matched = false;
            end
        end
    end

    if is_matched
        fprintf('check_bids_match: All %d file pairs matched (sub/ses/run/task).\n', ...
            length(file_list1));
    end
end

function entities = parse_bids_entities(filename)
% Parse BIDS key-value entities from a filename string.
%   e.g.  'sub-01_ses-02_task-rest_run-03_bold.nii' ->
%         struct('sub','01','ses','02','task','rest','run','03')

    entities = struct('sub','', 'ses','', 'run','', 'task','');

    keys = fieldnames(entities);
    for k = 1:length(keys)
        pattern = [keys{k}, '-([a-zA-Z0-9]+)'];
        tok = regexp(filename, pattern, 'tokens', 'once');
        if ~isempty(tok)
            entities.(keys{k}) = tok{1};
        end
    end
end