function R = PACKAGEVerifierNEWTON_fromContents(varargin)
%PACKAGEVerifierNEWTON_fromContents  Verify leigqNEWTON package integrity (public/private-safe).
%
% Why this exists:
%   - Helps you verify that the toolbox folder is self-consistent (no missing pieces).
%   - Detects path shadowing (multiple versions on the MATLAB path).
%   - Warns if root-level .m files are not listed in Contents.m.
%
% The verifier is designed to work for both public and private distributions:
%   - Metadata files (README/LICENSE/CITATION) may differ or be absent -> configurable.
%   - Docs folders may be doc/source or docs/source -> detected automatically.
%
% Usage:
%   R = PACKAGEVerifierNEWTON_fromContents
%   R = PACKAGEVerifierNEWTON_fromContents('MetaMode','auto')   % default
%   R = PACKAGEVerifierNEWTON_fromContents('MetaMode','strict') % public-release QA
%
% Options (name-value):
%   'RootDir'        : toolbox root. Default: folder of this file.
%   'Verbose'        : true/false (default true)
%   'CheckUnlisted'  : true/false (default true)
%   'ScanSubfolders' : true/false (default false)
%   'CheckDocs'      : true/false (default true)
%   'CheckCalls'     : true/false (default true)
%   'CheckMeta'      : true/false (default true)
%   'MetaMode'       : 'auto'|'warn'|'strict'|'off' (default 'auto')
%
% Notes on Contents.m parsing:
%   This verifier accepts both explicit entries
%     %   leigqNEWTON_refine_auto - ...
%   and wildcard entries
%     %   leigqNEWTON_refine_*    - ...
%   where '*' is interpreted as "any function starting with this prefix".
%
% See also: checkNEWTON, leigqNEWTON

t0 = tic;

p = inputParser;
p.FunctionName = mfilename;
addParameter(p,'RootDir','',@(s)ischar(s)||isstring(s));
addParameter(p,'Verbose',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'CheckUnlisted',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'ScanSubfolders',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'CheckDocs',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'CheckCalls',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'CheckMeta',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'MetaMode','auto',@(s)ischar(s)||isstring(s));
parse(p,varargin{:});
opt = p.Results;

Verbose        = logical(opt.Verbose);
CheckUnlisted  = logical(opt.CheckUnlisted);
ScanSubfolders = logical(opt.ScanSubfolders);
CheckDocs      = logical(opt.CheckDocs);
CheckCalls     = logical(opt.CheckCalls);
CheckMeta      = logical(opt.CheckMeta);
MetaMode       = lower(string(opt.MetaMode));

R = struct();
R.ok = true;
R.when = datestr(now);
R.rootDir = '';
R.contentsFile = '';
R.entries = struct();
R.unlisted = struct();
R.docs = struct();
R.calls = struct();
R.meta = struct();
R.errors = {};
R.warnings = {};
R.elapsed = NaN;

% Root
if strlength(string(opt.RootDir))==0
    rootDir = fileparts(mfilename('fullpath'));
else
    rootDir = char(opt.RootDir);
end
R.rootDir = rootDir;

% Contents
contentsFile = fullfile(rootDir,'Contents.m');
R.contentsFile = contentsFile;
if exist(contentsFile,'file')~=2
    R = local_err(R, sprintf('Missing Contents.m at: %s', contentsFile), Verbose);
    R.elapsed = toc(t0);
    return;
end

% Root .m files
[rootFns, scanned] = local_root_mfiles(rootDir, ScanSubfolders);
R.unlisted.scanned = scanned;

% Parse Contents entries (explicit + wildcard prefixes)
[explicit, prefixes] = local_parse_contents_entries(contentsFile);
R.entries.explicit = explicit;
R.entries.prefixes = prefixes;

if isempty(explicit) && isempty(prefixes)
    R = local_err(R, 'No entries found in Contents.m.', Verbose);
    R.elapsed = toc(t0);
    return;
end

% Expand wildcard prefixes using rootFns
listed = string(explicit);
for i=1:numel(prefixes)
    pref = string(prefixes{i});
    hits = string(rootFns(startsWith(string(rootFns), pref)));
    listed = [listed; hits(:)]; %#ok<AGROW>
end
listed = unique(listed,'stable');

% Keep only those that exist in root (avoid pulling in example/docs names)
listed = listed(ismember(listed, string(rootFns)));
fnList = cellstr(listed);

R.entries.root_functions = rootFns;
R.entries.used_for_checks = fnList;

if Verbose
    fprintf('=== NEWTON PACKAGE verifier from Contents ===\n');
    fprintf('Root    : %s\n', rootDir);
    fprintf('Contents: %s\n', contentsFile);
    fprintf('Found %d explicit entries and %d wildcard prefixes in Contents.\n', numel(explicit), numel(prefixes));
    fprintf('Root functions: %d, Covered by Contents: %d\n\n', numel(rootFns), numel(fnList));
end

% 1) Unlisted .m files in root
if CheckUnlisted
    unlisted = setdiff(string(rootFns), string(fnList), 'stable');
    unlisted = cellstr(unlisted);
    R.unlisted.names = unlisted;

    if Verbose && isfield(scanned,'excluded') && ~isempty(scanned.excluded)
        fprintf('Excluded (intentionally ignored for “unlisted” check):\n');
        for k=1:numel(scanned.excluded)
            fprintf('  - %s\n', scanned.excluded{k});
        end
        fprintf('\n');
    end

    if ~isempty(unlisted)
        R = local_warn(R, sprintf('Found %d .m files not listed in Contents.', numel(unlisted)), Verbose);
        if Verbose
            fprintf('Unlisted .m files:\n');
            for k=1:numel(unlisted)
                fprintf('  - %s\n', unlisted{k});
            end
            fprintf('\n');
        end
    end
end

% 2) Docs (doc/source or docs/source; accept 00_index.m or 00_index_NEWTON.m)
if CheckDocs
    docDir = local_find_doc_source_dir(rootDir);
    if ~isempty(docDir)
        R.docs.sourceDir = docDir;

        idxCandidates = { ...
            fullfile(docDir,'00_index.m'), fullfile(docDir,'00_index_NEWTON.m'), ...
            fullfile(docDir,'00_index.mlx'), fullfile(docDir,'00_index_NEWTON.mlx')};
        haveIdx = false;
        for k=1:numel(idxCandidates)
            if exist(idxCandidates{k},'file')==2
                haveIdx = true;
                R.docs.indexFile = idxCandidates{k};
                break;
            end
        end
        if ~haveIdx
            R = local_warn(R, 'Doc index missing (00_index.m or 00_index_NEWTON.m).', Verbose);
        end

        % For doc pages: check only functions covered by Contents (fnList)
        for i=1:numel(fnList)
            fname = fnList{i};
            cand = { ...
                fullfile(docDir, ['doc_', fname, '.m']), fullfile(docDir, ['doc_', fname, '.mlx']), ...
                fullfile(docDir, [fname, '_doc.m']), fullfile(docDir, [fname, '_doc.mlx'])};
            ok = false;
            for k=1:numel(cand)
                if exist(cand{k},'file')==2
                    ok = true; break;
                end
            end
            if ~ok
                R = local_warn(R, sprintf('Doc page missing for %s (expected doc_%s.m or equivalent).', fname, fname), Verbose);
            end
        end
    else
        R = local_warn(R, 'Docs source folder not found (skipping doc check).', Verbose);
    end
end

% 3) Which shadowing checks
if CheckCalls
    for i=1:numel(fnList)
        f = fnList{i};
        paths = local_which_all(f);
        R.calls.(f) = paths;
        if numel(paths)>1
            R = local_warn(R, sprintf('Shadowing detected for %s (multiple versions on path).', f), Verbose);
        end
    end
end

% 4) Metadata checks (README/LICENSE/CITATION), safe for public/private
if CheckMeta && MetaMode~="off"
    meta = local_check_meta(rootDir);
    meta.mode = char(MetaMode);

    doCheck = true;
    if MetaMode=="auto"
        doCheck = meta.any_present;
    end

    if doCheck
        R = local_apply_meta(R, meta, MetaMode, Verbose);
    end
    R.meta = meta;
end

R.elapsed = toc(t0);

if Verbose
    fprintf('=== Result ===\n');
    if R.ok
        fprintf('OK. Elapsed %.2fs\n', R.elapsed);
    else
        fprintf('FAILED. Elapsed %.2fs\n', R.elapsed);
    end
end
end

% -------------------------------------------------------------------------
function [explicit, prefixes] = local_parse_contents_entries(contentsFile)
txt = fileread(contentsFile);
lines = regexp(txt, '\r\n|\n|\r', 'split');

explicit = {};
prefixes = {};
for i=1:numel(lines)
    L0 = strtrim(lines{i});
    if ~startsWith(L0,'%'); continue; end

    % Remove leading '%'
    L = regexprep(L0,'^%','');
    % Accept list entries: indent >=2 spaces, then a token like name or prefix*
    tok = regexp(L, '^\s{2,}([A-Za-z]\w*\*?)\s+-\s*', 'tokens', 'once');
    if isempty(tok); continue; end

    name = tok{1};
    if endsWith(name,'*')
        prefixes{end+1,1} = extractBefore(name, strlength(name)); %#ok<AGROW>
    else
        explicit{end+1,1} = name; %#ok<AGROW>
    end
end

explicit = unique(explicit,'stable');
prefixes = unique(prefixes,'stable');
end

function [rootFns, scanned] = local_root_mfiles(rootDir, scanSubfolders)
if scanSubfolders
    files = dir(fullfile(rootDir,'**','*.m'));
else
    files = dir(fullfile(rootDir,'*.m'));
end

names = strings(0,1);
excluded = {};

rootPrefix = string(rootDir) + filesep;

for i=1:numel(files)
    fpath = fullfile(files(i).folder, files(i).name);
    rel = erase(string(fpath), rootPrefix);

    % Exclude docs/legacy trees
    if contains(rel, filesep + "docs" + filesep) || startsWith(rel, "docs"+filesep)
        excluded{end+1} = char(rel); %#ok<AGROW>
        continue;
    end
    if contains(rel, filesep + "doc" + filesep) || startsWith(rel, "doc"+filesep)
        excluded{end+1} = char(rel); %#ok<AGROW>
        continue;
    end
    if contains(rel, filesep + "legacy" + filesep) || startsWith(rel, "legacy"+filesep)
        excluded{end+1} = char(rel); %#ok<AGROW>
        continue;
    end

    base = string(files(i).name);
    baseLower = lower(base);

    % Exclude Contents and verifiers (these are intentionally not required to be listed)
    if base=="Contents.m"
        excluded{end+1} = char(base); %#ok<AGROW>
        continue;
    end
    if startsWith(baseLower, "packageverifiernewton_fromcontents")
        excluded{end+1} = char(base); %#ok<AGROW>
        continue;
    end

    fn = erase(base, ".m");
    names(end+1,1) = fn; %#ok<AGROW>
end

rootFns = cellstr(unique(names,'stable'));

scanned = struct();
scanned.count_scanned = numel(files);
scanned.count_kept = numel(rootFns);
scanned.excluded = excluded;
end

function docDir = local_find_doc_source_dir(rootDir)
candDirs = {fullfile(rootDir,'docs','source'), fullfile(rootDir,'doc','source')};
docDir = '';
for i=1:numel(candDirs)
    if exist(candDirs{i},'dir')==7
        docDir = candDirs{i}; return;
    end
end
end

function paths = local_which_all(fname)
w = which(fname,'-all');
if isempty(w)
    paths = {};
    return;
end
s = string(w);
s = splitlines(s);
s = strtrim(s);
s = s(s~="");
s = s(~startsWith(s,"built-in"));
paths = cellstr(s);
end

function meta = local_check_meta(rootDir)
meta = struct();
meta.rootDir = rootDir;

meta.readme.candidates  = {'README.md','README.txt','README_NEWTON_Examples.txt','README'};
meta.license.candidates = {'LICENSE','LICENSE.txt','COPYING','COPYING.txt'};
meta.citation.candidates= {'CITATION.cff','CITATION.bib','CITATION'};

meta.readme.present  = local_any_exist(rootDir, meta.readme.candidates);
meta.license.present = local_any_exist(rootDir, meta.license.candidates);
meta.citation.present= local_any_exist(rootDir, meta.citation.candidates);

meta.any_present = meta.readme.present || meta.license.present || meta.citation.present;
end

function tf = local_any_exist(rootDir, candidates)
tf = false;
for i=1:numel(candidates)
    if exist(fullfile(rootDir,candidates{i}),'file')==2
        tf = true; return;
    end
end
end

function R = local_apply_meta(R, meta, MetaMode, Verbose)
% README
if ~meta.readme.present
    msg = sprintf('Recommended file missing: %s', strjoin(meta.readme.candidates,', '));
    R = local_meta_msg(R, msg, MetaMode, Verbose);
end
% LICENSE
if ~meta.license.present
    msg = sprintf('Recommended file missing: %s', strjoin(meta.license.candidates,', '));
    R = local_meta_msg(R, msg, MetaMode, Verbose);
end
% CITATION
if ~meta.citation.present
    msg = sprintf('Recommended file(s) missing: %s', strjoin(meta.citation.candidates,', '));
    R = local_meta_msg(R, msg, MetaMode, Verbose);
end
end

function R = local_meta_msg(R, msg, MetaMode, Verbose)
if MetaMode=="strict"
    R = local_err(R, msg, Verbose);
else
    R = local_warn(R, msg, Verbose);
end
end

function R = local_warn(R, msg, Verbose)
R.warnings{end+1} = msg;
if Verbose
    fprintf('WARN : %s\n', msg);
end
end

function R = local_err(R, msg, Verbose)
R.errors{end+1} = msg;
R.ok = false;
if Verbose
    fprintf('ERROR: %s\n', msg);
end
end
