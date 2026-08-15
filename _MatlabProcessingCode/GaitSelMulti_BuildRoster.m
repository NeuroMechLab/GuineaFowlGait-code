function R = GaitSelMulti_BuildRoster(opts)
%GAITSELMULTI_BUILDROSTER  One roster of every guinea fowl LEVEL trial with raw
%   capture across the three studies in the LEVEL_collation archive.
%
%   Walks RVC_Daley_2008-2011, Blum_Surface_2012-06 and Blum_DropVsPothole_2012-02
%   for paired marker (.tsv) + force (.txt) recordings, dedupes RVC trials that
%   appear in both Level_running and LevelforObstacleRunning, parses
%   study/date/bird/trial from the path and filename, seeds body mass from
%   _metadata/AllBirds_metadata.csv (definitive mass is recomputed from force in
%   the batch), and writes _RAnalysis/data/trialRoster_multi.csv.
%
%   Bird IDs are study-tagged (rvc_/surf_/pot_ + colour code) so colour codes
%   reused across cohorts stay distinct; a conservative cross-study subject
%   grouping is layered on later for the statistics. The 2011-11-15 clutter
%   session is skipped (only .qtm projects, no exported .tsv markers).
%
%   See also: GaitSelMulti_ImportBout, GaitSelMulti_BatchProcess.

    if nargin < 1, opts = struct(); end
    P = projectPaths(); addpath(P.helpers);
    % Raw level-trial archive. projectPaths looks for LEVEL_collation under the project
    % root and then one level above it, accepting a candidate only if it holds the
    % OtherSpecies_byDate/Guinea fowl tree. Override with opts.gfRoot to point at that
    % tree directly; CROOT is then its grandparent, so the recorded paths stay relative
    % to the archive either way.
    GFROOT = getOpt(opts,'gfRoot', '');
    if isempty(GFROOT)
        if strlength(P.collation) == 0
            error('GaitSelMulti:noCollation', ...
                ['LEVEL_collation not found. Looked for a folder containing ' ...
                 'OtherSpecies_byDate/Guinea fowl at:\n    %s\n    %s\n' ...
                 'Put the archive in either place, or pass ' ...
                 'GaitSelMulti_BuildRoster(struct(''gfRoot'', <path to Guinea fowl>)).'], ...
                P.collationTried{1}, P.collationTried{2});
        end
        GFROOT = char(P.gfRoot);
        CROOT  = char(P.collation);
    else
        GFROOT = char(GFROOT);
        CROOT  = fileparts(fileparts(GFROOT));
    end
    fprintf('GaitSelMulti_BuildRoster: archive at %s\n', CROOT);
    metaCsv = fullfile(fileparts(fileparts(GFROOT)), '_metadata', 'AllBirds_metadata.csv');

    studies = {
        'RVC_Daley_2008-2011',        'rvc',  'RVC_Daley (morphometric ref)'
        'Blum_Surface_2012-06',       'surf', 'Blum_Surface_2012-06'
        'Blum_DropVsPothole_2012-02', 'pot',  'Blum_DropVsPothole_2012-02'};

    massMap = loadMasses(metaCsv);           % containers.Map: 'cohort|gf_code' -> kg

    rows = {}; seen = containers.Map('KeyType','char','ValueType','logical');
    for s = 1:size(studies,1)
        studyDir = fullfile(GFROOT, studies{s,1});
        tag = studies{s,2}; cohort = studies{s,3};
        dd = dir(fullfile(studyDir, '**', '*.tsv'));
        for i = 1:numel(dd)
            if contains(dd(i).folder, [filesep '_metadata']), continue; end
            tsv = fullfile(dd(i).folder, dd(i).name);
            txt = regexprep(tsv, '\.tsv$', '.txt');
            if ~exist(txt,'file'), continue; end                 % need force + markers
            log = regexprep(tsv, '\.tsv$', '_logfile.txt');
            if ~exist(log,'file'), log = ''; end

            [~, stem] = fileparts(dd(i).name);                   % gf_<code>_<lvl>_<NNNN>
            tok = regexp(stem, '^gf_([A-Za-z0-9]+)_(.+)$', 'tokens', 'once');
            if isempty(tok), continue; end
            code = lower(tok{1}); rest = tok{2};                 % rest e.g. L00cm_0017
            dateStr = regexp(dd(i).folder, '(\d{4})-(\d{2})-(\d{2})', 'tokens', 'once');
            if isempty(dateStr), continue; end
            dateCode = str2double([dateStr{1} dateStr{2} dateStr{3}]);
            num = regexp(rest, '(\d+)\s*$', 'tokens', 'once');
            trialNo = NaN; if ~isempty(num), trialNo = str2double(num{1}); end

            key = sprintf('%s|%d|%s|%s', tag, dateCode, code, rest);   % dedup key
            if isKey(seen,key), continue; end                    % RVC dup across level folders
            seen(key) = true;

            bird = sprintf('%s_%s', tag, code);
            mkey = sprintf('%s|gf_%s', cohort, code);
            massSeed = NaN; if isKey(massMap,mkey), massSeed = massMap(mkey); end

            r = struct();
            r.boutID      = string(sprintf('%s_%d_%s', bird, dateCode, rest));
            r.bird        = string(bird);
            r.colourCode  = string(individualCode(code));
            r.study       = string(studies{s,1});
            r.dateCode    = dateCode;
            r.trialNo     = trialNo;
            r.trialStem   = string(rest);
            r.bodyMass_kg = massSeed;
            r.windowMode  = "auto";
            % Recorded relative to the archive root, with forward slashes, so the roster
            % is portable and names no home directory. resolveRawPath rejoins them to
            % wherever LEVEL_collation is found on the reading machine.
            r.forceFile   = relToArchive(txt, CROOT);
            r.markerFile  = relToArchive(tsv, CROOT);
            r.logFile     = relToArchive(log, CROOT);
            rows{end+1} = r; %#ok<AGROW>
        end
    end

    R = struct2table([rows{:}]);
    R = sortrows(R, {'study','bird','dateCode','trialNo'});
    if ~exist(P.rData,'dir'); mkdir(P.rData); end
    outCsv = fullfile(P.rData, 'trialRoster_multi.csv');
    writetable(R, outCsv);

    fprintf('GaitSelMulti_BuildRoster: %d level trials with raw capture.\n', height(R));
    [g,gn] = findgroups(R.study);
    for k = 1:numel(gn)
        nb = numel(unique(R.bird(g==k)));
        fprintf('  %-28s %3d trials, %d birds\n', gn(k), sum(g==k), nb);
    end
    nMass = sum(~isnan(R.bodyMass_kg));
    fprintf('  mass seed found for %d/%d trials; wrote %s\n', nMass, height(R), outCsv);
end

% -------------------------------------------------------------------- helpers
function rel = relToArchive(absPath, croot)
%RELTOARCHIVE  Path relative to the LEVEL_collation root, forward-slashed.
    if isempty(absPath), rel = ""; return; end
    s = string(absPath); c = string(croot);
    if ~endsWith(c, filesep), c = c + filesep; end
    if startsWith(s, c)
        rel = extractAfter(s, strlength(c));
    else
        rel = s;                       % outside the archive: record as found
    end
    rel = replace(rel, "\", "/");
end

function m = loadMasses(csvFile)
    m = containers.Map('KeyType','char','ValueType','double');
    if ~exist(csvFile,'file'); warning('GaitSelMulti:noMeta','no %s', csvFile); return; end
    T = readtable(csvFile, 'TextType','string', 'VariableNamingRule','preserve');
    for i = 1:height(T)
        if ~strcmpi(strtrim(T.Species(i)),"Guinea fowl"), continue; end
        key = sprintf('%s|%s', strtrim(T.Cohort(i)), strtrim(T.BirdCode(i)));
        v = double(T.Mass_kg(i));
        if isfinite(v), m(char(key)) = v; end
    end
end

function c = individualCode(code)
%INDIVIDUALCODE  The colour code of the ANIMAL, where the recorded label names it differently.
%
%   The June 2012 surface study recorded one bird as `noc`, a bird carrying no colour band. It is
%   the blue bird, `blu`, and the two records that establish it are in the study's own
%   BirdBodyMasses.xlsx:
%
%   1. Its WeightRecords sheet lists the seven birds run, weighed before and after each day:
%      red, blue, black, yellow, green-yellow, red-red and green. The archive holds seven birds
%      for the same study, red, noc, bla, yel, gry, rre and g00, and the only difference between
%      the two lists is blue against noc. Two birds were run on 26 June, blue and black, and the
%      archive holds exactly two for that date, noc and bla.
%   2. Its CoM_offsets sheet gives the original analysis's own per-bird offsets for that same
%      seven-bird set, using the label `noc`. Their vertical offset follows the lab's anatomical
%      rule -0.0525*m^(1/3), so it inverts to the mass they used: every one of the seven returns
%      that bird's own weighed pre/post mean to the gram, and `noc` returns 1.3305 kg, which is
%      blue's mean on 27 June. No other bird in the record weighed that.
%
%   The recorded label is kept as the bird and session code, because it names the raw files. Only
%   the colour code changes, which is what the individual grouping and the random effect use, so
%   the surface trials join the same individual as the February 2012 blue bird.
    c = string(code);
    if c == "noc", c = "blu"; end
end

function v = getOpt(o,n,d)
    if isstruct(o) && isfield(o,n) && ~isempty(o.(n)); v = o.(n); else; v = d; end
end
