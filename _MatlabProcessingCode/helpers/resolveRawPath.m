function abspath = resolveRawPath(stored, collationRoot)
%RESOLVERAWPATH  Turn a roster file reference into a path on this machine.
%
%   abspath = resolveRawPath(stored, collationRoot) resolves one of the roster's
%   forceFile / markerFile / logFile entries against the LEVEL_collation archive
%   located by projectPaths.
%
%   The roster stores each raw file by its path relative to LEVEL_collation, for
%   example
%       OtherSpecies_byDate/Guinea fowl/Blum_Surface_2012-06/2012-06-25/gf_red_LE0cm_0006.txt
%   which is joined to collationRoot here. Two other forms are accepted so that a
%   roster written before this convention still resolves:
%
%     * a project-relative path beginning LEVEL_collation/, from which the leading
%       segment is dropped before joining; and
%     * an absolute path, used as given when it exists on this machine, and
%       otherwise re-resolved from its LEVEL_collation segment onward, which is
%       what makes a roster built on another machine usable here.
%
%   An empty or missing reference returns "" rather than erroring, because logFile
%   is legitimately absent for the trials that have no logfile.
%
%   See also: projectPaths, GaitSelMulti_BuildRoster, GaitSelMulti_ImportBout.

    abspath = "";
    if nargin < 1 || isempty(stored), return; end
    s = string(stored);
    if ismissing(s) || strlength(strtrim(s)) == 0, return; end
    s = strtrim(s);

    % Written on Windows, read on a POSIX machine or the reverse.
    s = replace(s, "\", "/");
    croot = "";
    if nargin >= 2 && ~isempty(collationRoot), croot = string(collationRoot); end

    isAbs = startsWith(s, "/") || ~isempty(regexp(char(s), '^[A-Za-z]:', 'once'));
    if isAbs
        native = fullfile(char(s));
        if isfile(native), abspath = string(native); return; end
        % Absolute but not on this machine: keep only the archive-relative tail.
        s = tailFromCollation(s);
        if strlength(s) == 0, abspath = string(native); return; end
    else
        s = tailFromCollation(s);
    end

    if strlength(croot) == 0
        abspath = string(fullfile(char(s)));    % no archive located; best effort
    else
        abspath = string(fullfile(char(croot), char(s)));
    end
end

% -------------------------------------------------------------------- helpers
function tail = tailFromCollation(s)
%TAILFROMCOLLATION  Drop everything up to and including a LEVEL_collation segment.
%   Leaves an already archive-relative path unchanged.
    tail = s;
    parts = split(s, "/");
    hit = find(parts == "LEVEL_collation", 1, 'last');
    if ~isempty(hit)
        if hit == numel(parts)
            tail = "";
        else
            tail = strjoin(parts(hit+1:end), "/");
        end
    end
end
