function P = projectPaths()
%PROJECTPATHS  Self-locating folder map for the GF gait-selection project.
%
%   P = projectPaths() returns a struct of absolute paths derived from the
%   location of this file, so the pipeline runs unchanged on any machine
%   (no hard-coded absolute paths).
%
%   Fields:
%     P.root      project root (parent of _MatlabProcessingCode)
%     P.code      _MatlabProcessingCode
%     P.helpers   _MatlabProcessingCode/helpers   (vendored kernels)
%     P.rData     _RAnalysis/data                  (tidy CSV hand-off to R)
%     P.saved     _MatlabProcessingCode/SavedResults (import cache, batch measures)
%     P.internal  _Internal/pipeline_qc            (QC diagnostics; absent from the
%                                                   publication bundle, so it is not
%                                                   created here)
%     P.collation LEVEL_collation, the raw force/marker archive, or "" when it is
%                 not present. See below.
%     P.gfRoot    the guinea fowl tree inside it, LEVEL_collation/OtherSpecies_byDate/
%                 Guinea fowl, or "" when P.collation is ""
%
%   LOCATING THE RAW ARCHIVE. LEVEL_collation is large and is not always kept
%   inside the project, so it is searched for in two places, in order: under
%   P.root, then one level above it (a sibling of the project folder, which is
%   how it is arranged when several projects share one archive). A candidate is
%   accepted only if it contains OtherSpecies_byDate/Guinea fowl, so a stray
%   empty directory of the right name is not mistaken for the archive.
%
%   P.collation is "" when neither candidate confirms, rather than an error,
%   because most of the pipeline runs from the tidy tables and never touches the
%   raw files. The two functions that do need it, GaitSelMulti_BuildRoster and
%   GaitSelMulti_ImportBout, raise their own error naming both candidate paths.
%
%   Every raw file is recorded in the roster by its path RELATIVE to
%   P.collation, and resolved at read time by resolveRawPath against whichever
%   of the two locations confirmed on this machine. Nothing downstream stores an
%   absolute path, so the roster is portable and carries no home directory.
%
%   See also: GaitSelMulti_RunFullWorkflow, GaitSelMulti_BuildRoster, resolveRawPath.

    here      = fileparts(mfilename('fullpath'));
    P.root    = fileparts(here);
    P.code    = here;
    P.helpers = fullfile(here, 'helpers');
    P.rData   = fullfile(P.root, '_RAnalysis', 'data');
    P.saved   = fullfile(P.code, 'SavedResults');
    P.internal= fullfile(P.root, '_Internal', 'pipeline_qc');

    for f = {'rData','saved'}
        if ~exist(P.(f{1}), 'dir'); mkdir(P.(f{1})); end
    end

    [P.collation, P.gfRoot, P.collationTried] = locateCollation(P.root);
end

% -------------------------------------------------------------------- helpers
function [croot, gfroot, tried] = locateCollation(projRoot)
%LOCATECOLLATION  Find LEVEL_collation under the project, else one level above.
    marker = fullfile('OtherSpecies_byDate', 'Guinea fowl');
    tried  = { fullfile(projRoot, 'LEVEL_collation'), ...
               fullfile(fileparts(projRoot), 'LEVEL_collation') };
    croot = ""; gfroot = "";
    for k = 1:numel(tried)
        if isfolder(fullfile(tried{k}, marker))
            croot  = string(tried{k});
            gfroot = string(fullfile(tried{k}, marker));
            return
        end
    end
end
