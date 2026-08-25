function P = projectPaths()
%PROJECTPATHS  Self-locating folder map for the guinea fowl gait-selection analysis.
%
%   P = projectPaths() returns a struct of absolute paths derived from the location of this
%   file, so the pipeline runs unchanged on any machine and nothing is hard-coded.
%
%   Fields:
%     P.root      project root (parent of _MatlabProcessingCode)
%     P.code      _MatlabProcessingCode
%     P.helpers   _MatlabProcessingCode/helpers    (biomechanics kernels)
%     P.rData     _RAnalysis/data                  (tidy CSV hand-off to R)
%     P.saved     _MatlabProcessingCode/SavedResults (per-step measures)
%     P.dryad     GaitSel_DryadPackage_AllGF       (the downloaded data package)
%     P.rerun     rerun                            (outputs of a re-run from the package)
%
%   THE DATA PACKAGE IS THE INPUT. Download it from Dryad, doi:10.5061/dryad.7m0cfxqc4, unzip
%   per_trial_timeseries.zip and mean_stepcycle_by_gait_speed.zip inside it, and put the package
%   at P.dryad. Every step of this pipeline reads from there; nothing writes back into it, so a
%   second run starts from the archived state.
%
%   See also: GaitSelMulti_BatchFromDryad.

    here      = fileparts(mfilename('fullpath'));
    P.root    = fileparts(here);
    P.code    = here;
    P.helpers = fullfile(here, 'helpers');
    P.rData   = fullfile(P.root, '_RAnalysis', 'data');
    P.saved   = fullfile(P.code, 'SavedResults');
    P.dryad   = fullfile(P.root, 'GaitSel_DryadPackage_AllGF');
    P.rerun   = fullfile(P.root, 'rerun');

    for f = {'rData','saved','rerun'}
        if ~exist(P.(f{1}), 'dir'); mkdir(P.(f{1})); end
    end
end
