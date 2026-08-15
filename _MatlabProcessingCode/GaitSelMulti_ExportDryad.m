function GaitSelMulti_ExportDryad(opts)
%GAITSELMULTI_EXPORTDRYAD  Compiled multi-study guinea fowl Dryad package.
%   Per-trial whole-bout time series for every level trial across the three
%   studies, written as a .mat and .csv per bout in a STANDARDIZED axis frame, plus
%   trial and step index tables and a README. Force and kinematics are oriented
%   bird-centered and right-handed: columns are [medio-lateral, fore-aft, vertical]
%   with fore-aft POSITIVE IN THE DIRECTION OF TRAVEL for every trial (trials where
%   the bird crossed the other way have their fore-aft and medio-lateral signs
%   flipped together, preserving handedness), and vertical positive up. Units SI
%   (force N, position m, velocity m/s, time s).
%
%   Reconstruction is re-run from the cached multi imports, reassembling each
%   bout's force with the per-(study,date) consensus fore-aft axis the batch used
%   (faAxisUsed + per-trial sign from the stored horizontal-impulse correlations).
%   Per-step gait class / outlier / dimensionless speed come from the tidy tables.
%
%   See also: GaitSelMulti_BatchProcess, GaitSelMulti_ExportCycleTraces, GaitSel_ExportDryadMeanCycles.
    if nargin < 1, opts = struct(); end
    P = projectPaths(); addpath(P.helpers); g = 9.81;
    root  = getOpt(opts,'root', fullfile(P.root, 'GaitSel_DryadPackage_AllGF'));
    % indexOnly rebuilds the index tables and leaves the per-trial time series in place. The time
    % series are raw reconstructions and do not depend on the gait classification, so a relabel
    % does not require rewriting them.
    indexOnly = getOpt(opts,'indexOnly', false);
    tsDir = fullfile(root, 'per_trial_timeseries');
    if ~exist(tsDir,'dir'), mkdir(tsDir); end

    R  = readtable(fullfile(P.rData,'trialRoster_multi.csv'), 'TextType','string', ...
                   'Delimiter',',');
    C  = load(fullfile(P.saved,'importCache_multi.mat'), 'bouts','calCx','calCy');
    M  = load(fullfile(P.saved,'perStepMeasures_multi.mat'), 'massBird','faAxisUsed');
    bouts = C.bouts; massBird = M.massBird; faAxisUsed = M.faAxisUsed;
    calCx = C.calCx; calCy = C.calCy;

    PS = readtable(fullfile(P.rData,'perStep_long_multi.csv'), 'TextType','string');
    MO = readtable(fullfile(P.rData,'morphology_multi.csv'), 'TextType','string');
    L0map = containers.Map(cellstr(string(MO.session)), num2cell(MO.L0_m));
    keyf  = @(b,s) sprintf('%s|%d', char(b), double(s));
    labG  = containers.Map('KeyType','char','ValueType','any');
    for i = 1:height(PS)
        labG(keyf(PS.boutID(i),PS.stepIndex(i))) = struct('gait',char(PS.gaitObjective(i)), ...
            'outlier',double(PS.outlier(i)), 'u',PS.meanSpeed_n(i));
    end
    % The gait label and the analysis gate live in R (02_clean.R and 04_analysis_sample.R), so the
    % package takes both from data/step_labels_analysis.csv rather than from this table's
    % gaitObjective column, which is the hodograph-based label rather than the classification the
    % paper uses, and rather than from
    % step_qcpass.csv alone, which flags the steps passing the force-based QC gate where the paper
    % analyses the smaller set that also has every gait-space descriptor defined. Without this the
    % package and the paper would disagree on both the labels and the sample size.
    labA = containers.Map('KeyType','char','ValueType','any');
    laf = fullfile(P.rData,'step_labels_analysis.csv');
    if exist(laf,'file')
        LA = readtable(laf,'Delimiter',',','TextType','string');
        for i = 1:height(LA)
            labA(keyf(LA.boutID(i),LA.stepIndex(i))) = struct( ...
                'gait',safechar(LA.gait(i)), 'gait4',safechar(LA.gait4(i)), ...
                'steadiness',safechar(LA.steadiness(i)), ...
                'inSample',double(LA.analysisSample(i)), 'why',safechar(LA.exclusionReason(i)));
        end
    else
        warning('step_labels_analysis.csv not found; run the R pipeline before exporting.');
    end
    % per-step QC-retention flag: read the SINGLE SOURCE OF TRUTH written by 02_clean.R
    % (step_qcpass.csv) so steps_index qcPass reproduces the analyzed set exactly, rather
    % than re-implementing the gate here.
    labQ = containers.Map('KeyType','char','ValueType','double');
    qcf = fullfile(P.rData,'step_qcpass.csv');
    if exist(qcf,'file')
        QC = readtable(qcf,'Delimiter',',','TextType','string');
        for i = 1:height(QC), labQ(keyf(QC.boutID(i),QC.stepIndex(i))) = double(QC.qcPass(i)); end
    end
    % per-step energy-based steadiness (from the R clean stage) and the reconciled
    % individual ID (random-effect grouping: RVC cohort, and the two 2012 Blum studies
    % pooled as one cohort by colour code; per-session mass/L0 are unchanged).
    labS = containers.Map('KeyType','char','ValueType','char');
    ssf = fullfile(P.rData,'step_steadiness.csv');
    if exist(ssf,'file')
        % explicit comma delimiter: this small file's boutIDs contain underscores, so
        % readtable's delimiter auto-detection otherwise mis-splits it on '_'.
        SS = readtable(ssf,'Delimiter',',','TextType','string');
        for i = 1:height(SS), labS(keyf(SS.boutID(i),SS.stepIndex(i))) = char(SS.accClass(i)); end
    end

    trialRows = {}; stepRows = {}; nOK = 0;
    for i = 1:numel(bouts)
        B = bouts{i}; if isempty(B), continue; end
        % massBird is keyed birdCode|YYYY-MM, the session key the batch, the offsets and L0 share
        bird = char(B.bird);
        sk = char(string(B.bird) + "|" + sessionMonth(B.dateCode));
        if ~isKey(massBird,sk), continue; end
        % reassemble force with the batch's (study,date) consensus fore-aft axis
        Fx = B.forceRawHoriz(:,1); Fy = B.forceRawHoriz(:,2); Fz = B.forceVert;
        if i<=numel(faAxisUsed) && faAxisUsed(i)=="X"
            B.force = [Fy, sgn(calCx(i))*Fx, Fz];
        else
            B.force = [Fx, sgn(calCy(i))*Fy, Fz];
        end
        try
            B.mass = massBird(sk);
            B = GaitSel_ReconstructCoM(B);
            E = GaitSel_DetectGaitEvents(B);
        catch ME
            warning('GaitSelMulti_ExportDryad: skip %s (%s)', char(B.boutID), ME.message); continue;
        end

        t = B.time(:); n = numel(t);
        F=B.force; com=B.com; comP=B.comProxy; vel=B.comVel;
        % Foot tracks are already spike-cleaned, once, in GaitSelMulti_ImportBout, so the series
        % published here are the same feet the analysis measured its virtual leg from.
        fR=B.footR; fL=B.footL;
        % standardize orientation: fore-aft positive in travel; flip ml with it to
        % keep a right-handed bird-centered frame. Vertical unchanged.
        travelDir = sign(com(end,2)-com(1,2)); if travelDir==0, travelDir = 1; end
        if travelDir < 0
            F(:,1:2)=-F(:,1:2); com(:,1:2)=-com(:,1:2); comP(:,1:2)=-comP(:,1:2);
            vel(:,1:2)=-vel(:,1:2); fR(:,1:2)=-fR(:,1:2); fL(:,1:2)=-fL(:,1:2);
        end

        si=E.stepIdx(:); nStep=max(0,numel(si)-1);
        stepOf=zeros(n,1); for k=1:nStep, stepOf(si(k):min(si(k+1)-1,n))=k; end
        bid=char(B.boutID); L0=0; if isKey(L0map,sk), L0=L0map(sk); end
        driftVmm=NaN; if numel(B.driftRMS)>=3, driftVmm=B.driftRMS(3)*1000; end

        Tt = table(t, F(:,1),F(:,2),F(:,3), com(:,1),com(:,2),com(:,3), ...
            comP(:,1),comP(:,2),comP(:,3), vel(:,1),vel(:,2),vel(:,3), ...
            fR(:,1),fR(:,2),fR(:,3), fL(:,1),fL(:,2),fL(:,3), stepOf, ...
            'VariableNames',{'t_s','F_ml_N','F_fa_N','F_vt_N', ...
              'com_ml_m','com_fa_m','com_vt_m','comProxy_ml_m','comProxy_fa_m','comProxy_vt_m', ...
              'comVel_ml_ms','comVel_fa_ms','comVel_vt_ms', ...
              'footR_ml_m','footR_fa_m','footR_vt_m','footL_ml_m','footL_fa_m','footL_vt_m','stepIndex'});
        if ~indexOnly, writetable(Tt, fullfile(tsDir,[bid '.csv'])); end

        sI=(1:nStep)'; sS=si(1:nStep); sEn=si(2:nStep+1);
        cc = char(R.colourCode(i));
        sid = sprintf('%s_%s', string(iff(contains(char(B.study),'RVC'),"rvc","blum12")), cc);
        sGait=strings(nStep,1); sOut=nan(nStep,1); sU=nan(nStep,1); sStead=strings(nStep,1);
        sQC=nan(nStep,1); sG4=strings(nStep,1); sIn=nan(nStep,1); sWhy=strings(nStep,1);
        for k=1:nStep
            kk=keyf(bid,k);
            if isKey(labG,kk), L=labG(kk); sOut(k)=L.outlier; sU(k)=L.u; end
            if isKey(labQ,kk), sQC(k)=labQ(kk); end       % qcPass from 02_clean.R (single source)
            if isKey(labS,kk), sStead(k)=string(labS(kk)); end
            % gait, four-cell gait and the analysis-sample flag come from the R label file
            if isKey(labA,kk)
                A=labA(kk); sGait(k)=string(A.gait); sG4(k)=string(A.gait4);
                sStead(k)=string(A.steadiness); sIn(k)=A.inSample; sWhy(k)=string(A.why);
            else
                sIn(k)=0; sWhy(k)="did not pass the QC gate";
            end
        end
        stepsTab = table(repmat(string(bid),nStep,1), sI, sS, sEn, t(sS), t(min(sEn,n)), sGait, sG4, sStead, sOut, sQC, sIn, sWhy, sU, ...
            'VariableNames',{'boutID','stepIndex','startSample','endSample','startTime_s','endTime_s','gait','gait4','steadiness','outlier','qcPass','analysisSample','exclusionReason','u'});

        trial = struct('boutID',bid,'bird',bird,'subjectID',sid,'study',char(B.study),'dateCode',B.dateCode, ...
            'colourCode',cc,'mass_kg',B.mass,'L0_m',L0,'bodyWeight_N',B.mass*g, ...
            'fs_force_Hz',B.fHz,'fs_kin_Hz',B.kinHz,'g',g,'nPlates',B.nPlates, ...
            'axisOrder','columns = [medio-lateral, fore-aft (=travel, positive forward), vertical (up)]', ...
            'foreaftAxis',char(faAxisUsed(min(i,numel(faAxisUsed)))),'driftRMS_vert_mm',driftVmm, ...
            'WE_identity_r',B.WEcheck.r,'time_s',t,'GRF_N',F,'com_m',com,'comProxy_m',comP, ...
            'comVel_ms',vel,'footR_m',fR,'footL_m',fL,'stepIndexPerSample',stepOf,'steps',stepsTab);
        if ~indexOnly, save(fullfile(tsDir,[bid '.mat']),'trial','-v7.3'); end

        stepRows{end+1}=stepsTab; %#ok<AGROW>
        trialRows{end+1}=table(string(bid),string(B.study),B.dateCode,string(bird),string(cc),string(sid), ...
            B.mass,L0,B.nPlates,n,t(end)-t(1),nStep,driftVmm,B.WEcheck.r, string(faAxisUsed(min(i,numel(faAxisUsed)))), ...
            'VariableNames',{'boutID','study','dateCode','bird','colourCode','subjectID','mass_kg','L0_m','nPlates', ...
              'nSamples','duration_s','nSteps','driftRMS_vert_mm','WE_identity_r','foreaftAxis'}); %#ok<AGROW>
        nOK=nOK+1;
        if mod(nOK,40)==0, fprintf('  exported %d\n', nOK); end
    end
    if ~isempty(trialRows), writetable(vertcat(trialRows{:}), fullfile(root,'trial_index.csv')); end
    if ~isempty(stepRows),  writetable(vertcat(stepRows{:}),  fullfile(root,'steps_index.csv')); end
    % copy the tidy analysis tables into the package
    for f = {'perStep_long_multi.csv','perStride_long_multi.csv','morphology_multi.csv','trialRoster_multi.csv'}
        src = fullfile(P.rData,f{1}); if exist(src,'file'), copyfile(src, fullfile(root,f{1})); end
    end
    % GaitSelMulti_BuildRoster records the raw files relative to LEVEL_collation. The package is
    % published, so its copy of the roster is normalised to the archive-relative form regardless of
    % how the roster was written; an absolute path would name a home directory.
    rosterOut = fullfile(root,'trialRoster_multi.csv');
    if exist(rosterOut,'file')
        RT = readtable(rosterOut, 'TextType','string', 'Delimiter',',');
        for c = ["forceFile","markerFile","logFile"]
            if ismember(c, string(RT.Properties.VariableNames))
                RT.(c) = arrayfun(@archiveRelative, RT.(c));
            end
        end
        writetable(RT, rosterOut);
    end
    fprintf('GaitSelMulti_ExportDryad: %d trials -> %s\n', nOK, tsDir);
end

function s = sgn(x), if ~isfinite(x)||x>=0, s=1; else, s=-1; end, end

function m = sessionMonth(dateCode)
%SESSIONMONTH  YYYY-MM from a YYYYMMDD date code; the session key mass and L0 share.
    t = regexprep(string(dateCode), '\D', '');
    assert(all(strlength(t) >= 6), 'GaitSelMulti:badDateCode', ...
           'a dateCode has no YYYYMM to key the session on');
    m = extractBetween(t,1,4) + "-" + extractBetween(t,5,6);
end

function rel = archiveRelative(p)
%ARCHIVERELATIVE  Keep only the part of a path below LEVEL_collation.
    rel = "";
    if ismissing(p) || strlength(strtrim(p)) == 0, return; end
    rel = replace(strtrim(p), "\", "/");
    parts = split(rel, "/");
    hit = find(parts == "LEVEL_collation", 1, 'last');
    if ~isempty(hit) && hit < numel(parts)
        rel = strjoin(parts(hit+1:end), "/");
    end
end
function o = iff(c,a,b), if c, o=a; else, o=b; end, end
function v = getOpt(o,n,d), if isstruct(o)&&isfield(o,n)&&~isempty(o.(n)), v=o.(n); else, v=d; end, end

function s = safechar(x)
% string -> char, with <missing> becoming empty rather than erroring
    if ismissing(x), s = ''; else, s = char(x); end
end
