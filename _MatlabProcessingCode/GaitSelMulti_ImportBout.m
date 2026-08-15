function B = GaitSelMulti_ImportBout(row, opts)
%GAITSELMULTI_IMPORTBOUT  Raw importer for the multi-study GF collation.
%
%   Reads one level trial's raw Kistler force (.txt, already in Newtons) and
%   Qualisys marker (.tsv) files and returns the bout struct B consumed by the
%   reconstruction / detection / measure core. It handles both rig generations,
%   RVC_Daley (5 plates + sync column, 7 markers incl. Head) and the Blum 2012
%   studies (6 plates, no sync column, 6 markers, no Head), by:
%
%     (1) active-plate count and column layout come from the paired _logfile.txt
%         when present (Plate N ON/OFF lines), else inferred as floor(nCol/8) with
%         a trailing sync column detected as nCol == nPlates*8 + 1;
%     (2) the CoM-proxy and foot markers are selected BY NAME from the .tsv
%         MARKER_NAMES header (Cranial+Caudal = CoM proxy; Right_TMP+Right_digit_3
%         and Left_TMP+Left_digit_3 = feet), so a present or absent Head marker
%         does not shift column indices;
%     (3) the analysis window is the whole tracked trial by default: the span
%         between the first and last frame where the CoM-proxy markers are tracked
%         (auto-window), rather than a curated StartFrame/EndFrame.
%
%   Per-plate channel order (8 cols): Fx12 Fx34 Fy14 Fy23 Fz1 Fz2 Fz3 Fz4. Net
%   force axes [ml fa vert]: Fx (X12+X34) -> ml, Fy (Y14+Y23) -> fa, Fz -> vert.
%
%   row fields: forceFile, markerFile, (logFile), bird, boutID, study, dateCode,
%   trialNo, bodyMass_kg, (fHz, kinHz), (dirML,dirFA,dirVT), (windowMode, startFrame,
%   endFrame). opts: forceCutoffHz (default 50), markerAxisMap, markerSign, etc.
%
%   See also: GaitSelMulti_BuildRoster, GaitSelMulti_BatchProcess, GaitSel_ReconstructCoM.

    if nargin < 2, opts = struct(); end
    axisMap = getOpt(opts,'markerAxisMap',[2 1 3]);   % marker (X,Y,Z) -> [ml fa vert]
    axisSgn = getOpt(opts,'markerSign',[1 1 1]);
    unitMode= getOpt(opts,'markerUnit','auto');
    smTol   = getOpt(opts,'smoothProxy',0);
    fCut    = getOpt(opts,'forceCutoffHz',50);

    row = table2struct2(row);
    P = projectPaths(); addpath(P.helpers);

    % ---- resolve the raw files against the archive on THIS machine ------
    % The roster records each file relative to LEVEL_collation, so it is portable;
    % resolveRawPath rejoins it to wherever projectPaths found the archive.
    row.forceFile  = resolveRawPath(getField(row,'forceFile',''),  P.collation);
    row.markerFile = resolveRawPath(getField(row,'markerFile',''), P.collation);
    row.logFile    = resolveRawPath(getField(row,'logFile',''),    P.collation);
    for f = ["forceFile","markerFile"]
        if strlength(row.(f)) == 0 || ~isfile(row.(f))
            if strlength(P.collation) == 0
                error('GaitSelMulti:noCollation', ...
                    ['%s: cannot resolve %s because LEVEL_collation was not found. ' ...
                     'Looked for a folder containing OtherSpecies_byDate/Guinea fowl at:' ...
                     '\n    %s\n    %s'], getField(row,'boutID','(bout)'), f, ...
                    P.collationTried{1}, P.collationTried{2});
            end
            error('GaitSelMulti:missingRaw', '%s: %s not found at %s', ...
                getField(row,'boutID','(bout)'), f, row.(f));
        end
    end

    % ---- resolve acquisition rates + plate layout -----------------------
    logFile = row.logFile;
    [nPlates, logHz] = parseLogfile(logFile);
    fHz  = firstFinite(getField(row,'fHz',NaN), firstFinite(logHz, 500));
    kinHz= firstFinite(getField(row,'kinHz',NaN), 250);
    sr   = fHz / kinHz;
    % Anatomical [ml fa vert] force axes are auto-calibrated to the kinematics
    % below (the plate-to-marker frame relationship differs by rig), so no manual
    % direction signs are needed here.

    % ---- force: read, determine plates, baseline, sum to net GRF --------
    raw = readmatrix(char(row.forceFile), 'FileType','text', 'Delimiter','\t');
    nF  = size(raw,1); nCol = size(raw,2);
    if ~isfinite(nPlates) || nPlates < 1
        nPlates = floor(nCol/8);                       % infer if no logfile
    end
    nChan = nPlates*8;
    if nChan > nCol, error('GaitSelMulti:force','%s: %d cols < %d plate channels', row.forceFile, nCol, nChan); end
    % baseline-zero each corner channel by its full-record median (plate unloaded
    % ~99% of the record). Any trailing column(s) beyond nChan (e.g. a sync column
    % in the 5-plate RVC files, nCol=41) are ignored.
    raw(:,1:nChan) = raw(:,1:nChan) - median(raw(:,1:nChan), 1, 'omitnan');
    FxP = zeros(nF,nPlates); FyP = zeros(nF,nPlates); FzP = zeros(nF,nPlates);
    for p = 1:nPlates
        b = (p-1)*8;
        FxP(:,p) = raw(:,b+1) + raw(:,b+2);            % x12 + x34 (plate X axis)
        FyP(:,p) = raw(:,b+3) + raw(:,b+4);            % y14 + y23 (plate Y axis)
        FzP(:,p) = raw(:,b+5) + raw(:,b+6) + raw(:,b+7) + raw(:,b+8);  % z1..z4
    end
    % Net force in the raw plate frame: horizontal X, horizontal Y, vertical Z. The
    % plate frame is oriented differently across rigs (the 2009 5-plate and 2012
    % 6-plate setups have their horizontal axes swapped and sign-flipped relative
    % to the Qualisys marker frame), so the anatomical [ml fa vert] assignment of
    % the two horizontal axes is NOT fixed. It is auto-calibrated below against the
    % kinematic travel direction. Vertical is unambiguous (Z, positive up).
    Fx = sum(FxP,2); Fy = sum(FyP,2); Fz = sum(FzP,2);
    if mean(Fz,'omitnan') < 0, Fz = -Fz; FzP = -FzP; end   % ensure GRF positive up
    horiz = [Fx, Fy];
    if isfinite(fCut) && fCut > 0 && fCut < fHz/2
        horiz = filterForce(horiz, fHz, struct('CutoffHz', fCut));
        Fz    = filterForce(Fz,    fHz, struct('CutoffHz', fCut));
        FzP   = filterForce(FzP,   fHz, struct('CutoffHz', fCut));
    end
    Fx = horiz(:,1); Fy = horiz(:,2);

    % ---- markers: read, resolve by name, build CoM proxy ----------------
    [M, names] = readTsv(char(row.markerFile));
    nM = size(M,1);
    iCran = markerIdx(names,{'Cranial','cranial'});
    iCaud = markerIdx(names,{'Caudal','caudal'});
    iRT   = markerIdx(names,{'Right_TMP','RightTMP','R_TMP','Right_TMT'});
    iRD   = markerIdx(names,{'Right_digit_3','Right_digit3','R_digit_3','Right_Digit_3'});
    iLT   = markerIdx(names,{'Left_TMP','LeftTMP','L_TMP','Left_TMT'});
    iLD   = markerIdx(names,{'Left_digit_3','Left_digit3','L_digit_3','Left_Digit_3'});
    if any(~isfinite([iCran iCaud]))
        error('GaitSelMulti:markers','%s: could not find Cranial/Caudal in [%s]', ...
              row.markerFile, strjoin(names,', '));
    end

    % ---- whole-trial auto-window (on-plate running bout) ----------------
    % The analysis window is the span where the bird is BOTH tracked and on the
    % force plates. F = m a reconstruction is only valid while the ground reaction
    % is captured, so the tracked-marker span alone is too broad on the longer
    % capture volumes (it includes off-plate approach/departure running, where the
    % summed GRF is ~0 but the bird is still accelerating). Intersect the tracked
    % CoM-proxy span with the force-active span (first-to-last vertical contact).
    windowMode = getField(row,'windowMode','auto');
    com1full = markerXYZ(M, iCran, 1:nM);
    com2full = markerXYZ(M, iCaud, 1:nM);
    tracked  = all(isfinite(com1full),2) & all(isfinite(com2full),2);
    if strcmpi(windowMode,'auto') || ~isfinite(getField(row,'startFrame',NaN))
        idx = find(tracked);
        if isempty(idx), error('GaitSelMulti:window','%s: no tracked CoM-proxy frames', row.markerFile); end
        kt0 = idx(1); kt1 = idx(end);                          % tracked span (kin frames)
        [fa0, fa1] = forceActiveSpan(Fz, fHz);                 % on-plate span (force samples)
        ka0 = round((fa0-1)/sr) + 1; ka1 = round((fa1-1)/sr) + 1;
        ks0 = max(kt0, ka0); ks1 = min(kt1, ka1);
        if ks1 - ks0 < 5                                       % fall back to force-active only
            ks0 = max(1, ka0); ks1 = min(nM, ka1);
        end
    else
        ks0 = max(1, double(row.startFrame));
        ks1 = min(nM, double(row.endFrame));
    end
    ksel = ks0:ks1;
    tm   = (0:numel(ksel)-1)'/kinHz;

    % force window mapped from the kinematic frame window
    fs0 = max(1,  round((ks0-1)*sr) + 1);
    fs1 = min(nF, round((ks1-1)*sr) + 1);
    fsel = fs0:fs1;
    FxW = Fx(fsel); FyW = Fy(fsel); FzW = Fz(fsel);
    FzPlate = FzP(fsel,:);
    t = (0:numel(fsel)-1)'/fHz;

    com1 = markerXYZ(M, iCran, ksel);
    com2 = markerXYZ(M, iCaud, ksel);
    comProxyRaw = (com1 + com2)/2;                     % (X,Y,Z) marker frame
    footR = footPair(M, iRT, iRD, ksel);
    footL = footPair(M, iLT, iLD, ksel);

    scale = unitScale(unitMode, comProxyRaw);
    comProxyRaw = comProxyRaw*scale; footR = footR*scale; footL = footL*scale;

    comProxy_k = applyAxis(comProxyRaw, axisMap, axisSgn);
    footR_k    = applyAxis(footR,       axisMap, axisSgn);
    footL_k    = applyAxis(footL,       axisMap, axisSgn);

    % ---- shift the marker midpoint to the estimated whole-body CoM -------
    % The Cranial-Caudal midpoint sits above and slightly behind the CoM, so a per-bird offset
    % (Birn-Jeffery) shifts it to the estimated whole-body CoM, bringing the virtual leg from CoM to
    % toe and L0 into line with her measured begin-stance leg length and standing hip height.
    % Applied here, in the anatomical frame and before any geometry
    % is taken, so leg length, leg angle and hip height all see one CoM. Energy fluctuations are
    % unaffected, because a constant shift cancels in a height CHANGE.
    %
    % [ml fa vert]. The fore-aft component is stored POSITIVE IN THE DIRECTION OF TRAVEL so that
    % one value per bird applies whichever way it crossed the runway; it is converted to the lab
    % frame here by this trial's own travel direction.
    comOff = getField(row,'comOffset',[0 0 0]);
    comOff = comOff(:)';
    if numel(comOff) == 3 && any(isfinite(comOff) & comOff ~= 0)
        comOff(~isfinite(comOff)) = 0;
        travelK = sign(median(diff(comProxy_k(:,2)),'omitnan'));
        if ~isfinite(travelK) || travelK == 0, travelK = 1; end
        comOff(2) = travelK * comOff(2);
        comProxy_k = comProxy_k + comOff;
    end

    % The two back markers are kept as well as their midpoint, for the marker-placement table and
    % the offset diagnostics.
    cran_k = applyAxis(markerXYZ(M, iCran, ksel)*scale, axisMap, axisSgn);
    caud_k = applyAxis(markerXYZ(M, iCaud, ksel)*scale, axisMap, axisSgn);

    comProxy_k = fillGaps(comProxy_k, tm);
    if smTol > 0, comProxy_k = smoothProxy(comProxy_k, tm, smTol); end

    comProxy = interpToForce(comProxy_k, tm, t);
    cranial  = interpToForce(fillGaps(cran_k, tm), tm, t);
    caudal   = interpToForce(fillGaps(caud_k, tm), tm, t);
    footR    = interpToForce(footR_k,    tm, t);
    footL    = interpToForce(footL_k,    tm, t);

    % ---- auto-calibrate horizontal force to the kinematic travel axis ----
    % Assign the horizontal force axis whose cumulative impulse best tracks the
    % kinematic fore-aft velocity to fore-aft (signed so the impulse follows the
    % velocity), the other horizontal axis to medio-lateral; vertical stays Z.
    [F, cal] = assembleForce(FxW, FyW, FzW, comProxy(:,2), t);

    B = struct();
    B.boutID = string(row.boutID);  B.bird = string(row.bird);
    B.gaitLabel = string(getField(row,'gaitLabel',""));
    B.gaitCode  = getField(row,'gaitCode',NaN);
    B.trialNo   = getField(row,'trialNo',NaN);
    B.dateCode  = getField(row,'dateCode',NaN);
    B.study     = string(getField(row,'study',""));
    B.force = F;  B.time = t;  B.comProxy = comProxy;  B.FzPlate = FzPlate;
    % raw windowed horizontal axes + vertical, kept so a per-study consensus can
    % override the per-trial fore-aft axis choice (see GaitSelMulti_BatchProcess).
    B.forceRawHoriz = [FxW, FyW];  B.forceVert = FzW;
    % Foot marker tracks are spike-cleaned ONCE, here, so every consumer sees the same feet:
    % the contact detector, the virtual-leg geometry, the offset fit and the published per-trial
    % series. An untracked marker extrapolates to physically impossible positions in all three
    % axes at the same samples, so the samples are rejected on the vertical channel and dropped
    % from all three together, then gap-filled and low-passed (cleanFootTrack).
    [B.footR, B.footL] = deal(cleanFootTrack(footR, fHz), cleanFootTrack(footL, fHz));
    B.cranial = cranial; B.caudal = caudal;
    B.comOffsetApplied = comOff;
    B.fHz = fHz;  B.kinHz = kinHz;  B.mass = getField(row,'bodyMass_kg',NaN);
    B.nPlates = nPlates;  B.forceCal = cal;
    B.forceFile = string(row.forceFile); B.markerFile = string(row.markerFile);
    B.frameWin = [ks0 ks1];
    B.axisMap = axisMap; B.axisSign = axisSgn;
end

% -------------------------------------------------------------------- helpers
function [nPlates, hz] = parseLogfile(logFile)
    nPlates = NaN; hz = NaN;
    if isempty(logFile) || all(ismissing(string(logFile))) || strlength(string(logFile))==0 ...
            || ~exist(char(string(logFile)),'file')
        return;
    end
    txt = readlines(char(logFile));
    on = 0;
    for k = 1:numel(txt)
        L = strtrim(txt(k));
        % a real plate line is "Plate <number>, ON/OFF, ..."; require the numeric
        % index and an exact ON token so "Plate Orientation, Movement along ..."
        % (whose "along" contains the letters o,n) is not miscounted.
        if ~isempty(regexp(L,'^Plate\s+\d+\s*,','once'))
            parts = split(L,",");
            if numel(parts) >= 2 && strcmpi(strtrim(parts(2)),"ON"); on = on + 1; end
        elseif startsWith(L,"Sample frequency","IgnoreCase",true)
            v = str2double(extractNumber(L)); if isfinite(v), hz = v; end
        end
    end
    if on > 0, nPlates = on; end
end

function s = extractNumber(L)
    parts = split(L,","); s = strtrim(parts(end));
end

% assembleForce and corrSafe are standalone helpers (helpers/assembleForce.m,
% helpers/corrSafe.m), on the path via projectPaths.

function [i0, i1] = forceActiveSpan(Fz, fHz)
    % First-to-last force sample where the bird is on the plates. Fz is the
    % baseline-subtracted summed vertical GRF (~0 off-plate / in flight, 2-3 body
    % weights at a footfall). Threshold at a fraction of the robust peak so contact
    % onset/offset are included; require the active region to be non-trivial.
    Fz = Fz(:);
    pk = prctile(Fz(isfinite(Fz)), 99);
    if ~isfinite(pk) || pk <= 0, i0 = 1; i1 = numel(Fz); return; end
    active = Fz > max(0.12*pk, 1.0);           % 12% of peak or 1 N, whichever larger
    idx = find(active);
    if isempty(idx), i0 = 1; i1 = numel(Fz); return; end
    i0 = idx(1); i1 = idx(end);
    % trim a leading/trailing isolated spike: keep the span holding >=95% of the
    % active samples (guards against a single stray sample far from the bout).
    if (i1-i0+1) > 1.5*numel(idx)
        lo = prctile(idx, 2.5); hi = prctile(idx, 97.5);
        i0 = floor(lo); i1 = ceil(hi);
    end
    i0 = max(1, i0); i1 = min(numel(Fz), i1);
end

function [M, names] = readTsv(f)
    % Qualisys .tsv: 10 header lines; MARKER_NAMES on the line starting that token.
    lines = readlines(f);
    hdr = 10; names = strings(0);
    for k = 1:min(20,numel(lines))
        if startsWith(strtrim(lines(k)),"MARKER_NAMES","IgnoreCase",true)
            toks = split(strtrim(lines(k)), sprintf('\t'));
            names = toks(2:end)';
            hdr = k;   % marker data begins on the line after MARKER_NAMES
        end
    end
    names = cellstr(names);
    M = readmatrix(f, 'FileType','text', 'Delimiter','\t', 'NumHeaderLines', hdr);
end

function i = markerIdx(names, aliases)
    i = NaN;
    for a = 1:numel(aliases)
        hit = find(strcmpi(names, aliases{a}), 1);
        if ~isempty(hit), i = hit; return; end
    end
end

function xyz = markerXYZ(M, mIdx, sel)
    c = (mIdx-1)*3 + (1:3);
    xyz = M(sel, c);
end

function xyz = footPair(M, i1, i2, sel)
    if isfinite(i1) && isfinite(i2)
        xyz = (markerXYZ(M,i1,sel) + markerXYZ(M,i2,sel))/2;
    elseif isfinite(i1)
        xyz = markerXYZ(M,i1,sel);
    elseif isfinite(i2)
        xyz = markerXYZ(M,i2,sel);
    else
        xyz = nan(numel(sel),3);
    end
end

function s = unitScale(mode, xyz)
    switch lower(char(mode))
        case 'mm', s = 1e-3;
        case 'm',  s = 1;
        otherwise
            v = max(abs(xyz(:)), [], 'omitnan');
            if isempty(v) || ~isfinite(v); s = 1; elseif v > 100; s = 1e-3; else; s = 1; end
    end
end

function y = applyAxis(x, map, sgn), y = x(:, map) .* sgn; end

function y = fillGaps(x, t)
    y = x;
    for k = 1:size(x,2)
        v = x(:,k); ok = isfinite(v);
        if nnz(ok) >= 2 && nnz(~ok) > 0
            y(~ok,k) = interp1(t(ok), v(ok), t(~ok), 'linear', 'extrap');
        end
    end
end

function y = smoothProxy(x, t, tol)
    y = x; hasCFT = exist('spaps','file') == 2;
    for k = 1:size(x,2)
        v = x(:,k); ok = isfinite(v);
        if nnz(ok) < 5, continue; end
        if hasCFT, sp = spaps(t(ok), v(ok), tol*nnz(ok)); y(:,k) = fnval(sp, t);
        else, y(:,k) = interp1(t(ok), v(ok), t, 'spline'); end
    end
end

function y = interpToForce(x, tm, tf)
    y = nan(numel(tf), size(x,2));
    for k = 1:size(x,2)
        ok = isfinite(x(:,k));
        if nnz(ok) >= 2, y(:,k) = interp1(tm(ok), x(ok,k), tf, 'spline', 'extrap'); end
    end
end

function s = table2struct2(row)
    if istable(row); s = table2struct(row(1,:)); else; s = row; end
end

function v = getField(s,n,d)
    v = d;
    if ~(isstruct(s) && isfield(s,n)), return; end
    x = s.(n);
    if isempty(x), return; end
    if (isstring(x)||ischar(x)) && all(ismissing(string(x))), return; end   % CSV empty -> <missing>
    if isnumeric(x) && all(isnan(x(:))), return; end
    v = x;
end

function v = firstFinite(x, dflt)
    x = double(x);
    if isscalar(x) && isfinite(x); v = x; else; v = dflt; end
end

function v = getOpt(o,n,d)
    if isstruct(o) && isfield(o,n) && ~isempty(o.(n)); v = o.(n); else; v = d; end
end
