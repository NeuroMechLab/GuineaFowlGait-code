function [T, factors] = scaleByDimension(T, dimMap, M, L, varargin)
%SCALEBYDIMENSION  Nondimensionalize table columns by DIMENSION (dynamic similarity).
%
%   [T, factors] = scaleByDimension(T, dimMap, M, L)
%   [T, factors] = scaleByDimension(T, dimMap, M, L, 'g', 9.81, 'suffix', '_n')
%
%   Dynamic-similarity normalization (Alexander & Jayes 1983) using base
%   quantities mass M, characteristic (leg) length L, and gravity g. Any
%   mechanical quantity has dimensions M^a * L^b * T^c; since T ~ sqrt(L/g) the
%   DIVIDING factor is  M^a * L^(b + c/2) * g^(-c/2)  and the dimensionless value
%   is raw/factor. Columns that share a DIMENSION share the factor automatically
%   (e.g. hip moment and ankle moment, both Force*Length) — you declare the
%   dimension once, never a per-column formula.
%
%   Scaling conventions follow Table 2 of Birn-Jeffery & Daley 2018
%   (DOI 10.1242/jeb.152538, as cited by the user/co-author). The derived
%   factors were checked against that table's values during the skill build.
%
%   INPUTS
%     T       - table with the raw columns.
%     dimMap  - containers.Map OR an N-by-2 cellstr {column, dimension; ...}.
%               dimension is a name (see DIMENSIONS below) or a 1x3 [a b c]
%               (M,L,T) exponent vector.
%     M, L    - body mass and characteristic length. Each may be a scalar, a
%               numeric column vector (per-row), or the NAME of a column in T.
%   NAME-VALUE
%     'g'      - gravity (default 9.81).
%     'suffix' - appended to new column names (default '_n').
%
%   OUTPUTS
%     T       - input table with added dimensionless columns.
%     factors - struct: factors.(column) = numeric factor(s) used.
%
%   Example
%     dimMap = {'peakForce','force'; 'Ank_PkMoment','moment'; ...
%               'Hip_PkMoment','moment'; 'kleg','stiffness'; 'stanceT','time'};
%     T = scaleByDimension(T, dimMap, 'bodyMass_kg', 'legLen_m');

p = inputParser;
p.addParameter('g', 9.81, @isnumeric);
p.addParameter('suffix', '_n', @(s)ischar(s)||isstring(s));
p.parse(varargin{:});
g = p.Results.g; suffix = char(p.Results.suffix);

% ---- resolve M and L to numeric column vectors -------------------------
Mv = local_resolve(T, M);
Lv = local_resolve(T, L);

% ---- normalize dimMap to an N-by-2 cell --------------------------------
if isa(dimMap, 'containers.Map')
    k = keys(dimMap); v = values(dimMap);
    dimCell = [k(:), v(:)];
else
    dimCell = dimMap;
end

factors = struct();
for i = 1:size(dimCell, 1)
    col = dimCell{i, 1};
    dim = dimCell{i, 2};
    if ~ismember(col, T.Properties.VariableNames)
        error('scaleByDimension:missingCol', 'column "%s" not in table', col);
    end
    e = local_factorExps(dim);           % [eM eL eg]
    factor = (Mv .^ e(1)) .* (Lv .^ e(2)) .* (g .^ e(3));
    T.([col suffix]) = T.(col) ./ factor;
    factors.(col) = factor;
end
end

% =========================================================================
function v = local_resolve(T, x)
if ischar(x) || (isstring(x) && isscalar(x))
    v = T.(char(x));
elseif isscalar(x)
    v = repmat(x, height(T), 1);
else
    v = x(:);
end
end

% =========================================================================
function e = local_factorExps(dim)
%LOCAL_FACTOREXPS  -> [eM eL eg] exponents of the dividing factor.
if isnumeric(dim)                 % raw [a b c] in M,L,T
    a = dim(1); b = dim(2); c = dim(3);
else
    [a, b, c] = local_dimension(char(dim));
end
e = [a, b + c/2, -c/2];
end

% =========================================================================
function [a, b, c] = local_dimension(name)
%LOCAL_DIMENSION  Named quantity -> (M,L,T) exponents. Mirror of kernel.py.
switch lower(name)
    case {'dimensionless','angle'},              d = [0 0 0];
    case 'mass',                                 d = [1 0 0];
    case 'length',                               d = [0 1 0];
    case 'area',                                 d = [0 2 0];
    case 'volume',                               d = [0 3 0];
    case 'time',                                 d = [0 0 1];
    case {'frequency','freq','angular_velocity'},d = [0 0 -1];
    case {'velocity','speed'},                   d = [0 1 -1];
    case 'acceleration',                         d = [0 1 -2];
    case 'force',                                d = [1 1 -2];
    case {'force_rate','load_rate','loadrate'},  d = [1 1 -3];
    case {'impulse','momentum'},                 d = [1 1 -1];
    case {'moment','torque','rot_stiffness'},    d = [1 2 -2];
    case {'work','energy'},                      d = [1 2 -2];
    case 'power',                                d = [1 2 -3];
    case {'stiffness','leg_stiffness'},          d = [1 0 -2];
    case 'pressure',                             d = [1 -1 -2];
    otherwise
        error('scaleByDimension:unknownDim', ...
              ['unknown dimension "%s". Use a known name or a raw ' ...
               '[a b c] (M,L,T) exponent vector.'], name);
end
a = d(1); b = d(2); c = d(3);
end
