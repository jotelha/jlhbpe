%%
% 
%  does regular expression replacements on a set of arrays
%  use:
%   prepTerms('z_id*D_id*F/RT',3,'z_id','D_id',{'z1';'z2';'z3'},{'D1','D2','D3'})
%  
%  returns:
%   'z1*D1*F/RT'
%   'z2*D2*F/RT'
%   'z3*D3*F/RT'
%   
% 
function t = prepTerm(format,varargin)
nv = nargin-1; % number of variable arguments
eStart = 1;
eEnd = nv/2;
rStart = nv/2+1;
rEnd = nv;

expressions = varargin(eStart:eEnd);
replacements = varargin(rStart:rEnd);

for i = 1:numel(replacements)
    if ischar(replacements{i})
        replacements{i} = replacements(i);
    end
end

n = numel(replacements{1});

%prepReplace = @(varargin) cellfun( @(arg) arg{i}, r);
t = arrayfun( ...
    @(i) regexprep( ...
        format, expressions, ...
        cellfun( @(arg) arg{i}, replacements, 'UniformOutput', false ) ), ...
    (1:n)' ,...
    'UniformOutput', false);
end