function h = maxhorsign(signrestriction)
%MAXHORSIGN Longest horizon used in a set of sign restrictions.
%
%   h = MAXHORSIGN(signrestriction)
%
%   Input
%     signrestriction  cell array of strings 'y(a,b,c)...' (see
%                      IRESPONSE_SIGN); b, the horizon, is read between the
%                      first and the second comma
%
%   Output
%     h                scalar, the largest horizon appearing in any b
%
%   Used to compute the responses only up to the horizon that the
%   restrictions need while searching for a valid rotation.
%
%   See also IRESPONSE_SIGN, IRESPONSE_SIGN_NARRATIVE, IRESPONSE_SIGN_HMOMENTS.

for ii = 1 : size(signrestriction,2)
    tmp    = (signrestriction{ii});
    index  = strfind(tmp,','); 
    eval(['tmp1 = ' tmp( index(1)+1 : index(2)-1 ) ';']);
    hh(ii) = max(tmp1);
end

h = max(hh);
