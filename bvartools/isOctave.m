function r = isOctave ()
%ISOCTAVE True when running in GNU Octave.
%
%   r = ISOCTAVE()
%
%   Output   r   1 in Octave, 0 in MATLAB (cached after the first call)

  persistent x;
  if (isempty (x))
    x = exist ('OCTAVE_VERSION', 'builtin');
  end
  r = x;
end