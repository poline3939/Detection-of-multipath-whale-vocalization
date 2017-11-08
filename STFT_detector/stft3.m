function [d1,tt,ff1] = stft3(x, f, w, h, sr)
% D = stft(X, F, W, H)                            Short-time Fourier transform.
%	Returns some frames of short-term Fourier transform of x.  Each 
%	column of the result is one F-point fft; each successive frame is 
%	offset by H points until X is exhausted.  Data is hamm-windowed 
%	at W pts, or rectangular if W=0, or with W if it is a vector.
%	See also 'istft.m'.
% dpwe 1994may05.  Uses built-in 'fft'
% $Header: /homes/dpwe/public_html/resources/matlab/pvoc/RCS/stft.m,v 1.2 2009/01/07 04:32:42 dpwe Exp $

% x = input signal
% f = length of fft
% w = length of window, defrault hanning window
% h = length of hop step
% d = output of STFT of x, square it to get the spectrogram image

if nargin < 5;  sr = 8000; end

% expect x as a row
if size(x,1) > 1
  x = x';
end
s = length(x);

if length(w) == 1
  if w == 0
    % special case: rectangular window
    win = ones(1,f)/sqrt(f);
  else
    if rem(w, 2) == 0   % force window to be odd-len
      w = w + 1;
    end
%     halflen = (w-1)/2;
%     halff = f/2;   % midpoint of win
%     halfwin = 0.54 + 0.46*cos( pi * (0:halflen)/halflen));
%     win = zeros(1, f);
%     acthalflen = min(halff, halflen);
%     win((halff+1):(halff+acthalflen)) = halfwin(1:acthalflen);
%     win((halff+1):-1:(halff-acthalflen+1)) = halfwin(1:acthalflen);
%     win=hamming(f)';
%     win=rectwin(f)'/sqrt(f);
    win=rectwin(f)';

  end
else
  win = w;
end

w = length(win);
% win2=win./sqrt(win*win');

c = 1;

% pre-allocate output array
d = zeros((1+f/2),1+fix((s-f)/h));

d1 = zeros(f,1+fix((s-f)/h));


for b = 0:h:(s-f)
  u = win.*x((b+1):(b+f));
  t = fft(u);
  d(:,c) = t(1:(1+f/2))';
  d1(:,c) = t(1:f)';
  c = c+1;
end

tt = [0:size(d,2)]*h/sr;
% ff = [0:size(d,1)]*sr/f;
ff1 = [1:size(d1,1)]*sr/f;

% If no output arguments, plot a spectrogram
if nargout == 0
%   tt = [0:size(d,2)]*h/sr;
%   ff = [0:size(d,1)]*sr/f;
  imagesc(tt,ff1,(abs(d1).^2));
  axis('xy');
  xlabel('time / sec');
  ylabel('freq / Hz')
  % leave output variable D undefined
else
  % Otherwise, no plot, but return STFT
  D = d1;
end