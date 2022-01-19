% RELAX_filtbutter() - filters the data using a zero-phase butterworth
%                   filter. Either a band pass, band stop, high pass, or low pass 
%                   filter can be implemented. The filter order is defined by the user.
%                   This function uses the matlab butter and filtfilt functions.
%
% Usage:
%   >>  EEG = RELAX_filtbutter( EEG, high, low, ord, type );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   high            - integer (required). Highpass frequency (allow frequencies
%                   above this value with bandpass or high pass)
%                   Example: 1
%   low             - integer (required). Lowpass frequency (allow
%                   frequencies below this value with bandpass or low pass)
%                   Example: 50
%   ord             - integer (required). Filter order.
%                   Example: 4 (designs a fourth order butterworth filter)
%   type            - 'str' (required). 'bandpass' | 'bandstop' | 'highpass' | 'lowpass'. Designs either
%                   a zero-phase band pass, band stop, high pass, or low pass butterworth filter
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
%
% Examples
%   EEG = RELAX_filtbutter( EEG, 1, 100, 4, 'bandpass' ); %zero-phase, 4th-order bandpass butterworth filter between 1-100 Hz
%   EEG = RELAX_filtbutter( EEG, 48, 52, 4, 'bandstop' ); %zero-phase, 4th-order bandstop butterworth filter between 48-52 Hz
%   EEG = RELAX_filtbutter( EEG, 1, [], 4, 'highpass' ); %zero-phase, 4th-order high pass butterworth filter allowing frequencies above 1 Hz
%   EEG = RELAX_filtbutter( EEG, [], 45, 4, 'lowpass' ); %zero-phase, 4th-order low pass butterworth filter allowing frequencies below 45 Hz
% 
% See also:
%   eegfiltnew

% Copyright (C) 2016  Nigel Rogasch, Monash University,
% nigel.rogasch@monash.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% Change log
% 20.6.2018 - included low and high pass filter options
% 20.6.2018 - changed to filtering across concatenated trials instead of
% individual trials to try and minimise edge artifacts and speed up script

function EEG = RELAX_filtbutter( EEG, high, low, ord, type )

if nargin < 5
	error('Not enough input arguments.');
end

%Check inputs
if ~isnumeric(high)
    error('Input ''high'' needs to be a number, not a string. e.g. 1, not ''1''.')
elseif ~isnumeric(low)
    error('Input ''low'' needs to be a number, not a string. e.g. 100, not ''100''.')
elseif high > low
    error('Input ''high'' needs to be less than input ''low''.')
elseif ~isnumeric(ord)
    error('Input ''ord'' needs to be a number, not a string. e.g. 4, not ''4''.')
elseif ~(strcmp(type,'bandpass') || strcmp(type,'bandstop') || strcmp(type,'highpass') || strcmp(type,'lowpass'))
    error('Input ''type'' needs to be either ''bandpass'', ''bandstop'', ''highpass'', or ''lowpass''.')
end

if strcmp(type,'highpass') && ~isempty(low)
    warning('Performing high pass filter - ignoring low pass input value.');
end
if strcmp(type,'lowpass') && ~isempty(high)
    warning('Performing low pass filter - ignoring high pass input value.');
end

Fs = EEG.srate;
ordIn = ord/2;

if strcmpi(type,'bandstop')
    type = 'stop';
end
               
if strcmp(type,'bandpass') || strcmp(type,'stop')
    [z1 p1] = butter(ordIn, [high low]./(Fs/2), type);
elseif strcmp(type,'highpass')
    [z1 p1] = butter(ordIn, high./(Fs/2), 'high');
elseif strcmp(type,'lowpass')
    [z1 p1] = butter(ordIn, low./(Fs/2), 'low');
end

data = double(EEG.data);
temp = NaN(size(data,1),EEG.pnts,size(data,3));
for x = 1:size(data,1) 
    dataIn = data(x,:);
    dataFilt1 = filtfilt(z1,p1,dataIn);
    temp(x,:,:) = reshape(dataFilt1,1,EEG.pnts,[]);
end 

EEG.data = temp;

%display message
if strcmp(type,'bandpass') || strcmp(type,'stop')
    fprintf('Data filtered using a %s zero-phase butterworth filter (order = %d) between %d and %d Hz.\n',type,ord,high,low);
elseif strcmp(type,'highpass')
    fprintf('Data filtered using a high pass zero-phase butterworth filter (order = %d) below %d Hz.\n',ord,high);
elseif strcmp(type,'lowpass')
    fprintf('Data filtered using a low pass zero-phase butterworth filter (order = %d) above %d Hz.\n',ord,low);
end
    

end
