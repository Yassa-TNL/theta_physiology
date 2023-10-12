function [csc_ts, csc, sFreq] = ReadCSC_ag(fname, Trange, output_sFreq)
%
%    function [csc_ts, csc] = ReadCSC_ag(fname, Trange, out_sFreq)
%
%    note: the main difference from 'ReadCR.m' is that the time is in micros.
%    This function allows for downsampling. This is based on
%    ReadCR_tsd, but is not structure based.
%
%    Read in continuously sampled data (csc file). 
%    Inputs:
%       fname: obvious
%       Trange: (optional) 2 element ts vec for reading file
%             [T_start,T_end] in microseconds (not 0.1ms as in others)!!!!!
%             if empty, then full record is returned.
%       output_sFreq: (optional) desired sample rate in samples/sec (1000 = 1kHz);
%             this must be less than original data. 
%
%     Outputs:
%       csc_rt: timestamps of data in ms (not 0.1ms)
%       csc:    sampled data in microvolts
%       sFreq:  sample rate in Hz
%
%
%  N.B. When output frequency is specified, the number of samples will
%  always be determined by the sample time range*output frequency.
%  Otherwise, small time ranges (WRST the sampling frequency) can cause the
%  return of different length vectors across epocs.
%
%     Aaron Gruber 2010_8_13
%
%       2012_2_7 - AG revised to improve resampling and return vector of
%       consistent length (as determined by desired sampling frequency and
%       the time interval). Also now can properly resample data that
%       contains gaps (recording was paused).

%  TO DO:
%           - add compatibility for OSx/UNIX

% read the data
    if nargin == 1 
        [ts, csc] = Nlx2MatCSC( fname , [1,0,0,0,1] , 0, 1 );     % read all
    elseif nargin == 2 || nargin == 3
        if(isempty(Trange))
            [ts, csc] = Nlx2MatCSC( fname , [1,0,0,0,1] , 0, 1 ); % read all if no time range specified (e.g. case of re-sampling)
        else
            [tsFirst] = Nlx2MatCSC( fname , [1,0,0,0,0] , 0, 3, [0 1]); % read first 2 data points to determine sampling rate of blocks 
            [ts, csc] = Nlx2MatCSC( fname , [1,0,0,0,1] , 0, 4, [Trange(1)-diff(tsFirst), Trange(2)+diff(tsFirst)] ); % read a range; start with block that was written when beginning of time interval occured
            if(Trange(1)<ts(1) || Trange(2)>ts(end))
                warning(['Trange is out of bounds for file: ', fname])
            end
        end
    elseif nargin > 3
        error('Too many inputs')
    end

% data are stored in blocks of 512 samples - reshape to be 1xn
    nBlocks = size(csc,1); 
    sFreq = nBlocks*1000000/median(diff(ts));  % compute sample rate in data
    csc=single(reshape(csc,1,length(csc(:))));  
    dT = 1000000/sFreq; % in tstamps
    blockSize = 512;    
    csc_ts = zeros(size(csc));

    for i = 1:nBlocks % -- assign Time Stamp for each bin; have to do it this way because recording may have time gaps
      csc_ts(i:nBlocks:end) = ts + dT * (i-1); 
    end
    clear ts;         % get rid of block ts just to avoid mistakes later

% Restrict the range with more resolution than by block (as given by Nlx2MatCSC).
% If no time range given, but re-sampling is specified, this sets the
% Trange as needed later.
    if nargin>1        
        if(isempty(Trange))
            Trange = [csc_ts(1), csc_ts(end)]; % assign here in case re-sampling is desired (the next block of code)
        else
            indx = find(csc_ts>=Trange(1) & csc_ts<(Trange(2)+dT)); % upper limit is the sample immediately following uppper range (this is so interpolation will not attempt to get a number > data chunk, which returns NaN)
            csc = csc(indx);
            csc_ts = csc_ts(indx);
        end
    end    

% Resample data so as to start on Trange(1) and occur with interval
% 1/output_sFreq. The last point will be the sample in the new time
% basis closest to Trange(2). This is done in a loop w/ concatonation in case 
% there are gaps in the recording (aquisition was paused), which throws off
% the relationship between the nuber of samples for a time range spanning
% the gap.
    if nargin==3
        Tjump_idx = find(diff(csc_ts)>median(diff(csc_ts))*4);  % find indicies of any gaps in recording (csc_ts is now restricted to Trange)
        Tjump_idx = [0, Tjump_idx, numel(csc_ts)];          % time epochs (defined by gaps in recording) 
        out_csc_ts = [];
        out_csc = [];
        for k=2:numel(Tjump_idx)                            % loop over epochs
            st = csc_ts(Tjump_idx(k-1)+1);                  % set start and end times for epoch
            et = csc_ts(Tjump_idx(k));
            n_out = floor(((et-st)/1000000)*output_sFreq);  % determine the number of samples; floor makes sure it is within the limits of the data.
            etNew = st+n_out*1000000/output_sFreq;          % compute new end time - this is the re-sample time closest to the desired Trange 

            part_ts = linspace(st,etNew,n_out);             % the new time basis
            out_csc_ts =  [out_csc_ts, part_ts];
            out_csc = [out_csc, interp1(csc_ts, csc, part_ts)];  % concatonate interpolated points starting at beginning of each epoch            
        end
        csc=out_csc;
        csc_ts = out_csc_ts;
        sFreq = output_sFreq; 
    end

% set the proper scale for the signal (based on digitizer scale)
    header = Nlx2MatCSC( fname , [0,0,0,0,0] , 1, 3, [0]);  % grab headder
    ADC_char=header{~cellfun(@isempty, strfind(header, '-ADBitVolts'))}; % get line with hash tag
    bitRes = str2num(ADC_char(strfind(ADC_char, ' '):end)); % read data after whitespace and cast to numeric
    csc = csc.*bitRes.*1e6;  

