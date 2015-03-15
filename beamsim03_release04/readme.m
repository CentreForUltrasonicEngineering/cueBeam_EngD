% readme.m
% release 4 notes:
%
% *-> i have finally been able to compile these for 64-bit windows, 64-bit
%     matlab
%
% *-> fixed a bug causing the first line of the calculated image to be
%     overwritten by last line. verification now confirms -120.4dB accuracy
%     even with "fast math" forced enabled
%
% *-> changed block size back to <16,24> (optimised for 9xxx, 2xx cards) 
%     to optimize for Fermi cards change block size to <16,26>
%     may need to change to <16,16> or lower to run on 8xxx series cards, 
%     not tested.
%
% note: official project name is "cuBEAM", i'll rename the files later

% run benchmark
testbed
% run live demo
rtdemo

% note: b
% beamsim03_* generation assume that array is in fluid and has pointlike
% elements (ominidirectional).
% into-solid directivity and material interfaces to come in upcoming
% versions, subject to funding availability

