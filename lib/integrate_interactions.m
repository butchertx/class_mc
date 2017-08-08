function [ rs_avg_int ] = integrate_interactions( interaction_file )
%INTEGRATE_INTERACTIONS give a real-space vector of time-averaged
%interactions
%   "interaction_file" is the filename of a csv defining interactions in
%   one spacial and one imaginary time dimension.  This will read those
%   interactions and integrate them along the imaginary time direction to
%   give a real-space, 0-frequency interaction strength
raw_data = csvread(interaction_file);
rs_avg_int = sum(raw_data, 2);


end

