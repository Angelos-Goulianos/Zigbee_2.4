
% This is the initialization script for Zigbee hf end to end transmitter
% receiver.When shape_method is set to 'shaped', the algorithm implements
% an optimum matched filtering half-sinusiodal receiver after the
% demodulator. When shape_method is set to notshaped a sub-optimum
% receiver is implemented which compares the received signal to the OQPSK
% reference constellation, according to MATLAB comm. toolbox objects.
% -Number of bytes should not exceed 127 according to the Zigbee standard
% -SNR values should be given in dB
% -example to run :[PER,BER]=zigbee_hf_phy(100,100,[-20:5:10],'shaped');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Angelos Goulianos%
% University of Bristol/Sphere project%
% January 2014%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PER,BER]=zigbee_hf_phy(No_packets,No_bytes,SNR,shape_method);

if No_bytes>127
    error('The maximum number of bytes allowwed is 127')
end

switch shape_method
    case 'shaped'
        zigbee_hf_shape % Zigbee hf (high frequency) (2.4 GHz) shaped
    case 'notshaped'
        zigbee_hf_noshape % Zigbee hf not shaped
end

