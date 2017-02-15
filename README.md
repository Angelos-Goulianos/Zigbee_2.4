# Zigbee_2.4

This is the initialization script for Zigbee hf end to end transmitter
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
%


