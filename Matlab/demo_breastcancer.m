function demo_breastcancer

% FUNCTION demo_breastcancer
% run demo_breastcancer to see the performance of DiscretizeCGH on the
% Berkeley Breast Cancer Cell Line Database.

% *****************************************
% * Copyright (c) Charalampos Tsourakakis *
% *****************************************


fprintf('Loading Breast Cancer data .. \n')
load BreastCancerCellLinesBerkeley %<-- load data
fprintf(' Data Loaded! \n')
fprintf('Testing CGHTRIMMER and CBS on BT474, HS578T, MCF10A, T47D for Chromosomes 3 and 17 \n');
fprintf(' Default value for Lambda = 0.2 \n');

% 7 BT-474, 29 HS578T, 31 MCF10A, 52 T47D
cols = [7 29 31 52];
chrtotest=[3 17];
lambda=0.2;
% filenames
trimmerfilename='trimmer_output.txt'; cbsfilename='csb_ouput.txt';

% run algorithm
[trimmertime cbstime] = CollectiveBreastCancerCellLinesAnalysis(BreastCancerCellLinesBerkeley,cols,chrtotest,lambda,trimmerfilename,cbsfilename);
fprintf('Speedups per chromosome per cell line \n')
cbstime./trimmertime
