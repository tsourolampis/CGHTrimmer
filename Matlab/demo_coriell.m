function time=demo_coriell(lambda,filename)

% *************************************************************************
% * FUNCTION demo_coriell
% * This demo shows the performance of DiscretizeCGH function on the  
% * standard Coriell Dataset, widely considered to be a golden dataset.  
% * It takes as input a value for the regularization parameter
% * lambda and a file name which will automatically be created to APPEND
% * the output. The file runs CGHtrimmer on every chromosome of every cell 
% * line of the Coriell. 
% * INPUT
% * lambda: the regularization parameter 
% * filename: output file 
% * OUTPUT 
% *time: a matrix #cell lines x # chromosomes containing the running times
% *************************************************************************

% *****************************************
% * Copyright (c) Charalampos Tsourakakis *
% *****************************************


%% input check 
if nargin==0
    lambda = 0.2;
    filename='CoriellAnalysis4Report.analysis';
end
if nargin ==1
    filename='CoriellAnalysis4Report.analysis';
end


fid = fopen(filename,'a+');

%% load input 
load coriell_baccgh;
data=coriell_data;
fish = data.FISHMap;
genomicposition = data.GenomicPosition;
clear coriell_data

n = length(data.Sample);

time = zeros(n,22);

%% trimmer
for sample = 1:n
    for chr = 1:22

        fprintf(fid,'*******************************************************\n');
        fish = data.FISHMap;
        Chromosome=data.Chromosome;
        chrindex = find(Chromosome==chr);
        data2= data.Log2Ratio(chrindex,sample);
        [dataclean indices]=  removeNaN(data2);
        indices = chrindex(indices);
        tic
        OPT = DiscretizeCGH(dataclean, lambda,  'log',2,genomicposition);
        time(sample,chr)=toc;
        fprintf(fid,'Cell line %s \n', char(data.Sample(sample)));
        fprintf(fid,'Chromosome id %d \n',chr);
        fprintf(fid, 'discretization took %f seconds for %s and Chromosome %d\n',time(sample,chr),char(data.Sample(sample)),chr);
        fprintf(fid,'=====================================================\n');
        fprintf(fid,'                       GAIN                          \n');
        fprintf(fid,'=====================================================\n');
        trisomies= find(OPT(:,5)>.3);
        tmp=indices(trisomies);
        for i = 1 : length(trisomies)
            fprintf(fid,'Chromosome position : %s  Genomic position %d Fitted Value (log) %f Variance of segment %f \n',fish{tmp(i)},genomicposition(tmp(i)),OPT(trisomies(i),5),  OPT(trisomies(i),6) );
        end
        fprintf(fid,'=====================================================\n');
        fprintf(fid,'                      LOSS                           \n');
        fprintf(fid,'=====================================================\n');
        monosomies=find(OPT(:,5)<-0.3);%<-- typical thresholding values in papers
        tmp=indices(monosomies);
        for i = 1 : length(monosomies)
              fprintf(fid,'Chromosome position : %s   Genomic position %d Fitted Value (log) %f Variance of segment %f \n',fish{tmp(i)},genomicposition(tmp(i)),OPT(monosomies(i),5),  OPT(monosomies(i),6) );
        end
        fprintf(fid,'*******************************************************\n');
    
    end

end
