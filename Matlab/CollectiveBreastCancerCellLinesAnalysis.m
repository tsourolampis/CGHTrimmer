function [trimmertime cbstime]= CollectiveBreastCancerCellLinesAnalysis(BreastCancerCellLinesBerkeley,cols,chrtotest,lambda,trimmerfilename,cbsfilename)

% FUNCTION CollectiveBreastCancerCellLinesAnalysis
% INPUT
% BreastCancerCellLinesBerkeley: matrix containing the Berkeley Breast
% Cancer Cell Lines
% cols: column ids corresponding to the cell lines to be analyzed
% chrtotest: chromosomes to be analyzed for every cell line 
% lambda: regularization parameter 
% trimmerfilename: file name for the output of CGHTRimmer 
% cbsfilename: file name for the output of 
% OUTPUT 
% trimmertime,cbstime: timings
% Note: for the CGHseg algorithm, please download the code from Picard's
% web site:
% http://www.agroparistech.fr/mia/doku.php?id=productions:logiciels 

% *****************************************
% * Copyright (c) Charalampos Tsourakakis *
% *****************************************


if nargin ~= 6
   error(' Usage: [trimmertime cbstime]= CollectiveBreastCancerCellLinesAnalysis(BreastCancerCellLinesBerkeley,cols, chrtotest,lambda,trimmerfilename,cbsfilename) \n');
end




textdata=BreastCancerCellLinesBerkeley.textdata;
data=BreastCancerCellLinesBerkeley.data;
clear BreastCancerCellLinesBerkeley 

trimmertime = zeros(length(cols),length(chrtotest));
cbstime     = zeros(length(cols),length(chrtotest));
fid1 = fopen(trimmerfilename,'a+');
fid2 = fopen(cbsfilename,'a+');


for c=1:length(cols) %for each cell line 
    fprintf(' Analyzing cell line %s \n',char(textdata(1,cols(c)+1)));
    close all %<-- save memory for plotting, if you want to see all plots together comment this line
    for j = 1 : length(chrtotest) %<-- for every chromosome
        chromosome = chrtotest(j);
        indices = find(data(:,1)==chromosome); %<-- chromosome indices
        tmp = data(indices,cols(c));
        [dataclean ind]= removeNaN(tmp);
        indices = indices(ind);
        genomicposition  = data(indices,2);%<-- genomic positions
        dataclean = data(indices,cols(c));
        tic
        OPT  = DiscretizeCGH(dataclean, lambda, 'log',2);
        trimmertime(c,j)=toc; % c th row corresponds to cols(c) cell line
        fprintf(fid1,'*******************************************************\n');
        fprintf(fid1,'Cell line %s \n', char(textdata(1,cols(c)+1)));
        fprintf(fid1,'Chromosome id %d \n',chromosome);
        fprintf(fid1, 'discretization took %f seconds for %s and Chromosome %d\n',trimmertime(c,j), char(textdata(1,cols(c)+1)),chromosome);
        fprintf( 'CGHTRIMMER: discretization took %f seconds for %s and Chromosome %d\n',trimmertime(c,j), char(textdata(1,cols(c)+1)),chromosome);
        fprintf(fid1,'=====================================================\n');
        fprintf(fid1,'                        GAIN                         \n');
        fprintf(fid1,'=====================================================\n');
        trisomies= find(OPT(:,5)>.3);%<-- typical thresholding values in papers, see reference in our paper
        for i = 1 : length(trisomies)
            fprintf(fid1,'Value %f Genomic position %d Fitted Value (log) %f Variance of segment %f \n',dataclean(trisomies(i)), genomicposition( trisomies(i)),OPT(trisomies(i),5),  OPT(trisomies(i),6) );
            %fprintf('Trimmer: Value %f Genomic position %d Fitted Value (log) %f Variance of segment %f \n',dataclean(trisomies(i)),genomicposition( trisomies(i)),OPT(trisomies(i),5),  OPT(trisomies(i),6) );
            %<-- uncomment above command to see the segmentation info on matlab display screen
        end
        
        fprintf(fid1,'=====================================================\n');
        fprintf(fid1,'                        LOSS                         \n');
        fprintf(fid1,'=====================================================\n');
        monosomies=find(OPT(:,5)<-0.3);%<-- typical thresholding values in papers
        for i = 1 : length(monosomies)
            fprintf(fid1,'Value %f  Genomic position %f Fitted Value (log) %f Variance of segment %f \n',dataclean(monosomies(i)),genomicposition(monosomies(i)),OPT(monosomies(i),5),  OPT(monosomies(i),6) );
            %fprintf('Trimmer: Value %f Genomic position %f Fitted Value (log) %f Variance of segment %f \n',dataclean(monosomies(i)),genomicposition(monosomies(i)),OPT(monosomies(i),5),  OPT(monosomies(i),6) );
            %<-- uncomment above command to see the segmentation info on matlab display screen
        end
        fprintf(fid1,'*******************************************************\n');
        %plot uncomment/comment to not plot/plot
        K = max(OPT(:,4));
        figure
        hold on
        d = OPT(:,1); %<-- column 1 contains the points
        tmp1 = max(d(:))+1; %<-- define the largest and smallest placeholders for y axis
        tmp2 = min(d(:))-1;
        
        
        for i = 1 : K
            ind = find( OPT(:,4) == i );
            st = min(ind);
            en = max(ind);
            %line([genomicposition(en) genomicposition(en)],[tmp1 tmp2],'color','r','MarkerSize',6,'Tag','Breakpoint');
            plot(genomicposition(ind),repmat( mean(OPT(ind,1)), 1 ,length(st:1:en)) ,'color','red','LineWidth',4,'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);
        end
        plot(genomicposition, OPT(:,1),'.b');
        ylims = [min(min(dataclean), -1), max(max(dataclean), 1)];
        ylim(gca, ylims)
        xlabel('Genomic Position','fontsize',20,'fontweight','b')
        ylabel('Log2Ratio','fontsize',20,'fontweight','b')
        title(strcat(' Chromosome ', int2str(chromosome)),'fontsize',30,'fontweight','b') ;
        hs_cytobands = cytobandread('hs_cytoBand.txt');
        chromosomeplot(hs_cytobands, chromosome, 'addtoplot', gca,  'unit', 2);
        
        %% COMPETITOR OLSHEN ET AL
        if length(genomicposition)-length(dataclean)~=0
            error('debug')
        end
        mydummystruct=[chromosome*ones(length(genomicposition),1) genomicposition dataclean];
        tic;a=cghcbs(mydummystruct,'showplot',chromosome);  %<- prefer this to show breakpoints
        cbstime(c,j)=toc;
        title(strcat(' Chromosome ', int2str(chromosome)),'fontsize',30,'fontweight','b') ;
        hs_cytobands = cytobandread('hs_cytoBand.txt');
        chromosomeplot(hs_cytobands, chromosome, 'addtoplot', gca,  'unit', 2);
        fprintf(fid2,'*******************************************************\n');
        fprintf(fid2,'Cell line %s \n', char(textdata(1,cols(c)+1)));
        fprintf(fid2,'Chromosome id %d \n',chromosome);
        fprintf(fid2, 'discretization took %f seconds for %s and Chromosome %d\n',cbstime(c,j), char(textdata(1,cols(c)+1)),chromosome);
        fprintf('CBS: discretization took %f seconds for %s and Chromosome %d\n',cbstime(c,j), char(textdata(1,cols(c)+1)),chromosome);
        fprintf(fid2,'=====================================================\n');
        fprintf(fid2,'                        GAIN                         \n');
        fprintf(fid2,'=====================================================\n');
        trisomies= find(a.SegmentData.Mean>.3);
        
        for i = 1 : length(trisomies)
            ind1 =  find( genomicposition== a.SegmentData.End(trisomies(i)));
            ind2 =  find( genomicposition== a.SegmentData.Start(trisomies(i)));
            for k = 1: (ind1(end) - ind2(1) +1)
                fprintf(fid2,'Genomic position : %f  Fitted Value (log) %f \n',genomicposition(ind2(1)+k-1),a.SegmentData.Mean(trisomies(i)));
                %fprintf('CBS: Genomic position : %f  Fitted Value (log) %f \n',genomicposition(ind2(1)+k-1),a.SegmentData.Mean(trisomies(i)));
                %<-- uncomment above command to see the segmentation info on matlab display screen
            end
        end
        
        fprintf(fid2,'=====================================================\n');
        fprintf(fid2,'                    LOSS                             \n');
        fprintf(fid2,'=====================================================\n'); 
        monosomies= find(a.SegmentData.Mean<-.3);
        for i = 1 : length(monosomies)
           % ind1 =  find( genomicposition == a.SegmentData.End(monosomies(i)));
            ind2 =  find( genomicposition== a.SegmentData.Start(monosomies(i)));
            for k = 1: find( genomicposition== a.SegmentData.End(monosomies(i))) - find( genomicposition== a.SegmentData.Start(monosomies(i))) +1
                fprintf(fid2,'CBS: Chromosome position : %f  Fitted Value (log) %f \n',genomicposition(ind2(1)+k-1),a.SegmentData.Mean(monosomies(i)));
                %fprintf('CBS: Chromosome position : %f  Fitted Value (log) %f \n',genomicposition(ind2(1)+k-1),a.SegmentData.Mean(monosomies(i)));
                %<-- uncomment above command to see the segmentation info on matlab display screen
            end
        end
        fprintf(fid2,'*******************************************************\n'); 
    end
end

fclose(fid1);
fclose(fid2);
