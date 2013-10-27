function [OPT time] = DiscretizeCGH(data, lambda,  typeofdata,base, genomicposition)

% FUNCTION DiscretizeCGH
% INPUT
% data: a one dimensional vector, the input data
% lambda: value of regularization parameter
% typeofdata: a string, either 'lin' or 'log' specifying if we are in the
% log or linear domain
% base: if we are in the log domain this should be equal to the base of the
% log being used
% genomicposition: a vector containing the genomic positions, i.e., data(1)
% is the observed value in probe/genomicposition(1) etc. This is assumed to be
% ordered in a natural way
% For a demo, run timings=demo_coriell(0.2,'demo_output.txt');

% *****************************************
% * Copyright (c) Charalampos Tsourakakis *
% *****************************************

%% check input arguments
[n m]= size(data);

if m==1
    data=data'; %<-- data is a row vector
end

if nargin == 6
    if length(genomicposition)~=length(data)
        error('data and genomicposition vector must have equal lengths \n')
    end
end


%% if the data is in the log domain
if strcmp(typeofdata,'log')
    data = base.^data; %<-- exponentiate the data...
end
%% computing all pairs of sum squared errors

tic
A = zeros(n);
A(1:n+1:n^2) = data;
for i = 1 : n
    for j = i+1:n
        A(i,j) = (j-i)/(j-i+1)*A(i,j-1)+ 1/(j-i+1)*data(j);
    end
end
e = zeros(n);
for i = 1 : n
    for j = i+1:n
        e(i,j) = e(i,j-1)+ (j-i)/(j-i+1)*(data(j)-A(i,j-1))^2;
    end
end


OPT = zeros(n,5);
OPT(:,1)=data; %<-- first column contains the data

for k = 1 : n
    breakatjerror = zeros(k,1);
    for j = 1 : k %<-- %find where the minimum occurs
        if( j == 1)
            breakatjerror(j)=e(j,k)+ 2*lambda;
        else
            breakatjerror(j) = OPT(j-1,3)+e(j,k)+ 2*lambda;
        end
        
    end
    %find where the minimum occurs, for ties we favor the index with the
    %longest line
    breakind = find( breakatjerror==min(breakatjerror));
    breakind = breakind(length(breakind));
    OPT(k,2) = breakind; %<breaking indices
    OPT(k,3) = breakatjerror(breakind); %<-- third column contains the errors
    %end
end


%% RECURSE BACK and FIND the segments along with the fitted values

if strcmp(typeofdata,'log') %<-- we go back to the log domain    
    mytmp = n;
    mycounter=1; 
    while(true)
        OPT( mytmp:-1:OPT(mytmp,2),4) = mycounter;
        mycounter=mycounter+1;
        OPT( mytmp:-1:OPT(mytmp,2),5) = log( mean(OPT(mytmp:-1:OPT(mytmp,2),1)) )/log(base);
        mytmp =  OPT(mytmp,2) - 1;
        if ~mytmp
            break
        end
    end
    OPT(:,1) = log(data)/log(base);
    mytmp = n;
    while(true)
        OPT(  mytmp:-1:OPT(mytmp,2), 6) = var( OPT( mytmp:-1:OPT(mytmp,2), 1)); %<-- variances
        mytmp =  OPT(mytmp,2) - 1;
        if ~mytmp
            break
        end
    end
    
elseif strcmp(typeofdata,'lin')
    
    mytmp = n;
    mycounter=1;   
    while(true)
        OPT( mytmp:-1:OPT(mytmp,2),4) = mycounter;
        mycounter=mycounter+1;
        OPT( mytmp:-1:OPT(mytmp,2),5) =  mean(OPT(mytmp:-1:OPT(mytmp,2),1));
        OPT(  mytmp:-1:OPT(mytmp,2), 6) = var( OPT( mytmp:-1:OPT(mytmp,2), 1)); %<-- variances
        mytmp =  OPT(mytmp,2) - 1;
        if ~mytmp
            break
        end
    end
else
    error('Type of data can be either lin or log\n');
end
time= toc;









