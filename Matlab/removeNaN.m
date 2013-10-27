function [data ind] = removeNaN(data)

% This function takes as input a table, recognizes for every column the
% rows that contain a NaN and then removes all those lines! 
% This turns out to be a typical task when working with CGH data. 
% INPUT
% data: a 2d table
% OUTPUT
% data: without any nans

% ************************************************************************
% TODO: none
% DEPENDENCIES: none 
% AUTHOR: Charalampos Tsourakakis (ctsourak@cs.cmu.edu)
% DATE: Wednesday 18 November 2009
% ************************************************************************

[rows cols] = size(data); 
toberemoved=[];
for i = 1 : cols
    toberemoved = [toberemoved;  find(isnan(data(:,i))==1) ];
end

ind = setdiff(1:rows, toberemoved);
data(unique(toberemoved),:)=[];
