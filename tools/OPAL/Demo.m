 
% Demo.m - demo to run OPAL and PROMIS predictors
% Created on: 01/04/2017
% Author: Ronesh Sharma
%

clear all;
close all;

%PROMIS
%cd PROMIS
%run_PROMIS('sequence'); % if only query sequence given
%run_PROMIS('sequence','pssm'); % if pssm available for query seq
%cd ..


%OPAL 
run_OPAL('sequence'); % if only query sequence given
%run_OPAL('sequence','pssm'); % if pssm available for query seq
