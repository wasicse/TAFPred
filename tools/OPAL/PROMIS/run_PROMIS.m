function [PROMIS1]=run_PROMIS(sequence,pssm)

% run_PROMISS  - call bigram and windows model, and combine scores
% Created on: 01/04/2017
% Author: Ronesh Sharma
%
% input - query seqeunce and (pssm file if available)
%
% output - scores saved in scores_PROMIS.txt 
'NO ?'
if nargin < 2
unix('bash run_spider2_seq.sh'); %run spider2 to get structural information given input as sequence
else
unix('bash run_spider2_pssm.sh'); %run spider2 to get structural information given input as pssm
end
path1= importdata('tr3.txt'); % get and add path for LIBSVM
pathlibsvm = path1{1,1};
addpath(pathlibsvm);
path2= importdata('tr2.txt'); % get and add path for spider2 
pathspider = path2{1,1};
%addpath(pathspider);

%seq1 = fastaread('input.txt'); %load seq file 
%seq =seq1.Sequence;
fid = fopen('input.txt', 'r');
tline = fgetl(fid);
seq1.Header= tline;
seq1.Sequence =fgetl(fid);
seq= seq1.Sequence;
fclose(fid);

HESa1 = importdata(['' pathspider '/1a1xA.hsa2']); %load hsa2 file 
sp_hsa2= HESa1.data(:,3:5);
HESb1 = importdata(['' pathspider '/1a1xA.hsb2']); %load hsb2 file 
sp_hsb2 = HESb1.data(:,3:5);
spd31 = importdata(['' pathspider '/1a1xA.spd3']); %load spd3 file 
sp_spd3= spd31.data(:,4:11);

profile_bigram_m1_m2 = [ sp_hsa2/100 sp_spd3(:,1)/100 sp_spd3(:,4)/100] ; % structural info for bigram
profile_windows_m_1 = [ sp_hsa2(:,2)/100 ] ; % structural info for windows method
ab=0;
hh=size(profile_windows_m_1,1); % check if protein length is less than total Flanking size
if size(profile_windows_m_1,1)< 42
    ab=1;
    profile_windows_m1 =[ profile_windows_m_1 ; zeros(42-hh,size(profile_windows_m_1,2))];
else
     profile_windows_m1 = profile_windows_m_1;
end
    
%call prediction function
MoRFscores_m1_bi = test_bigram_method(profile_bigram_m1_m2); 
MoRFscores_m1_wn = test_windows_method(profile_windows_m1);
if ab==1  % reset if zeros were added due to protein length less then Flanking size
  MoRFscores_wn = MoRFscores_m1_wn(1:(size(MoRFscores_m1_wn,1)-(42-hh)),:) ;
else
     MoRFscores_wn = MoRFscores_m1_wn;
end
MoRFbi_1=MoRFscores_m1_bi(:,1); % score outputs from each model
MoRFbi_2=MoRFscores_m1_bi(:,2);
MoRFwin_1=MoRFscores_wn;

%Post process each score using score conservation equation 2 from paper
MoRFbi_1p= Process_score(MoRFbi_1,20);
MoRFbi_2p= Process_score(MoRFbi_2,20);
MoRFwin_1p= Process_score(MoRFwin_1,12);

PROMIS= (MoRFbi_1p+MoRFbi_2p+MoRFwin_1p)/3 ;  %common averaging
fileID = fopen('Scores_PROMIS.txt','w');  % save scores in txt file scores_PROMIS.txt
fprintf(fileID,'%3s  %3s  %6s %3s  %3s %3s \n','No:', 'residues','PROMIS','MoRFbi-1' ,'MoRFbi-2','MoRFwin');
for rry=1:size(MoRFbi_1,1)
fprintf(fileID,'%0.1f  %3s  %f %f %f %f \n',rry,seq(rry),PROMIS(rry,1),MoRFbi_1p(rry,1),MoRFbi_2p(rry,1),MoRFwin_1p(rry,1) );
end
fclose(fileID);
%remove file 
unix('bash removefiles.sh'); 
PROMIS1=PROMIS;
end 
