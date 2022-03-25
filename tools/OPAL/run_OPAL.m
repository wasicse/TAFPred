function run_OPAL(sequence,pssm)
% run_OPAL.m  - call MoRFchibi and PROMISS, and combine scores
% Created on: 01/04/2017
% Author: Ronesh Sharma

%seq1 = fastaread('input.fasta'); %load seq file 
%seq =seq1.Sequence;

fid = fopen('input.fasta', 'r');
tline = fgetl(fid);
seq1.Header= tline;
seq1.Sequence =fgetl(fid);
seq= seq1.Sequence;
fclose(fid);

%run MoRFchibi
'processing MoRFchibi'
unix('bash run_MoRFchibi.sh'); % run MoRFchibi to get scores
path3= importdata('tr4.txt'); % get path for MoRFchibi
pathMoRFchibi= path3{1,1};
%addpath(pathMoRFchibi);
%morfchibi1 = importdata('output.txt'); %MoRFchibi scores
morfchibi1 = importdata(['' pathMoRFchibi '/output1.txt']); %MoRFchibi scores
%morfchibi_scores =morfchibi1.data(:,1);
morfchibi_scores =morfchibi1.data(:,3);

      
%run PROMISS
'processing PROMIS'
path4= importdata('tr5.txt'); % get path for PROMIS
pathPROMISS= path4{1,1};
addpath(pathPROMISS);
cd PROMIS
'where'
if nargin<2
PROMIS1=run_PROMIS(sequence); %if sequence available only
else
PROMIS1=run_PROMIS(sequence,pssm); %if pssm file available
end

%PROMIS1 = importdata('Scores_PROMIS.txt'); %promiss score
%PROMIS_scores =PROMIS1.data(:,1);
PROMIS_scores =PROMIS1(:,1);

MoRFchi_p = Process_score(morfchibi_scores,4); %Process MoRFchibi score
'Combining MoRFchibi and PROMIS'
OPAL_ = (MoRFchi_p+ PROMIS_scores)/2; %combine scores for MoRFchibi and PROMIS
OPAL = Process_score(OPAL_,12); %Process score
cd ..
fileID = fopen('Scores_OPAL.txt','w');  % save scores in txt file scores_OPAL.txt
fprintf(fileID,'%3s  %6s  %6s %6s  %6s \n','No:', 'residues','OPAL','PROMISS','MoRFchibi');
for rry=1:size(PROMIS_scores,1)
fprintf(fileID,'%0.1f  %3s  %f %f %f \n',rry,seq(rry),OPAL(rry,1),PROMIS_scores(rry,1),MoRFchi_p(rry,1) );
end
fclose(fileID);
'OPAL scores saved'
%remove file 
unix('bash removefiles_1.sh'); 
end
