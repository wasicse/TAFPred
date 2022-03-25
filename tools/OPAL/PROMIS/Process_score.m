function Scores_P = Process_score(res_scores,hh)
%  Process_score - score conservation (post process each score)
%
% Input
% Predicted scores of query sequence
%
% Output
% MoRFscore: processed scores of query sequence
%
% Ronesh Sharma, FNU, Fiji. 
% Email: sharmaronesh@yahoo.com
% Ref. Sharma et al., 2017

we=size(res_scores,1);
%hh=12; %Flank size
new_w1=zeros(1,we+(hh*2));
new_w1(hh+1:we+hh)=res_scores';
wt=1;
for wr=hh+1:size(res_scores,1)+hh
   res_scores_P(wt,1) = (median(new_w1(wr-hh:wr+hh)) +max(new_w1(wr-hh:wr+hh)) )/2;
   wt=wt+1;
end
Scores_P  = res_scores_P;
res_scores_P=[];
end
%##############################