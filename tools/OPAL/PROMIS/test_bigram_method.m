function [MoRFscores] = test_bigram_method(profile)
% MoRF prediction algorithm (bigram  method)
%
% Input
% structural attribute of size L by n, where L is length of protein seq and n is number of attributes: 
%
% Output
% MoRFscore: scores for each residue of the query sequence
%
% Ronesh Sharma, FNU, Fiji. 
% Email: sharmaronesh@yahoo.com
% Ref. Sharma et al., 2017

load SVM_model_1_bi_1to2; % load trained models bi m1
load SVM_model_2_bi_1to1; % load trained models bi m2
  
%---------------------------------------------------------------------------
kk=1;win_flank_siz= 20;
mat3=profile ; % profile
T_len=size(mat3,1);
%AZ1=71; AZ2= 24; AZ3=23; %6:24
AZ1=83; AZ2= 30; AZ3=29;   %30:30
if T_len <AZ1 % condition if length of protein seq is less than 71 (due to 6 to 30 varying window size)
    ma_3=zeros(AZ1,size(mat3,2));
    C=round((AZ1-T_len)/2);
    A = C+1;
    B =T_len+C;
    ma_3(A:B,:)= mat3;
    mat3= ma_3;
else
    A=1;
    B=T_len;
end
  for e=A:B
  if e < AZ2 + win_flank_siz  % seq at start
%--------------------------------------------------------------------------
 %morf = 6:24;
  morf = 30:30;
      f1=1;
      for e1=1:length(morf)
          g1=morf(e1);
          if e<g1
              g2= e;
          else
              g2=g1;
          end
          for e2=1:g2 
              res_seq = mat3(e-e2+1:e-e2+1+morf(e1)-1,: );
              if e==1  
                l_f_seq=[zeros(win_flank_siz,size(mat3,2)) ]; 
              elseif  e-e2+1-1==0
                l_f_seq=[zeros(win_flank_siz,size(mat3,2)) ];  
              elseif e-e2+1-1 - win_flank_siz+1<1
                l_f_seq =  [zeros(win_flank_siz-(e-e2+1-1) ,size(mat3,2)); mat3(1:e-e2+1-1,:)];
              else
                  l_f_seq =   mat3(e-e2+1-1 - win_flank_siz+1:e-e2+1-1,:); 
              end
            if e-e2+1+morf(e1)-1+1+win_flank_siz-1 > size(mat3,1)
            r_f_seq =[mat3(e-e2+1+morf(e1)-1+1:end,:) ;zeros((e-e2+1+morf(e1)-1+1+win_flank_siz-1)-size(mat3,1),size(mat3,2))  ]; 
            else
            r_f_seq =mat3(e-e2+1+morf(e1)-1+1:e-e2+1+morf(e1)-1+1+win_flank_siz-1,:);
            end
            sample_d1_seq= [l_f_seq; res_seq ;r_f_seq  ];
            matZ=[ sample_d1_seq];
 %bigram
 %%{
%****************************
for h=1:size(matZ,2)-1 % since last attribute for m2
        for h1=1:size(matZ,1)-1
        Bi(h1,:)= matZ(h1,h) * matZ(h1+1,1:4);
        end
        Bigram(h,:)= sum(Bi,1)/size(matZ,1);
end
Bigram_1= Bigram(:)';
Bi=[];
Bigram=[];
%****************************
        for h1=1:size(matZ,1)-1
        Bi(h1,:)= matZ(h1,5) * matZ(h1+1,5);
        end
        Bigram = sum(Bi,1)/size(matZ,1);
Bigram_2= Bigram(:)';
Bi=[];
Bigram=[];
%****************************
F1(f1,:) = [Bigram_1];
F2(f1,:) = [Bigram_2];
%%}
            %--------------------
            f1=f1+1;
          end
      end
%---------------------------------------------------------------------------
    elseif e> T_len-AZ3-win_flank_siz %seq at end
       %  morf = 6:24;
         morf = 30:30;
         f1=1;
      for e1=1:length(morf)
          g1=morf(e1);
          if e> size(mat3,1)- g1+1
              g2= e-(size(mat3,1)- g1) ;
          else
              g2=1;
          end
          for e2=g2:morf(e1)
              res_seq = mat3(e-e2+1:e-e2+1+morf(e1)-1,: );
              l_f_seq =mat3(e-e2+1-win_flank_siz:e-e2+1-1,:); 
              if e==size(mat3,1) 
                r_f_seq=[zeros(win_flank_siz,size(mat3,2))];  
              elseif  e-e2+1+morf(e1)-1 ==size(mat3,1)
              r_f_seq=[zeros(win_flank_siz,size(mat3,2))];  
              elseif e-e2+1+morf(e1)-1+1+win_flank_siz-1 > size(mat3,1)
             r_f_seq = [  mat3(e-e2+1+morf(e1)-1+1:end,:) ; zeros(win_flank_siz-(size(mat3,1)-(e-e2+1+morf(e1)-1)) ,size(mat3,2)) ] ;
              else
              r_f_seq =   mat3(e-e2+1+morf(e1)-1+1:e-e2+1+morf(e1)-1+1+win_flank_siz-1,: ) ;
              end
            sample_d1_seq= [l_f_seq; res_seq ;r_f_seq ];
            matZ=[ sample_d1_seq];
 %bigram
 %%{
 %****************************
for h=1:size(matZ,2)-1 % since last attribute for m2
        for h1=1:size(matZ,1)-1
        Bi(h1,:)= matZ(h1,h) * matZ(h1+1,1:4);
        end
        Bigram(h,:)= sum(Bi,1)/size(matZ,1);
end
Bigram_1= Bigram(:)';
Bi=[];
Bigram=[];
%****************************
        for h1=1:size(matZ,1)-1
        Bi(h1,:)= matZ(h1,5) * matZ(h1+1,5);
        end
        Bigram = sum(Bi,1)/size(matZ,1);
Bigram_2= Bigram(:)';
Bi=[];
Bigram=[];
%****************************
F1(f1,:) = [Bigram_1];
F2(f1,:) = [Bigram_2];
 %%}
            %--------------------
            f1=f1+1;
          end
      end
%---------------------------------------------------------------------------------
    else % seq in middle
      %morf = 6:24;
      morf = 30:30;
      f1=1;
      for e1=1:length(morf)
          for e2=1:morf(e1)
              res_seq = mat3(e-e2+1:e-e2+1+morf(e1)-1,: ); 
              l_f_seq =mat3(e-e2+1-win_flank_siz:e-e2+1-1,:);
              r_f_seq =mat3(e-e2+1+morf(e1)-1+1:e-e2+1+morf(e1)-1+1+win_flank_siz-1,:);
            sample_d1_seq= [l_f_seq; res_seq; r_f_seq  ] ;
           matZ=[ sample_d1_seq];
 %bigram
 %%{
 %****************************
for h=1:size(matZ,2)-1 % since last attribute for m2
        for h1=1:size(matZ,1)-1
        Bi(h1,:)= matZ(h1,h) * matZ(h1+1,1:4);
        end
        Bigram(h,:)= sum(Bi,1)/size(matZ,1);
end
Bigram_1= Bigram(:)';
Bi=[];
Bigram=[];
%****************************
        for h1=1:size(matZ,1)-1
        Bi(h1,:)= matZ(h1,5) * matZ(h1+1,5);
        end
        Bigram = sum(Bi,1)/size(matZ,1);
Bigram_2= Bigram(:)';
Bi=[];
Bigram=[];
%****************************
F1(f1,:) = [Bigram_1];
F2(f1,:) = [Bigram_2];
    %%}
            %--------------------
            f1=f1+1;
          end
      end
    end
    
feature_m1=[F1]; %feature for m1
feature_m2=[F2]; %feature for m2
class=ones(size(feature_m1,1),1);
score1=  indp_LibSVM_fg(feature_m1, class, SVM_model_1_bi_1to2 ); %
MoRFscores_1(kk,1)= max(score1);
score2=  indp_LibSVM_fg(feature_m2, class, SVM_model_2_bi_1to1 ); %
MoRFscores_2(kk,1)= max(score2);
kk=kk+1;
clear sample_d1_seq;
clear feature_seq1;clear feature_seq2;clear feature_seq3;clear feature_seq4;clear F1,clear F2;clear feature_seq5;
end
mat3=[];
MoRFscores=[ MoRFscores_1 MoRFscores_2];
end
%##############################
