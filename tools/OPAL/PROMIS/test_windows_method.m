function MoRFscores = test_windows_method(profile)
% MoRF prediction algorithm (window method)
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

load SVM_model_1_win_1to1; %load trained models
%--------------------------------------------------------------------------
kk=1;win_flank_siz= 20;  mat4=profile ;T_len=size(mat4,1);
 for e=1:T_len
    mf_1=zeros(((win_flank_siz*2)+1),size(mat4,2)); 
    mf_2=zeros(((win_flank_siz*2)+1),size(mat4,2)); 
    if e<win_flank_siz+1 % seq at start
      if e>1  
      mf_1((win_flank_siz+2)-e:win_flank_siz,:)=  mat4(1:e-1,:); 
      mf_1((win_flank_siz+1):((win_flank_siz*2)+1),:)=  mat4(e:e+win_flank_siz,:);
      end
      if e==1
      mf_1((win_flank_siz+1):((win_flank_siz*2)+1),:)=  mat4(1:e+win_flank_siz,:) ;
      end 
    sample_d= mf_1; 
    elseif e> T_len-win_flank_siz %seq at end
        mf_2(1:(win_flank_siz+1),:)= mat4(e-win_flank_siz:e,:);
        if e~=T_len
        mf_2((win_flank_siz+2):(win_flank_siz+1)+(T_len-e),:)= mat4(e+1:end,:);
        end        
    sample_d= mf_2;
    else % seq in middle
    sample_d= mat4(e-win_flank_siz:e+win_flank_siz,:);
    end
    
 flanksreq=20; ALL_wn41= sample_d(win_flank_siz-flanksreq+1:win_flank_siz+flanksreq+1,:); %win 41
 feature_pp41(kk,:)=  ALL_wn41(:)';
 kk=kk+1;
 end
mat4=[];
%-------------------------------------------------------------------------
% Predict and combine scores 
[predict_label_L, accuracy_L1, dec_values_1] = svmpredict( ones(size(feature_pp41,1),1), feature_pp41, SVM_model_1_win_1to1,['-q -b 1']);
MoRFscores = [ dec_values_1(:,1)] ;
end
%##############################
