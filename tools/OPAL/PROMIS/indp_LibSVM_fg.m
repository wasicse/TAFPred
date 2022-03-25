
function indp_result =  indp_LibSVM_fg(Test_m,Test_c,model)
% indp_LibSVM_fg  - call LIBSVM to predict
% Created on: 01/04/2017
% Author: Ronesh Sharma
%
[predict_label_L, accuracy_L, dec_values_L] = svmpredict( Test_c,Test_m,model,['-b 1 -q']);
indp_result = [ dec_values_L(:,1)'  ];
end
%##############################