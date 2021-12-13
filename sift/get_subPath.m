function [subPath,allTrainSamples] = get_subPath(outexProb)
% subpath information
if strcmp(outexProb,'TC10_Train')
    subPath = '000/train.txt';
    allTrainSamples = 480;
elseif strcmp(outexProb,'TC10_Test')
    subPath = '000/test.txt';
    allTrainSamples = 3840;
elseif strcmp(outexProb,'TC22_Train')
    subPath = '000/train.txt';
    allTrainSamples = 1360;
elseif strcmp(outexProb,'TC22_Test')
    subPath = '000/test.txt';
    allTrainSamples = 1360 ;
elseif strcmp(outexProb,'TC23_Train') || strcmp(outexProb,'TC23a_Train') ...
        || strcmp(outexProb,'TC23b_Train') || strcmp(outexProb,'TC23c_Train') ...
        || strcmp(outexProb,'TC23d_Train') || strcmp(outexProb,'TC23f_Train') ...
        || strcmp(outexProb,'TC23g_Train') || strcmp(outexProb,'TC23h_Train') ...
        || strcmp(outexProb,'TC23i_Train') || strcmp(outexProb,'TC23j_Train')...
        || strcmp(outexProb,'TC23o_Train') || strcmp(outexProb,'TC23r_Train') ...
        || strcmp(outexProb,'TC23p_Train') || strcmp(outexProb,'TC23s_Train') ...
        || strcmp(outexProb,'TC23q_Train') || strcmp(outexProb,'TC23t_Train') ...
        || strcmp(outexProb,'TC23u_Train') || strcmp(outexProb,'TC23v_Train') ...
        || strcmp(outexProb,'TC23w_Train') || strcmp(outexProb,'TC23z_Train') ...
        || strcmp(outexProb,'TC23x_Train') || strcmp(outexProb,'TC23k_Train') ...
        || strcmp(outexProb,'TC23y_Train') || strcmp(outexProb,'TC23l_Train')
    subPath = '000/train.txt';
    allTrainSamples = 1360;
elseif strcmp(outexProb,'TC24_Train')
    subPath = '002/train.txt';
    allTrainSamples = 1360;
elseif strcmp(outexProb,'TC24_Test')
    subPath = '002/test.txt';
    allTrainSamples = 2720;
elseif strcmp(outexProb,'TC23_Test') || strcmp(outexProb,'TC23a_Test') ...
        || strcmp(outexProb,'TC23b_Test') || strcmp(outexProb,'TC23c_Test') ...
        || strcmp(outexProb,'TC23d_Test') || strcmp(outexProb,'TC23f_Test') ...
        || strcmp(outexProb,'TC23g_Test') || strcmp(outexProb,'TC23h_Test') ...
        || strcmp(outexProb,'TC23i_Test') || strcmp(outexProb,'TC23j_Test') ...
        || strcmp(outexProb,'TC23o_Test') || strcmp(outexProb,'TC23r_Test') ...
        || strcmp(outexProb,'TC23p_Test') || strcmp(outexProb,'TC23s_Test') ...
        || strcmp(outexProb,'TC23q_Test') || strcmp(outexProb,'TC23t_Test')...
        || strcmp(outexProb,'TC23u_Test') || strcmp(outexProb,'TC23v_Test') ...
        || strcmp(outexProb,'TC23w_Test') || strcmp(outexProb,'TC23z_Test') ...
        || strcmp(outexProb,'TC23x_Test') || strcmp(outexProb,'TC23k_Test') ...
        || strcmp(outexProb,'TC23y_Test') || strcmp(outexProb,'TC23l_Test')
    subPath = '000/test.txt';
    allTrainSamples = 1360 ;
elseif strcmp(outexProb,'TC11n_Train') || strcmp(outexProb,'TC11a_Train') ...
        || strcmp(outexProb,'TC11b_Train') || strcmp(outexProb,'TC11c_Train') ...
        || strcmp(outexProb,'TC11d_Train') || strcmp(outexProb,'TC11e_Train')...
        || strcmp(outexProb,'TC11f_Train') || strcmp(outexProb,'TC11g_Train') ...
        || strcmp(outexProb,'TC11h_Train') || strcmp(outexProb,'TC11i_Train') ...
        || strcmp(outexProb,'TC11j_Train') || strcmp(outexProb,'TC11o_Train') ...
        || strcmp(outexProb,'TC11o_Train') || strcmp(outexProb,'TC11r_Train') ...
        || strcmp(outexProb,'TC11p_Train') || strcmp(outexProb,'TC11s_Train') ...
        || strcmp(outexProb,'TC11q_Train') || strcmp(outexProb,'TC11t_Train') ...
        || strcmp(outexProb,'TC11u_Train') || strcmp(outexProb,'TC11v_Train') ...
        || strcmp(outexProb,'TC11x_Train') || strcmp(outexProb,'TC11w_Train') ...
        || strcmp(outexProb,'TC11y_Train') || strcmp(outexProb,'TC11k_Train') ...
        || strcmp(outexProb,'TC11z_Train') || strcmp(outexProb,'TC11l_Train')
    subPath = '000/train.txt';
    allTrainSamples = 480;
elseif strcmp(outexProb,'TC11n_Test') || strcmp(outexProb,'TC11a_Test') ...
        || strcmp(outexProb,'TC11b_Test') || strcmp(outexProb,'TC11c_Test') ...
        || strcmp(outexProb,'TC11d_Test') || strcmp(outexProb,'TC11e_Test')...
        || strcmp(outexProb,'TC11f_Test') || strcmp(outexProb,'TC11g_Test') ...
        || strcmp(outexProb,'TC11h_Test') || strcmp(outexProb,'TC11i_Test') ...
        || strcmp(outexProb,'TC11j_Test') || strcmp(outexProb,'TC11o_Test')...
        || strcmp(outexProb,'TC11o_Test') || strcmp(outexProb,'TC11r_Test') ...
        || strcmp(outexProb,'TC11p_Test') || strcmp(outexProb,'TC11s_Test') ...
        || strcmp(outexProb,'TC11q_Test') || strcmp(outexProb,'TC11t_Test') ...
        || strcmp(outexProb,'TC11u_Test') || strcmp(outexProb,'TC11v_Test') ...
        || strcmp(outexProb,'TC11x_Test') || strcmp(outexProb,'TC11w_Test') ...
        || strcmp(outexProb,'TC11y_Test') || strcmp(outexProb,'TC11k_Test') ...
        || strcmp(outexProb,'TC11z_Test') || strcmp(outexProb,'TC11l_Test')
    subPath = '000/test.txt';
    allTrainSamples = 480 ;
elseif strcmp(outexProb,'TC12_Prob000_Train')
    subPath = '000/train.txt';
    allTrainSamples = 480;
elseif strcmp(outexProb,'TC12_Prob000_Test')
    subPath = '000/test.txt';
    allTrainSamples = 4320;
elseif strcmp(outexProb,'TC12_Prob001_Train')
    subPath = '001/train.txt';
    allTrainSamples = 480;
elseif strcmp(outexProb,'TC12_Prob001_Test')
    subPath = '001/test.txt';
    allTrainSamples = 4320;
else
    error('No such problem!');
end
% ==================================================
end