function [dsquaref_mat] = dsquaref(dsquaref_reg1_mat, dsquaref_reg2_mat, dsquaref_reg3_val)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

dsquaref_temp1 = [dsquaref_reg1_mat dsquaref_reg2_mat];
dsquaref_temp2 = [dsquaref_reg2_mat dsquaref_reg3_val];
dsquaref_mat = [dsquaref_temp1; dsquaref_temp2];

end

