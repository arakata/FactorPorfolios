function [output] = mixing(wghat_vec,whhat_vec,C)
	output = wghat_vec + C*whhat_vec;
end