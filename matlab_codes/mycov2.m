function [output] = mycov2(x_mat)
[NTIME, NASSET] = size(x_mat);
	
auxacc = ones(NASSET,NASSET)*0;
emu_vec = mean(x_mat);

for i=1:NTIME 
	aux=(x_mat(i,:)'-emu_vec')*((x_mat(i,:)'-emu_vec')');
	auxacc=auxacc+aux;
end

output = auxacc/NTIME; 
output = (NTIME/(NTIME-NASSET-1))*output;
end
