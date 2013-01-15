% simple script to plot the candidate function for optimization using GAS
LOW = -1;
HI = 2;
% pointer to candidate function
fct=@cand_func;
fplot(fct, [LOW HI]);
