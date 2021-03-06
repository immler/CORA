function des_mat= full_fact_mod(levels)
% Full Fact: gives full factorial design matrix for any levels  
% more than 2 of any number of variables (minimum 2)
% usage: des_mat= full_fact(x1,x2,x3);
% usage: des_mat = full_fact([-1 1],[100 200 300],[1:4]);
% 
%  arguments: (input)
%  x1,x2, x3-  variable levels either in row or column vector. These are
%  not number of levels but the levels itself
%  
% arguments: (output)
%  des_mat - mxn m= total number of designs (product of all levels) and  
%           n is number of variables  
%            The first column shows all the first variable levels, second 
%            column shows second variable levels and so on. 
%  
% Example usage:
%  x1=[-1 1];x2=[100:100:300];
%  des_mat = full_fact(x1,x2) % OR
%  des_mat = full_fact([-1 1],[100:100:300])
%  des_mat = 
%     -1   100
%     -1   200
%     -1   300
%      1   100
%      1   200
%      1   300
% 
%  Bhaskar Dongare
%  bhaskar_dongare@yahoo.com
%  Nov-17-2008
%  modified: 08-April-2016, Matthias Althoff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~all(levels >1)
    error 'Each variables should have minimum 2 levels'
end
% Total number of design points  
total=prod(levels);
%Initilization of output matrix
des_mat=zeros(0,1);
% Loop for full factorial points
for i=1:length(levels)
    if i~=1 && i~=length(levels)
        temp=zeros(0,1);
        for j=1:prod(levels(1:i-1))
            temp1=repmat(1:levels(i),prod(levels(i+1:end)),1);
            temp1=sortrows(temp1);
            temp=[temp; temp1];
        end
    elseif i==length(levels)
        temp=repmat((1:levels(i))',total/levels(i),1);
    else
       temp=repmat((1:levels(i))',total/levels(i),1);
       temp=sortrows(temp);
    end
    des_mat=[des_mat temp];
end;
