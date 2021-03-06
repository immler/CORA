function res = test_tan(~)
% test_tan - unit_test_function comparison to IntLabV6
%
% Syntax:  
%    res = test_tan(intVal)
%
% Inputs:
%    intVal - interval object
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Daniel Althoff
% Written:      03-November-2015
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
tol = 1e-9;
res = true;

int0 = tan(interval(0,0));
if abs(int0.inf-0)>tol || abs(int0.sup-0)>tol
	res = false;
	disp('test_tan failed');
	return;
end

int0 = tan(interval(-1,0));
if abs(int0.inf+1.55740772465491)>tol || abs(int0.sup-0)>tol
	res = false;
	disp('test_tan failed');
	return;
end

int0 = tan(interval(-2,0));
if ~isinf(int0.inf) || ~isinf(int0.sup)
	res = false;
	disp('test_tan failed');
	return;
end

%% ToDo check this
tol_big=35;
int0 = tan(interval(-pi/2,0));
if abs(int0.inf + 1.0e+016*1.63312393531954) >tol_big || abs(int0.sup-0)>tol
	res = false;
	disp('test_tan failed');
	return;
end
 

int0 = tan(interval(-pi/2-1e-9,0));
if  ~isinf(int0.inf) || ~isinf(int0.sup)
	res = false;
	disp('test_tan failed');
	return;
end

int0 = tan(interval(1.0,1.0));
if  abs(int0.inf-1.55740772465491)>tol || abs(int0.sup-1.55740772465491)>tol
	res = false;
	disp('test_tan failed');
	return;
end

int0 = tan(interval(-pi/2,+pi/2));
if  ~isinf(int0.inf) || ~isinf(int0.sup)
	res = false;
	disp('test_tan failed');
	return;
end


int0 = tan(interval(-pi/2,pi/2-1e-9));
if  abs(int0.inf+1.0e+016*1.63312393531954)>tol_big || abs(int0.sup-1.000000078071901e+09)>tol*1000
	res = false;
	return;
end

disp('test_tan successful');
return;

% tan(infsup(0,0)) = 0;
% tan(infsup(-1,0)) = [  -1.55740772465491,   0.00000000000000] 
% tan(infsup(-2,0)) = +/- inf
% tan(infsup(-pi/2,0)) =  1.0e+016 * [  -1.63312393531954,   0.00000000000000]
% tan(infsup(-pi/2-1e-9,0)) = +/-Inf 
% tan(infsup(1,1)) = 1.55740772465491
% tan(infsup(-pi/2,pi/2)) = +/-Inf
% tan(infsup(-pi/2,pi/2-1e-9)) =  1.0e+016 * [  -1.63312393531954,   0.00000010000001] 







% res = interval();
% 
% res.inf = tan(x.inf);
% res.sup = tan(x.sup);
%------------- END OF CODE --------------

