function [rodriguesRotMatrix] = rotrodrigues(angle,vector)
%% Return a rotation matrix for rotation around an arbitrary vector.
%
% Input: angle in radians.
%
% References:
% https://math.stackexchange.com/questions/142821
% https://mathworld.wolfram.com/RodriguesRotationFormula.html

% Guarantee normalization of vector into an unit vector.
vector = vector / norm(vector);

I = eye(3); % [3x3] identity matrix.

W = [       0  -vector(3)  vector(2);
     vector(3)         0  -vector(1);
    -vector(2)  vector(1)         0];

rodriguesRotMatrix = I + sin(angle)*W + 2*(sin(angle/2)^2)*(W^2);

end % Function rotrodrigues()
