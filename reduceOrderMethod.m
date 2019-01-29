function [SolutionX] = reduceOrderMethod(MatrixA, VectorQ)
% This is the matlab script for "Ximing Fang, Zhijun Qiao, The reduced
% order method for solving the linear complementarity problem with an M-matrix."
% x^T (Ax + q) =  0, x >= 0, Ax + q >=0                  (1)
% input varies:  MatrixA : M-matrix
%                 VectorQ
% output varies: SolutionX : solution of equation (1)


%% Pseudocode 
% Input：MatrixA, VectorQ
% output：SolutionX

% function SolutionX = reduceOrderMethod(MatrixA, VectorQ)
%	if all elements of VectorQ are negative
%		SolutionX <-- inv(MatrixA) * (- VectorQ);
%		return;
%	end if 
%	if all elements of VectorQ are not negative
%		SolutionX <-- 0
%		return;
%	end if  %the upper codes regarding to Lema 2.1
%
%	extract MatrixASnSn, MatrixASpSn, MatrixASpSp and MatrixASnSp from MatrixA according to the locaitons of negative elements in VectorQ
%	extract VectorQ1 and VectorQ2 from VectorQ according to the locations of negative elements in VectorQ;
%		
%	VectorV1 <-- inv(MatrixASnSn) * (-VectorQ1);	
%	VectorV2 <-- VectorQ2 + MatrixASpSn * VectorV1;
%		
%	if all elements of VectorV2 are negative
%		SolutionXPart1 <-- VectorV1;
%		SolutionXPart2 <-- 0;
%		return
%	end if  % those upper codes regarding to Theorem 2.1
%
%	VectorQ <-- VectorV2;
%	MatrixA <-- MatrixASpSp - MatrixASpSn / MatrixASnSn *MatrixASnSp;
%	
%	SolutionXPart2	<-- reduceOrderMethod (MatrixA, VectorQ);
%	VectorB2 <-- -VectorQ1 - MatrixASnSp * SolutionXPart2
%	SolutionXPart1 <-- inv(MatrixASnSn) * VectorB2;
%	return;  % those upper codes regarding to Theorem 2.2


%% Matlab Code
LenP = length(VectorQ);

[FlagAllNegative, FlagAllNotNegative, NegativeLocation, NotNegativeLocation] = isAllNegative(VectorQ);


if FlagAllNegative
    SolutionX = inv(MatrixA) * (- VectorQ);
    return;
end

if FlagAllNotNegative
    SolutionX = zeros(length(VectorQ), 1);
    return;
end  % the upper codes regarding to Lema 2.1



[MatrixASnSn, MatrixASpSn, MatrixASpSp, MatrixASnSp] = extractElements(MatrixA, NegativeLocation, NotNegativeLocation);
VectorQ1 = VectorQ(NegativeLocation);
VectorQ2 = VectorQ(NotNegativeLocation);

VectorV1 = (MatrixASnSn) \ (-VectorQ1);
VectorV2 = VectorQ2 + MatrixASpSn * VectorV1;

[~, FlagAllNotNegative, ~, ~] = isAllNegative(VectorV2);
if FlagAllNotNegative
    SolutionX = assignSolution(VectorV1, NegativeLocation, length(VectorQ));
    return;
end  % those upper codes regarding to Theorem 2.1


VectorQ = VectorV2;
MatrixA = MatrixASpSp - MatrixASpSn / MatrixASnSn * MatrixASnSp;

SolutionX2 = reduceOrderMethod(MatrixA, VectorQ);

VectorB1 = -VectorQ1 - MatrixASnSp * SolutionX2;
SolutionX1 =  MatrixASnSn \ VectorB1 ;

SolutionX = assignSolution2(SolutionX1, SolutionX2, NegativeLocation, LenP); 
return;% those upper codes regarding to Theorem 2.2

 end



function output = assignSolution2(Vector1, Vector2, Location, Len)
output = zeros(Len,1);
for ForI = 1 : length(Location)
    output(Location(ForI)) = Vector1(ForI);
end

Location2 = setdiff(1 : Len, Location);
for ForI = 1 : length(Location2)
    output(Location2(ForI)) = Vector2(ForI);
end


end



function output = assignSolution(Vector, Location, Len)
output = zeros(Len, 1);
for ForI = 1 : length(Location)
    output(Location(ForI)) = Vector(ForI);
end
end


function [FlagAllNegative, FlagAllNotNegative, NegativeLocation, NotNegativeLocation] = isAllNegative(Vector)
FlagAllNegative = 0;
FlagAllNotNegative = 0;

NegativeLocation = find(Vector < 0);

if length(NegativeLocation) == length(Vector)
    FlagAllNegative = 1;
end

if isempty(NegativeLocation)
    FlagAllNotNegative = 1;
end

NotNegativeLocation = setdiff(1 : length(Vector), NegativeLocation);
end

function [MatrixA, MatrixB, MatrixC, MatrixD] =  extractElements(InputMatrix, LocationA, LocationB)

LenA = length(LocationA);
MatrixA = zeros(LenA, LenA);


for ForI = 1 : LenA
    for ForJ = 1 : LenA
        MatrixA(ForI, ForJ) = InputMatrix(LocationA(ForI), LocationA(ForJ));
    end
end

LenB = length(LocationB);
MatrixB = zeros(LenB, LenA);
for ForI = 1 : LenB
    for ForJ = 1 : LenA
        MatrixB(ForI, ForJ) = InputMatrix(LocationB(ForI), LocationA(ForJ));
    end
end

MatrixC = zeros(LenB, LenB);
for ForI = 1 : LenB
    for ForJ = 1 : LenB
        MatrixC(ForI, ForJ) = InputMatrix(LocationB(ForI), LocationB(ForJ));
    end
end


MatrixD = zeros(LenA, LenB);
for ForI = 1 : LenA
    for ForJ = 1 : LenB
        MatrixD(ForI, ForJ) = InputMatrix(LocationA(ForI), LocationB(ForJ));
    end
end

end