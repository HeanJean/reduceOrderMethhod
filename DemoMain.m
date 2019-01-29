format rat;

%% example 1
A = [
    1 0 -1 0 -1
    -1 2 0 -1 -1
    0 -1 3 -1 0
    -1 0 -1 4 -1
    0 -1 -1 0 5];

VectorQ1 = [-1 2 -1 2 1]';
disp(' solutions of example 1: ')
SolutionX = reduceOrderMethod(A, VectorQ1)

%% example 2
VectorQ2 = [-1 1 -1 0 1]';
disp(' solutions of example 2: ')
SolutionX = reduceOrderMethod(A, VectorQ2)



%% example 3
Dim = 100;
ZerosMarix = zeros(Dim, Dim);
M1 = diag(-1 * ones(Dim -1, 1), 1);
M2 = diag(3 * ones(Dim, 1));
Matrix = M1 + M2 + M1';
Q = -1 * ones(Dim, 1);
for ForI = 1 : Dim
    Q(ForI) = Q(ForI) ^(ForI);
end

disp(' solutions of example 3: ')
SolutionX = reduceOrderMethod(Matrix, Q)
