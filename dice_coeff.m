function dice = dice_coeff(A,B)
% function to calculat the Dice coefficient beetween 2 matrices
% output is a single value
    dice = 2*nnz(A & B)/(nnz(A)+nnz(B));
end
