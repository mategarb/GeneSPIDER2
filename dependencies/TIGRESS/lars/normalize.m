function normmat=normalize(inpmat)

n=size(inpmat,1);
inpmat=(eye(n)-1/n*ones(n))*inpmat;

n22 = sqrt(diag(inpmat'*inpmat));
normmat= inpmat./(ones(size(inpmat,1),1)*n22');