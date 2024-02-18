function res = func(A_Omega, X, Omega)
	res = 0.5*norm(A_Omega - X(Omega))^2;
end