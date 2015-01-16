ndimgauss.o: ndimgauss_scalar.f GaussNDim.f
	f77 -c -recursive ndimgauss_scalar.f

GaussNDim.f: GaussPointsNDim.m
	math -run "Points = {6, 6, 6, 8, 8, 8, 8}; << GaussPointsNDim.m; Exit[]"

