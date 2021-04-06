
FiniteMatrix finiteObj;

//The link coefficients
FiniteMatrix::finiteMat AE(NI, vector<vector<FiniteMatrix> >(NJ, vector<FiniteMatrix> (NK)));
FiniteMatrix::finiteMat AW(NI, vector<vector<FiniteMatrix> >(NJ, vector<FiniteMatrix> (NK)));
FiniteMatrix::finiteMat AN(NI, vector<vector<FiniteMatrix> >(NJ, vector<FiniteMatrix> (NK)));
FiniteMatrix::finiteMat AS(NI, vector<vector<FiniteMatrix> >(NJ, vector<FiniteMatrix> (NK)));
FiniteMatrix::finiteMat AT(NI, vector<vector<FiniteMatrix> >(NJ, vector<FiniteMatrix> (NK)));
FiniteMatrix::finiteMat AB(NI, vector<vector<FiniteMatrix> >(NJ, vector<FiniteMatrix> (NK)));
FiniteMatrix::finiteMat AP(NI, vector<vector<FiniteMatrix> >(NJ, vector<FiniteMatrix> (NK)));

// The source terms
FiniteMatrix::finiteMat SU(NI, vector<vector<FiniteMatrix> >(NJ, vector<FiniteMatrix> (NK)));
FiniteMatrix::finiteMat SV(NI, vector<vector<FiniteMatrix> >(NJ, vector<FiniteMatrix> (NK)));
FiniteMatrix::finiteMat SW(NI, vector<vector<FiniteMatrix> >(NJ, vector<FiniteMatrix> (NK)));

// The link coefficients for boundaries
FiniteMatrix::finiteMat APU(NI, vector<vector<FiniteMatrix> >(NJ, vector<FiniteMatrix> (NK)));
FiniteMatrix::finiteMat APV(NI, vector<vector<FiniteMatrix> >(NJ, vector<FiniteMatrix> (NK)));
FiniteMatrix::finiteMat APW(NI, vector<vector<FiniteMatrix> >(NJ, vector<FiniteMatrix> (NK)));


 

