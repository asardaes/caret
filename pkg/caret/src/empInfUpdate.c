#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>

SEXP empInfUpdate(SEXP matrix, SEXP nrow, SEXP id_rows, SEXP id_columns, SEXP values)
{
	double *p_matrix, *p_values;
	int *p_rows, *p_columns;
	int nrows, i, j, k, I, J;
	
	I = LENGTH(id_rows);
	J = LENGTH(id_columns);
	nrows = asInteger(nrow);
	
	p_matrix = REAL(matrix);
	p_rows = INTEGER(id_rows);
	p_columns = INTEGER(id_columns);
	p_values = REAL(values);
	
	k = 0;
	for(i = 0; i < I; i++)
	{
		for(j = 0; j < J; j++)
		{
			p_matrix[(p_rows[i] - 1) + nrows * (p_columns[j] - 1)] = p_matrix[(p_rows[i] - 1) + nrows * (p_columns[j] - 1)] + p_values[k];
			k++;
		}
	}
	
	return R_NilValue;
}