// Input/Output functions
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// 16-Jun-2011

#include "resampTest.h"
#include <stdio.h>
#include <string.h>

int vector_filesize(FILE *f)
{
	float n;
	int size = 0;

	while (EOF != fscanf(f, "%g", &n)) size++;

	rewind(f);
        return size;
}

void matrix_filesize(FILE *f, int * row, int * col)
{
	char line[MAX_LINE_LENGTH];

	fgets(line, MAX_LINE_LENGTH-1, f);
	*row = 1;
	strtok(line, " \t");
	*col = 0;

        // Note that here is a potential bug:
	// The line must end with an empty space, otherwise 
	// the last number will not be counted in
	while (NULL != strtok(NULL, " \t")) (*col)++;             

	while (NULL != fgets(line, MAX_LINE_LENGTH-1, f)) (*row)++;
      
//        printf("col=%d, row=%d\n", *(col), *(row));	

	rewind(f);

}


gsl_matrix * load_m(const char * file)
{
	FILE * f = fopen(file, "r");
	int row, col;
	gsl_matrix * out;

	matrix_filesize(f, &row, &col);

//        printf("matrix file load size: %d rows, %d cols\n", row, col);

	out = gsl_matrix_alloc(row, col);

	gsl_matrix_fscanf(f, out);

	fclose(f);

	return out;
}

gsl_vector * load_v(const char * file)
{
	FILE * f = fopen(file, "r");
	int size = vector_filesize(f);
	gsl_vector * out;

//        printf("size=%d.\n", size);
	out = gsl_vector_alloc(size);

	gsl_vector_fscanf(f, out);
	fclose(f);

	return out;
}

void displaymatrix(gsl_matrix * m, const char * name)
{
	size_t i, j;

	printf("%s =\n", name);
	for (i = 0; i < m->size1; i++)
	{
	    for (j = 0; j < m->size2; j++)
		printf("%.2f\t", gsl_matrix_get(m, i, j));
	    printf("\n");
	}

}

void displayvector(gsl_vector * v, const char * name)
{
	printf("%s =\n", name);
	for (size_t i = 0; i < v->size; i++)
	    printf("%.2f ", gsl_vector_get(v, i));
	printf("\n");

}

int getBootID (mv_Method *tm, char *fname, gsl_matrix *bootID)
{
    if (tm->reprand == TRUE ) {
       bootID = load_m(fname);
       if ((tm->resamp!=SCOREZ) & (tm->resamp!=SCOREBOOT))
          gsl_matrix_add_constant(bootID, -1.0); // Matlab id ->C id 
    }
    return 0;
}

