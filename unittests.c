#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tests/utest.h"

#include "linalg.h"
#include "mtxio.h"

double norm_difference(double *x, double *x_sol, int size)
{
    double diff = 0;
    for (int n = 0; n < size; n++)
    {
	diff += pow(x[n] - x_sol[n], 2);
    }

    return sqrt(diff);
}

UTEST(mtxio, read_2d_mtx)
{

    double A_sol[MAX_SIZE][MAX_SIZE] = {
	{40, 35, 47, 36},
	{2, 42, 62, 14},
	{31, -28, 50, 23},
	{29, 11, -3, 25}
    };
    
    char filename[] = "tests/A_io.mtx";
    double A[MAX_SIZE][MAX_SIZE];

    int size_2d = read_2d_mtx(filename, A);

    ASSERT_EQ(size_2d, 4);
    
    double diff = 0;
    for (int i = 0; i < size_2d; i++)
    {
	for (int j = 0; j < size_2d; j++)
	{
	    diff = pow(A[i][j] - A_sol[i][j], 2);
	}
    }
    double norm = sqrt(diff);
    
    ASSERT_LT(norm, 10e-6);    
}

UTEST(mtxio, read_1d_mtx)
{

    double b_sol[MAX_SIZE] = {2, 29, 12, 6};
    
    char filename[] = "tests/b_io.mtx";
    double b[MAX_SIZE];

    int size_1d = read_1d_mtx(filename, b);

    ASSERT_EQ(size_1d, 4);
    
    double diff = 0;
    for (int i = 0; i < size_1d; i++)
    {
	diff += pow(b[i] - b_sol[i], 2);
    }
    double norm = sqrt(diff);

    ASSERT_LT(norm, 10e-6);
    
}

UTEST(mtxio, write_1d_mtx)
{
    char filename_in[] = "tests/b_io.mtx";
    double b[MAX_SIZE];
    int size_1d = read_1d_mtx(filename_in, b);

    char filename_out[] = "tests/b_tmp.mtx";
    write_1d_mtx(filename_out, b, size_1d);

    char command[50];
    sprintf(command, "diff -w %s %s ",filename_in, filename_out);
    int ret = system(command);
    ASSERT_TRUE(ret==0);
    
    sprintf(command, "rm %s",filename_out);
    ret = system(command);
    (void)ret;
}

UTEST(linalg, det)
{
    int ret = system("python tests/python_linalg.py > tmp.txt");
    char filename_in[] = "A_test.mtx";

    
    double A[MAX_SIZE][MAX_SIZE];
    int size_2d = read_2d_mtx(filename_in, A);
    
    double det_A = det(A, size_2d);

    FILE *det_file = fopen("tmp.txt","r");
    double det_A_sol;

    if (fscanf(det_file,"%lf",&det_A_sol) == -1)
    {
	printf("ERROR: failed to read determinant file\n");
	exit(1);
    }

    fclose(det_file);
  
    double diff = fabs((det_A - det_A_sol)/det_A_sol);
    ASSERT_LT(diff, 10e-12);

    ret = system("rm tmp.txt A_test.mtx b_test.mtx x_test.mtx");
    (void)ret;
    
}

UTEST(linalg, solver)
{
    int ret = system("python tests/python_linalg.py > tmp.txt");
    char filename_in_A[] = "A_test.mtx";
    char filename_in_b[] = "b_test.mtx";
    char filename_in_x[] = "x_test.mtx";

    double A[MAX_SIZE][MAX_SIZE];
    read_2d_mtx(filename_in_A, A);

    double b[MAX_SIZE];
    read_1d_mtx(filename_in_b, b);

    double x_sol[MAX_SIZE];
    int size_x_sol = read_1d_mtx(filename_in_x, x_sol);

    double x[MAX_SIZE];
    solve(x, A, b, size_x_sol);

    double diff = 0;
    for (int i = 0; i < size_x_sol; i++)
    {
	diff += pow(x[i] - x_sol[i], 2);
    }
    double norm = sqrt(diff);

    ASSERT_LT(norm, 10e-6);

    ret = system("rm tmp.txt A_test.mtx b_test.mtx x_test.mtx");
    (void)ret;
}

UTEST(main, usage)
{
    int out = system("./solver > tmp_1.txt");

    ASSERT_EQ(out, 0);

    FILE *fp = fopen("tmp_2.txt","w");

    fprintf(fp,"Usage:\n    $ ./solver <A_file> <b_file> "
	    "<x_file>\n");
    fclose(fp);

    out = system("diff -w tmp_1.txt tmp_2.txt");

    ASSERT_EQ(out, 0);

    int ret = system("rm tmp_1.txt tmp_2.txt");
    (void)ret;

}


UTEST(main, solution)
{
    int ret = system("python tests/python_linalg.py > tmp_1.txt");
    int out = system("./solver A_test.mtx b_test.mtx "
		     "x.mtx > tmp_2.txt");

    ASSERT_EQ(out, 0);
    
    char filename_in_x[] = "x.mtx";    
    double x[MAX_SIZE]; 
    int size = read_1d_mtx(filename_in_x, x);

    char filename_in_x_sol[] = "x_test.mtx";    
    double x_sol[MAX_SIZE]; 
    read_1d_mtx(filename_in_x_sol, x_sol);

    double diff = norm_difference(x, x_sol, size);
    
    ASSERT_LT(diff, 10e-12);
    
    ret = system("rm A_test.mtx b_test.mtx x_test.mtx "
		     "tmp_1.txt tmp_2.txt x.mtx");
    (void)ret;
}

UTEST(main, output)
{

    int ret = system("python tests/python_linalg.py > tmp_1.txt");
    int out = system("./solver A_test.mtx b_test.mtx "
		     "x_sol.mtx > tmp_2.txt");

    ASSERT_EQ(out, 0);

    char mat_size_str[10];
    int mat_size;
    
    FILE *fp = fopen("tmp_2.txt","r");
    if (fscanf(fp,"%*s %s %*s",mat_size_str) == -1)
    {
	printf("ERROR:  fscanf() failed to read file.\n");
	exit(1);
    }

    mat_size = atoi(mat_size_str);
    
    double time;
    if (fscanf(fp,"%*s %lf %*s",&time) == -1)
    {
	printf("ERROR:  fscanf() failed to read file.\n");
	exit(1);
	    
    }
    fclose(fp);
    
    ASSERT_LT(time, 3);

    FILE *fout = fopen("tmp_3.txt","w");
    fprintf(fout,"solving %ix%i system\n",
	    mat_size,mat_size);
    fprintf(fout,"done...took %f seconds",time);
    fclose(fout);

    out = system("diff -w tmp_2.txt tmp_3.txt");
    ASSERT_EQ(out, 0);

    ret = system("rm A_test.mtx b_test.mtx x_test.mtx "
		     "tmp_1.txt tmp_2.txt tmp_3.txt x_sol.mtx");
    (void)ret;

}



UTEST_MAIN();
