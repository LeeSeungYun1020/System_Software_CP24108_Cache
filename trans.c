/* Lee Seung Yun, 201645825 */
/*
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */
#include <stdio.h>
#include "cachelab.h"

void transpose_submit(int M, int N, int A[N][M], int B[M][N]);
void transpose_N32(int idxRow, int idxCol, int N, int A[N][N], int B[N][N]);
void transpose_N64(int idxRow, int idxCol, int N, int A[N][N], int B[N][N]);
void transpose_NM(int idxRow, int idxCol, int M, int N, int A[N][M], int B[M][N]);
void transpose_additional(int M, int N, int A[N][M], int B[M][N]);
void trans(int M, int N, int A[N][M], int B[M][N]);
void registerFunctions();
int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/*
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded.
 */
char transpose_submit_desc[] = "Transpose submission";

void transpose_submit(int M, int N, int A[N][M], int B[M][N]) {
	int i, j;
	if (N == M) {
		for (i = 0; i < N; i += 8) {
			for (j = 0; j < N; j += 8) {
				if (N == 32) {
					transpose_N32(j, i, N, A, B);
				} else if (N == 64) {
					transpose_N64(j, i, N, A, B);
				}

			}
		}
	} else {
		for (i = 0; i < M; i += 16) {
			for (j = 0; j < N; j += 16) {
				transpose_NM(j, i, M, N, A, B);
			}
		}
	}
}

void transpose_N32(int idxRow, int idxCol, int N, int A[N][N], int B[N][N]) {
	int value, idx;
	int i, j;
	if (N == 32) {
		for (i = idxRow; i < idxRow + 8; i++) {
			for (j = idxCol; j < idxCol + 8; j++) {
				if (i == j) {
					value = A[i][j];
					idx = i;
				} else
					B[j][i] = A[i][j];
			}
			if (idxRow == idxCol) {
				B[idx][idx] = value;
			}
		}
	}
}

void transpose_N64(int idxRow, int idxCol, int N, int A[N][N], int B[N][N]) {
	int num1, num2, num3, num4, num5, num6, num7, num8;
	int i, j, idx;
	for (i = idxRow; i < idxRow + 4; i++) {
		num1 = A[i][idxCol];
		num2 = A[i][idxCol + 1];
		num3 = A[i][idxCol + 2];
		num4 = A[i][idxCol + 3];
		num5 = A[i][idxCol + 4];
		num6 = A[i][idxCol + 5];
		num7 = A[i][idxCol + 6];
		num8 = A[i][idxCol + 7];
		B[idxCol + 0][i] = num1;
		B[idxCol + 0][i + 4] = num5;
		B[idxCol + 1][i] = num2;
		B[idxCol + 1][i + 4] = num6;
		B[idxCol + 2][i] = num3;
		B[idxCol + 2][i + 4] = num7;
		B[idxCol + 3][i] = num4;
		B[idxCol + 3][i + 4] = num8;
	}
	for (j = idxCol; j < idxCol + 4; j++) {
		idx = j + 4;
		num1 = B[j][idxRow + 4];
		num2 = B[j][idxRow + 5];
		num3 = B[j][idxRow + 6];
		num4 = B[j][idxRow + 7];
		B[j][idxRow + 4] = A[idxRow + 4][j];
		B[j][idxRow + 5] = A[idxRow + 5][j];
		B[j][idxRow + 6] = A[idxRow + 6][j];
		B[j][idxRow + 7] = A[idxRow + 7][j];

		B[idx][idxRow] = num1;
		B[idx][idxRow + 1] = num2;
		B[idx][idxRow + 2] = num3;
		B[idx][idxRow + 3] = num4;
		B[idx][idxRow + 4] = A[idxRow + 4][idx];
		B[idx][idxRow + 5] = A[idxRow + 5][idx];
		B[idx][idxRow + 6] = A[idxRow + 6][idx];
		B[idx][idxRow + 7] = A[idxRow + 7][idx];
	}
}

void transpose_NM(int idxRow, int idxCol, int M, int N, int A[N][M], int B[M][N]) {
	int i, j;
	int value = 0, index = 0;
	for (i = idxRow; i < idxRow + 16 && i < N; i++) {
		for (j = idxCol; j < idxCol + 16 && j < M; j++) {
			if (i == j) {
				value = A[i][j];
				index = i;
			} else
				B[j][i] = A[i][j];
		}
		if (idxRow == idxCol) {
			B[index][index] = value;
		}
	}
}

/*
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started.
 */
char transpose_additional_desc[] = "Transpose ADTL";

void transpose_additional(int M, int N, int A[N][M], int B[M][N]) {
	int block = sizeof(int);
	int i, j, b;
	for (i = 0; i < N; i += block) {
		for (j = 0; j < M; ++j) {
			for (b = 0; b < block && i + b < N; ++b) {
				B[j][i + b] = A[i + b][j];
			}
		}
	}
}
/*
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Transpose BASIC";

void trans(int M, int N, int A[N][M], int B[M][N]) {
	int i, j, tmp;

	for (i = 0; i < N; i++) {
		for (j = 0; j < M; j++) {
			tmp = A[i][j];
			B[j][i] = tmp;
		}
	}

}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions() {
	/* Register your solution function */
	registerTransFunction(transpose_submit, transpose_submit_desc);

	/* Register any additional transpose functions */
	registerTransFunction(trans, trans_desc);
	registerTransFunction(transpose_additional, transpose_additional_desc);
}

/*
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N]) {
	int i, j;
	for (i = 0; i < N; i++) {
		for (j = 0; j < M; ++j) {
			if (A[i][j] != B[j][i]) {
				return 0;
			}
		}
	}
	return 1;
}
