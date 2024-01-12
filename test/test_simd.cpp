#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#include <iostream>
#include <chrono>
#include <cblas.h>
#include <immintrin.h>
#include <algorithm>
#include <omp.h>
#include <Eigen/Dense> // 添加Eigen头文件


struct AE{
  union { __m256 v; struct {double x, y, z, w; }; double d[4];};
};


// 普通矩阵相乘
void matrix_multiply_normal(const double* A, const double* B, double* C, int rows, int cols, int common) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double sum = 0.0f;
            for (int k = 0; k < common; ++k) {
                sum += A[i * common + k] * B[j * cols + k];
            }
            C[i * cols + j] = sum;
        }
    }
}

// SIMD 矩阵相乘
void matrix_multiply_simd(const double* A, const double* B, double* C, int rows, int cols, int common) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            __m256d sum = _mm256_setzero_pd();
            for (int k = 0; k < common; k += 4) {
                __m256d a = _mm256_loadu_pd(&A[i * common + k]);
                __m256d b = _mm256_loadu_pd(&B[k * cols + j]);
                sum = _mm256_add_pd(sum, _mm256_mul_pd(a, b));
            }
            double result[4];
            _mm256_storeu_pd(result, sum);
            C[i * cols + j] = result[0] + result[1] + result[2] + result[3];
        }
    }
}

void matrix_multiply_blas(const double* A, const double* B, double* C, int rows, int cols, int common) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, cols, common, 1.0, A, common, B, cols, 0.0, C, cols);
}

// Eigen 矩阵相乘
void matrix_multiply_eigen(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, Eigen::MatrixXd& C) {
    C = A * B;
}

// 生成随机矩阵
void generate_random_matrix(double* matrix, int rows, int cols) {
    for (int i = 0; i < rows * cols; ++i) {
        matrix[i] = static_cast<double>(rand()) / RAND_MAX;
    }
}

int main() {
    const uint32_t rows = 100;
    const uint32_t cols = 100;
    const uint32_t common = 100;
    const uint32_t cnum = 10000;

    uint32_t nA = rows * common;
    uint32_t nB = cols * common;
    uint32_t nC = cols * rows;

    double* A = new double[cnum*nA];
    double* B = new double[cnum*nB];
    double* C_normal = new double[cnum*rows * cols];
    double* C_simd = new double[cnum*rows * cols];
    double* C_blas = new double[cnum*rows * cols];

    // Eigen 矩阵
    std::vector<Eigen::MatrixXd> eigen_A(cnum, Eigen::MatrixXd(rows, common));
    std::vector<Eigen::MatrixXd> eigen_B(cnum, Eigen::MatrixXd(common, cols));
    std::vector<Eigen::MatrixXd> eigen_C(cnum, Eigen::MatrixXd(rows, cols));

    // 生成随机矩阵
    auto start_gen = std::chrono::high_resolution_clock::now();
    //#pragma omp parallel for
    for(uint32_t i = 0; i < cnum; i++)
    {
      generate_random_matrix(&A[i*nA], rows, common);
      generate_random_matrix(&B[i*nB], common, cols);
      eigen_A[i] = Eigen::Map<Eigen::MatrixXd>(&A[i*nA], rows, common);
      eigen_B[i] = Eigen::Map<Eigen::MatrixXd>(&B[i * nB], common, cols);
    }
    auto end_gen = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_gen = end_gen - start_gen;
    std::cout << "Normal matrix multiplication time: " << elapsed_gen.count() << " seconds\n";

    // 普通矩阵相乘
    auto start_normal = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for(uint32_t i = 0; i < cnum; i++)
      matrix_multiply_normal(&A[i*nA], &B[i*nB], &C_normal[i*nC], rows, cols, common);
    auto end_normal = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_normal = end_normal - start_normal;
    std::cout << "Normal matrix multiplication time: " << elapsed_normal.count() << " seconds\n";

    // SIMD 矩阵相乘
    auto start_simd = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for(uint32_t i = 0; i < cnum; i++)
      matrix_multiply_simd(&A[i*nA], &B[i*nB], &C_simd[i*nC], rows, cols, common);
    auto end_simd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_simd = end_simd - start_simd;
    std::cout << "SIMD matrix multiplication time: " << elapsed_simd.count() << " seconds\n";

    // BLAS 矩阵相乘
    auto start_blas = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for(uint32_t i = 0; i < cnum; i++)
      matrix_multiply_blas(&A[i*nA], &B[i*nB], &C_blas[i*nC], rows, cols, common);
    auto end_blas = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_blas = end_blas - start_blas;
    std::cout << "BLAS matrix multiplication time: " << elapsed_blas.count() << " seconds\n";

    // BLAS 矩阵相乘
    auto start_eigen = std::chrono::high_resolution_clock::now();
    //#pragma omp parallel for
    //for(uint32_t i = 0; i < cnum; i++)
    //{
      //matrix_multiply_eigen(eigen_A[i], eigen_B[i], eigen_C[i]);
      //eigen_A[i] = eigen_A[i].transpose();
    //}
    for(uint32_t i = 0; i < cnum; i++)
    {
      eigen_A[i] = eigen_A[i].transpose();
    }
    auto end_eigen = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_eigen = end_eigen - start_eigen;
    std::cout << "Eigen matrix multiplication time: " << elapsed_eigen.count() << " seconds\n";

    // 验证结果是否相同
    // for (uint32_t i = 0; i < cnum*rows * cols; ++i) {
    //     if (C_normal[i] != eigen_C[i]) {
    //         std::cerr << "Results do not match!\n";
    //         break;
    //     }
    // }
    AE ae;
    std::cout << ae.x << " " << ae.y << std::endl;

    // 释放内存
    delete[] A;
    delete[] B;
    delete[] C_normal;
    delete[] C_simd;

    return 0;
}

