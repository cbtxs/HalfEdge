#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

using CSRMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;

int main(int, char * argv[]) 
{
    // 读取COO矩阵文件
    std::ifstream cooFile(argv[1]);
    if (!cooFile.is_open()) 
    {
        std::cerr << "Error opening COO matrix file." << std::endl;
        return 1;
    }

    // 读取COO矩阵数据
    int rows = 170187, cols = 170187, nnz = 6103179;
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nnz);

    std::string line;
    while (std::getline(cooFile, line)) 
    {
      std::istringstream lineStream(line);
      double value;
      int row, col;
      lineStream >> value >> row >> col;
      value = std::stod(lineStream.str().substr(0, lineStream.str().find(" ")));
      col = std::stoi(lineStream.str().substr(lineStream.str().rfind(" ") + 1));

      tripletList.push_back(Eigen::Triplet<double>(row, col, value));
    }
    cooFile.close();

    // 使用Triplet表示COO矩阵
    CSRMatrix cooMatrix(rows, cols);
    cooMatrix.setFromTriplets(tripletList.begin(), tripletList.end());

    // 读取一维数组文件
    std::ifstream arrayFile(argv[2]);
    if (!arrayFile.is_open()) {
        std::cerr << "Error opening one-dimensional array file." << std::endl;
        return 1;
    }

    // 读取一维数组数据
    Eigen::VectorXd b(rows);
    for (int i = 0; i < rows; ++i) {
        arrayFile >> b[i];
    }
    arrayFile.close();

    // 使用BiCGStab求解方程组
    clock_t s = clock();
    Eigen::ConjugateGradient<CSRMatrix, Eigen::Lower, Eigen::IncompleteCholesky<double>> solver;
    solver.compute(cooMatrix);
    //solver.setTolerance(1e-18); 

    if (solver.info() != Eigen::Success) {
        std::cerr << "Eigen BiCGStab decomposition failed." << std::endl;
        return 1;
    }

    Eigen::VectorXd x = solver.solve(b);
    clock_t t = clock();
    std::cout << " time : " << (double)(t-s)/CLOCKS_PER_SEC << std::endl;

    if (solver.info() != Eigen::Success) {
        std::cerr << "Eigen BiCGStab solve failed." << std::endl;
        return 1;
    }

    // 打印解向量
    //std::cout << "Solution vector (x):" << std::endl << x << std::endl;

    return 0;
}

