#include <vector>
#include <span>
#include <iostream>

// 定义一个不规则二维数组的类
class JaggedArray {
public:
    // 添加一行数据
    void addRow(const std::vector<int>& row) {
        startPositions.push_back(data.size());
        data.insert(data.end(), row.begin(), row.end());
        rowSizes.push_back(row.size());
    }

    // 获取指定行的视图
    std::span<int> getRow(size_t rowIndex) {
        if (rowIndex >= startPositions.size()) {
            throw std::out_of_range("Row index out of range");
        }
        return std::span<int>(data.data() + startPositions[rowIndex], rowSizes[rowIndex]);
    }

    // 打印所有行数据
    void print() {
        for (size_t i = 0; i < startPositions.size(); ++i) {
            auto row = getRow(i);
            std::cout << "Row " << i << ": ";
            for (auto& elem : row) {
                std::cout << elem << " ";
            }
            std::cout << std::endl;
        }
    }

private:
    std::vector<int> data;                 // 存储所有数据
    std::vector<size_t> startPositions;    // 存储每行的起始位置
    std::vector<size_t> rowSizes;          // 存储每行的大小
};
