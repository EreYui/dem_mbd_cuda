//  matrixtrans.h
//  自定义矩阵和向量类库，用于替代Eigen库的功能

#pragma once
#include <cmath>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <functional>
#include <ostream>

// 定义顺序：先定义 vector2d，（这里插入vector3i的简单声明）再定义 vector3d，再定义 vector4d，再定义 vector3i，最后定义 matrix3d
//---------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------- vector2d ----------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------
//二维双精度浮点向量类
class vector2d {
public:
    double data[2];  ///< 存储向量数据的数组
    
    // 构造函数
    vector2d() = default;
    vector2d(double x, double y) : data{x, y} {}
    
    // 获取指定索引处元素的引用
    double& operator[](int i) { return data[i]; }
    const double& operator[](int i) const { return data[i]; }
    
    vector2d& operator*=(double s) { data[0] *= s; data[1] *= s; return *this; }
    vector2d& operator+=(const vector2d& b) { data[0] += b[0]; data[1] += b[1]; return *this; }
    vector2d& operator-=(const vector2d& b) { data[0] -= b[0]; data[1] -= b[1]; return *this; }
    vector2d operator-() const { return {-data[0], -data[1]}; } // 一元负号运算符（向量取反）

    double norm() const { return std::sqrt(data[0]*data[0] + data[1]*data[1]); }
    
    //计算向量的平方模长
    double squaredNorm() const { return data[0]*data[0] + data[1]*data[1]; }
    
    // 归一化向量（单位向量）
    vector2d normalize() const {
        double len = norm();
        if (len > 1e-10) { vector2d r = *this; r *= (1.0/len); return r; }
        return {0,0};
    }
    
    //返回归一化后的向量副本，原向量不被修改
    vector2d normalized() const {
        double len = norm();
        if (len > 1e-10) { return {data[0] / len, data[1] / len}; }
        return {0,0}; 
    }
    
    //将向量设置为零向量
    void setZero() { data[0]=0; data[1]=0; }
    
    //创建零向量
    static vector2d Zero() { return {0,0}; }
};

inline vector2d operator+(const vector2d& a, const vector2d& b) { return {a[0]+b[0], a[1]+b[1]}; }
inline vector2d operator-(const vector2d& a, const vector2d& b) { return {a[0]-b[0], a[1]-b[1]}; }
inline vector2d operator*(const vector2d& a, double b) { return {a[0]*b, a[1]*b}; }
inline vector2d operator*(double b, const vector2d& a) { return {a[0]*b, a[1]*b}; }
inline vector2d operator/(const vector2d& a, double b) { return {a[0]/b, a[1]/b}; }

// Forward declaration for vector3i since it's defined after vector3d
class vector3i;
//---------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------- vector3d -----------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------
// 三维双精度浮点向量类
class vector3d {
public:
    double data[3];  // 存储向量数据的数组
    
    // 构造函数
    vector3d() = default;
    vector3d(double x, double y, double z) : data{x, y, z} {}
    
    // 获取指定索引处元素的引用
    double& operator[](int i) { return data[i]; }
    const double& operator[](int i) const { return data[i]; }
    
    vector3d& operator*=(double s) {data[0] *= s; data[1] *= s; data[2] *= s;return *this;}
    vector3d& operator+=(const vector3d& b) {data[0] += b[0]; data[1] += b[1]; data[2] += b[2];return *this;}
    vector3d& operator-=(const vector3d& b) {data[0] -= b[0]; data[1] -= b[1]; data[2] -= b[2];return *this;}
    vector3d operator-() const { return {-data[0], -data[1], -data[2]}; } // 一元负号运算符（向量取反）

    // 计算两个向量的点积
    double dot(const vector3d& b) const {return data[0] * b[0] + data[1] * b[1] + data[2] * b[2];}
    
    // 计算两个向量的叉积
    vector3d cross(const vector3d& b) const {return {data[1] * b[2] - data[2] * b[1],data[2] * b[0] - data[0] * b[2],data[0] * b[1] - data[1] * b[0]};}
    
    // 计算向量的模长（长度）
    double norm() const {return std::sqrt(this->dot(*this));}

    // 求值运算符
    const vector3d& eval() const { return *this; }
    
    // 计算向量的平方模长
    double squaredNorm() const { return this->dot(*this); }

    // 归一化向量（单位向量）
    vector3d normalize() const {
        double len = norm();
        if (len > 1e-14) {
            vector3d result = *this;
            result *= (1.0 / len);
            return result;
        }
        return {0, 0, 0};
    }
    
    // 返回归一化后的向量副本
    vector3d normalized() const { return normalize(); }

    // 计算与另一个向量的逐元素最小值
    vector3d min(const vector3d& b) const {return {std::fmin(data[0], b[0]), std::fmin(data[1], b[1]), std::fmin(data[2], b[2])};}
    
    // 计算与另一个向量的逐元素最大值
    vector3d max(const vector3d& b) const {return {std::fmax(data[0], b[0]), std::fmax(data[1], b[1]), std::fmax(data[2], b[2])};}

    // 将向量设置为零向量
    void setZero() { data[0] = 0; data[1] = 0; data[2] = 0; }
    
    // 创建零向量
    static vector3d Zero() { return {0, 0, 0}; }
    
    // 创建所有元素都为相同值的向量
    static vector3d Constant(double value) { return {value, value, value}; }

    // 转置操作（返回副本）
    vector3d transpose() const {return *this;}
    
    // 逐元素大于等于比较运算符
    bool operator>=(const vector3d& other) const {return data[0] >= other[0] && data[1] >= other[1] && data[2] >= other[2];}

    // 逐元素大于比较运算符
    bool operator>(const vector3d& other) const {return data[0] > other[0] && data[1] > other[1] && data[2] > other[2];}
    
    // 逐元素小于等于比较运算符
    bool operator<=(const vector3d& other) const {return data[0] <= other[0] && data[1] <= other[1] && data[2] <= other[2];}
    
    // 逐元素小于比较运算符
    bool operator<(const vector3d& other) const {return data[0] < other[0] && data[1] < other[1] && data[2] < other[2];}    

    // 检查向量中是否存在非零分量
    bool any() const {return data[0] != 0 || data[1] != 0 || data[2] != 0;}
    
    // 检查是否存在分量大于等于另一个向量的对应分量
    bool any_ge(const vector3d& other) const {return data[0] >= other[0] || data[1] >= other[1] || data[2] >= other[2];}
    
    // 检查是否存在分量小于等于另一个向量的对应分量
    bool any_le(const vector3d& other) const {return data[0] <= other[0] || data[1] <= other[1] || data[2] <= other[2];}
    
    // 检查是否存在分量等于另一个向量的对应分量
    bool any_eq(const vector3d& other) const {return data[0] == other[0] || data[1] == other[1] || data[2] == other[2];}
    
    // 检查是否存在分量不等于另一个向量的对应分量
    bool any_ne(const vector3d& other) const {return data[0] != other[0] || data[1] != other[1] || data[2] != other[2];}
    
    // 转换为整数数组, 向上取整
    vector3i ceiled() const;

    // 转换为整数数组, 向下取整
    vector3i floored() const;
};
// Define necessary operators before rotateVector function
inline vector3d operator+(const vector3d& a, const vector3d& b) {return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};}
inline vector3d operator-(const vector3d& a, const vector3d& b) {return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};}
inline vector3d operator*(const vector3d& a, double b) {return {a[0] * b, a[1] * b, a[2] * b};}
inline vector3d operator*(double b, const vector3d& a) {return {a[0] * b, a[1] * b, a[2] * b};}
inline vector3d operator/(const vector3d& a, double b) {return {a[0] / b, a[1] / b, a[2] / b};}
inline vector3d operator+(const vector3d& a, double b) {return {a[0] + b, a[1] + b, a[2] + b};}
inline vector3d operator-(const vector3d& a, double b) {return {a[0] - b, a[1] - b, a[2] - b};}
inline vector3d operator+(double b, const vector3d& a) {return {a[0] + b, a[1] + b, a[2] + b};}
inline vector3d operator-(double b, const vector3d& a) {return {b - a[0], b - a[1], b - a[2]};}


//---------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------- vector4d -----------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------
// 4D double precision floating point vector class
class vector4d {
public:
    double data[4];  ///< 存储向量数据的数组
    
    // 构造函数
    vector4d() = default;
    vector4d(double x, double y, double z, double w) : data{x, y, z, w} {}
    
    // 获取指定索引处元素的引用
    double& operator[](int i) { return data[i]; }
    const double& operator[](int i) const { return data[i]; }
    
    vector4d& operator*=(double s) { data[0] *= s; data[1] *= s; data[2] *= s; data[3] *= s; return *this; }
    vector4d& operator+=(const vector4d& b) {data[0] += b[0]; data[1] += b[1]; data[2] += b[2]; data[3] += b[3];return *this;}
    vector4d& operator-=(const vector4d& b) {data[0] -= b[0]; data[1] -= b[1]; data[2] -= b[2]; data[3] -= b[3];return *this;}

    // 计算向量的模长（长度）
    double norm() const { return std::sqrt(data[0] * data[0] + data[1] * data[1] + data[2] * data[2] + data[3] * data[3]); }
    
    // 计算向量的平方模长
    double squaredNorm() const { return data[0] * data[0] + data[1] * data[1] + data[2] * data[2] + data[3] * data[3]; }
    
    // 归一化向量（单位向量）
    vector4d normalize() const {
        double len = norm();
        if (len > 1e-10) {
            vector4d result = *this;
            result *= (1.0 / len);
            return result;
        }
        return {0, 0, 0, 0};
    }
    
    // 返回归一化后的向量副本
    vector4d normalized() const { return normalize(); }
    
    // 将向量设置为零向量
    void setZero() { data[0] = 0; data[1] = 0; data[2] = 0; data[3] = 0; }
    
    // 创建零向量
    static vector4d Zero() { return {0, 0, 0, 0}; }
};

inline vector4d operator+(const vector4d& a, const vector4d& b) {return {a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]};}
inline vector4d operator-(const vector4d& a, const vector4d& b) {return {a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]};}
inline vector4d operator*(const vector4d& a, double b) {return {a[0] * b, a[1] * b, a[2] * b, a[3] * b};}
inline vector4d operator*(double b, const vector4d& a) {return {a[0] * b, a[1] * b, a[2] * b, a[3] * b};}
inline vector4d operator/(const vector4d& a, double b) {return {a[0] / b, a[1] / b, a[2] / b, a[3] / b};}
//---------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------- vector3i -----------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------
// 三维整数向量类
class vector3i {
public:
    int data[3];  ///< 存储向量数据的数组
    
    // 构造函数
    vector3i() = default;
    vector3i(int x, int y, int z) : data{x, y, z} {}
    
    // 获取指定索引处元素的引用
    int& operator[](int i) { return data[i]; }
    const int& operator[](int i) const { return data[i]; }

    // 计算两个向量的点积
    int dot(const vector3i& b) const { return data[0] * b[0] + data[1] * b[1] + data[2] * b[2]; }
    
    // 计算两个向量的叉积
    vector3i cross(const vector3i& b) const {return {data[1] * b[2] - data[2] * b[1], data[2] * b[0] - data[0] * b[2], data[0] * b[1] - data[1] * b[0]};}
    
    // 计算与另一个向量的逐元素最小值
    vector3i fmin(const vector3i& b) const { return {std::min(data[0], b[0]), std::min(data[1], b[1]), std::min(data[2], b[2])}; }
    
    // 计算与另一个向量的逐元素最大值
    vector3i fmax(const vector3i& b) const { return {std::max(data[0], b[0]), std::max(data[1], b[1]), std::max(data[2], b[2])}; }
    
    // 将向量设置为零向量
    void setZero() { data[0] = 0; data[1] = 0; data[2] = 0; }
    
    // 创建零向量
    static vector3i Zero() { return {0, 0, 0}; }

    // 相等比较运算符
    bool operator==(const vector3i& other) const {return data[0] == other.data[0] && data[1] == other.data[1] && data[2] == other.data[2];}
    
    // 不等比较运算符
    bool operator!=(const vector3i& other) const {return !(*this == other);}
    
    // Transpose operation (returns copy)
    vector3i transpose() const { return *this; }
};

inline vector3i operator+(const vector3i& a, const vector3i& b) {return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};}
inline vector3i operator-(const vector3i& a, const vector3i& b) {return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};}
inline vector3i operator*(const vector3i& a, int b) {return {a[0] * b, a[1] * b, a[2] * b};}
inline vector3i operator*(int b, const vector3i& a) {return {a[0] * b, a[1] * b, a[2] * b};}
inline vector3i operator+(const vector3i& a, int b) {return {a[0] + b, a[1] + b, a[2] + b};}
inline vector3i operator-(const vector3i& a, int b) {return {a[0] - b, a[1] - b, a[2] - b};}
inline vector3i operator+(int b, const vector3i& a) {return {a[0] + b, a[1] + b, a[2] + b};}
inline vector3i operator-(int b, const vector3i& a) {return {b - a[0], b - a[1], b - a[2]};}



//---------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------- matrix3d ----------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------

// 3x3双精度浮点矩阵类
class matrix3d {
public:
    double data[3][3];  // 存储矩阵数据的二维数组
    
    /** 默认构造函数 */
    matrix3d() = default;
    
    /**
     * 带参数的构造函数
     * @param m00 第一行第一列元素
     * @param m01 第一行第二列元素
     * @param m02 第一行第三列元素
     * @param m10 第二行第一列元素
     * @param m11 第二行第二列元素
     * @param m12 第二行第三列元素
     * @param m20 第三行第一列元素
     * @param m21 第三行第二列元素
     * @param m22 第三行第三列元素
     */
    matrix3d(double m00, double m01, double m02,
             double m10, double m11, double m12,
             double m20, double m21, double m22) {
        data[0][0] = m00; data[0][1] = m01; data[0][2] = m02;
        data[1][0] = m10; data[1][1] = m11; data[1][2] = m12;
        data[2][0] = m20; data[2][1] = m21; data[2][2] = m22;
    }

    // 获取指定行的指针
    double* operator[](int i) { return data[i]; }
    const double* operator[](int i) const { return data[i]; }
    
    // 矩阵标量乘法赋值运算符
    matrix3d& operator*=(double s) {
        for(int i=0; i<3; ++i)
            for(int j=0; j<3; ++j)
                data[i][j] *= s;
        return *this;
    }

    // 矩阵加法赋值运算符
    matrix3d& operator+=(const matrix3d& b) {
        for(int i=0; i<3; ++i)
            for(int j=0; j<3; ++j)
                data[i][j] += b[i][j];
        return *this;
    }
    
    // 矩阵减法赋值运算符
    matrix3d& operator-=(const matrix3d& b) {
        for(int i=0; i<3; ++i)
            for(int j=0; j<3; ++j)
                data[i][j] -= b[i][j];
        return *this;
    }

    // 矩阵转置运算符
    matrix3d transpose() const {
        return matrix3d(
            data[0][0], data[1][0], data[2][0],
            data[0][1], data[1][1], data[2][1],
            data[0][2], data[1][2], data[2][2]
        );
    }
    
    // 矩阵置零运算符
    void setZero() {
        for(int i=0; i<3; ++i)
            for(int j=0; j<3; ++j)
                data[i][j] = 0;
    }
    
    // 创建零矩阵
    static matrix3d Zero() {return matrix3d(0,0,0, 0,0,0, 0,0,0);}
    
    // 创建单位矩阵
    static matrix3d identity() {return matrix3d(1,0,0, 0,1,0, 0,0,1);}
    
    // 创建所有元素都为相同值的矩阵
    static matrix3d Constant(double value) {return matrix3d(value, value, value, value, value, value, value, value, value);}
    
    // 获取指定行作为vector3d
    vector3d row(int i) const {return vector3d(data[i][0], data[i][1], data[i][2]);}
    
    // 获取指定列作为vector3d
    vector3d col(int j) const {return vector3d(data[0][j], data[1][j], data[2][j]);}

    // 获取对角线元素作为vector3d
    vector3d diagonal() const {return vector3d(data[0][0], data[1][1], data[2][2]);}
    
    // 设置对角线元素
    void setDiagonal(const vector3d& diag) {data[0][0] = diag[0]; data[1][1] = diag[1]; data[2][2] = diag[2];}
};

// 矩阵加法运算符
inline matrix3d operator+(const matrix3d& a, const matrix3d& b) {
    matrix3d res;
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            res[i][j] = a[i][j] + b[i][j];
    return res;
}

// 矩阵减法运算符
inline matrix3d operator-(const matrix3d& a, const matrix3d& b) {
    matrix3d res;
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            res[i][j] = a[i][j] - b[i][j];
    return res;
}

// 矩阵标量乘法运算符
inline matrix3d operator*(const matrix3d& a, double s) {matrix3d res = a;res *= s;return res;}
inline matrix3d operator*(double s, const matrix3d& a) {return a * s;}

// 矩阵向量乘法运算符
inline vector3d operator*(const matrix3d& m, const vector3d& v) {
    return {
        m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2],
        m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2],
        m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2]
    };
}

// 矩阵矩阵乘法运算符
inline matrix3d operator*(const matrix3d& a, const matrix3d& b) {
    matrix3d res;
    for(int i=0; i<3; ++i) {
        for(int j=0; j<3; ++j) {
            res[i][j] = a[i][0]*b[0][j] + a[i][1]*b[1][j] + a[i][2]*b[2][j];
        }
    }
    return res;
}



//------------------------------------------------------------Output stream operators------------------------------------------------------------
// Output stream operators
inline std::ostream& operator<<(std::ostream& os, const vector2d& v) {
    os << v[0] << " " << v[1];
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const vector3d& v) {
    os << v[0] << " " << v[1] << " " << v[2];
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const vector4d& v) {
    os << v[0] << " " << v[1] << " " << v[2] << " " << v[3];
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const vector3i& v) {
    os << v[0] << " " << v[1] << " " << v[2];
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const matrix3d& m) {
    for(int i=0; i<3; ++i) {
        for(int j=0; j<3; ++j) {
            os << m[i][j];
            if(j < 2) os << " ";
        }
        if(i < 2) os << "\n";
    }
    return os;
}

inline vector3i vector3d::ceiled() const {return {(int)std::ceil(data[0]), (int)std::ceil(data[1]), (int)std::ceil(data[2])};}
inline vector3i vector3d::floored() const { return {(int)std::floor(data[0]), (int)std::floor(data[1]), (int)std::floor(data[2])}; }




/*
 Quat2DCM: transfer quaternion to direction cosine matrix
 Quat2IvDCM: transfer quaternion to the inverse of direction cosine matrix
           DCM - direction cosine matrix
           the input q must be unit quaternion
*/

/**
 * 将四元数转换为方向余弦矩阵的逆矩阵
 * @param q 输入的单位四元数
 * @return 方向余弦矩阵的逆矩阵
 */
matrix3d Quat2IvDCM(vector4d q);

/**
 * 将四元数转换为方向余弦矩阵
 * @param q 输入的单位四元数
 * @return 方向余弦矩阵
 */
matrix3d Quat2DCM(vector4d q);

/**
 * 对四元数进行归一化
 * @param q 要归一化的四元数（就地修改）
 */
void QuatNorm(vector4d& q);

/**
 * 将三维索引(I,J,K)转换为一维索引
 * @param I 第一个维度的索引
 * @param J 第二个维度的索引
 * @param K 第三个维度的索引
 * @param D0Sz 第一个维度的大小
 * @param D1Sz 第二个维度的大小
 * @param D2Sz 第三个维度的大小
 * @return 一维索引
 */
long long IJK2N(int I, int J, int K, int D0Sz, int D1Sz, int D2Sz);

/**
 * 将一维索引转换为三维索引(I,J,K)
 * @param N 一维索引
 * @param IJK 用于存储结果的三维索引数组
 * @param D0Sz 第一个维度的大小
 * @param D1Sz 第二个维度的大小
 * @param D2Sz 第三个维度的大小
 */
void N2IJK(int N, int IJK[3], int D0Sz, int D1Sz, int D2Sz);

/**
 * 计算伽马函数
 * @param orien 方向向量
 * @param omg 四元数向量
 * @return 伽马函数结果
 */
vector4d gama(vector4d orien, vector3d omg);



