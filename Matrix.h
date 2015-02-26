#ifndef MATRIX_H
#define MATRIX_H

//#include <algorithm>
#include <array>
#include <iomanip>
#include <sstream>

namespace Linear {

template <typename T, size_t N, bool C>
class Vec;

template <typename T, size_t M, size_t N = M>
class Matrix {
private:
    std::array<std::array<T, N>, M> v;

public:
    Matrix() {
        static_assert(M > 0 && N > 0, "Matrix cannot have a zero dimension");
    }

    template <typename U>
    Matrix(const Matrix<U, M, N>& other) {
        static_assert(M > 0 && N > 0, "Matrix cannot have a zero dimension");
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                v[i][j] = other[i][j];
            }
        }
    }

    template <typename U>
    Matrix(std::initializer_list<std::initializer_list<U>> list) {
        static_assert(M > 0 && N > 0, "Matrix cannot have a zero dimension");
        size_t i = 0;
        for (auto it = list.begin(); it != list.end() && i < M; ++it, ++i) {
            size_t j = 0;
            for (auto sub_it = (*it).begin(); sub_it != (*it).end() && j < N;
                    ++sub_it, ++j) {
                v[i][j] = *sub_it;
            }
        }
    }

    // Column vector
    Matrix(const Vec<T, M>& other) {
        static_assert(M > 0, "Matrix cannot have a zero dimension");
        static_assert(N == 1, "Matrix does not share the vector's dimensions");
        for (size_t i = 0; i < M; ++i) {
            v[i][0] = other[i];
        }
    }

    // Row vector
    Matrix(const Vec<T, N, false>& other) {
        static_assert(N > 0, "Matrix cannot have a zero dimension");
        static_assert(M == 1, "Matrix does not share the vector's dimensions");
        for (size_t i = 0; i < N; ++i) {
            v[0][i] = other[i];
        }
    }

    ~Matrix() {
    }

    static Matrix zeros() {
        static_assert(M > 0 && N > 0, "Matrix cannot have a zero dimension");
        Matrix<T, M, N> rval;
        rval.zero();
        return rval;
    }

    static Matrix identity() {
        static_assert(M > 0 && N > 0, "Matrix cannot have a zero dimension");
        static_assert(M == N, "Identity matrix must be square");
        Matrix<T, M, N> rval;
        rval.zero();
        for (size_t i = 0; i < M; ++i) {
            rval[i][i] = T(1);
        }
        return rval;
    }

    void zero() {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                v[i][j] = T();
            }
        }
    }

    Vec<T, N> row(size_t m) const {
        Vec<T, N> rval;
        for (size_t i = 0; i < N; ++i) {
            rval[i] = v[m][i];
        }
        return rval;
    }

    Vec<T, M> col(size_t n) const {
        Vec<T, M> rval;
        for (size_t i = 0; i < M; ++i) {
            rval[i] = v[i][n];
        }
        return rval;
    }

    Matrix<T, N, M> transpose() const {
        Matrix<T, N, M> rval;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                rval[j][i] = v[i][j];
            }
        }
        return rval;
    }

    size_t height() const {
        return M;
    }

    size_t width() const {
        return N;
    }

    std::array<T, N>& operator[](size_t i) {
        return v[i];
    }

    const std::array<T, N>& operator[](size_t i) const {
        return v[i];
    }

    Matrix operator+(const Matrix<T, M, N>& rhs) const {
        Matrix<T, M, N> rval;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                rval[i][j] = v[i][j] + rhs[i][j];
            }
        }
        return rval;
    }

    Matrix operator-(const Matrix<T, M, N>& rhs) const {
        Matrix<T, M, N> rval;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                rval[i][j] = v[i][j] - rhs[i][j];
            }
        }
        return rval;
    }

    Matrix& operator+=(const Matrix<T, M, N>& rhs) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                v[i][j] += rhs[i][j];
            }
        }
        return *this;
    }

    Matrix& operator-=(const Matrix<T, M, N>& rhs) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                v[i][j] -= rhs[i][j];
            }
        }
        return *this;
    }

    Vec<T, M> operator*(const Vec<T, N>& rhs) const {
        Vec<T, M> rval;
        for (size_t i = 0; i < M; ++i) {
            rval[i] = T();
            for (size_t j = 0; j < N; ++j) {
                rval[i] += v[i][j] * rhs[j];
            }
        }
        return rval;
    }

    Matrix operator*(T factor) const {
        Matrix<T, M, N> rval;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                rval[i][j] = v[i][j] * factor;
            }
        }
        return rval;
    }

    template <size_t O>
    Matrix<T, M, O> operator*(const Matrix<T, N, O>& rhs) const {
        Matrix<T, M, O> rval;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < O; ++j) {
                rval[i][j] = T();
                for (size_t k = 0; k < N; ++k) {
                    rval[i][j] += v[i][k] * rhs[k][j];
                }
            }
        }
        return rval;
    }

    Matrix operator/(T factor) const {
        Matrix<T, M, N> rval;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                rval[i][j] = v[i][j] / factor;
            }
        }
        return rval;
    }

    Matrix& operator*=(T factor) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                v[i][j] *= factor;
            }
        }
        return *this;
    }
    
    Matrix& operator*=(const Matrix<T, M, N>& rhs) {
        static_assert(M == N, "Matrix dimensions must agree for in place "
                              "matrix multiply");
        return *this = *this * rhs;
    }

    Matrix& operator/=(T factor) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                v[i][j] /= factor;
            }
        }
        return *this;
    }

    // This is really slow for large matrices (and would also be very prone
    // to stack overflows)
    T determinant() const {
        static_assert(M == N, "Determinant defined only for square matrices");
        return _determinant(*this);
    }

    // A very slow method of solving a linear equation. It uses the determinant
    // to perform its calclations
    Matrix<T, N, 1> slow_solve(const Matrix<T, N, 1>& rhs) const {
        static_assert(M == N, "Only square matrices are solvable");
        Matrix<T, N, 1> rval;
        T d = determinant();
        for (size_t j = 0; j < N; ++j) {
            Matrix<T, N> temp = *this;
            for (size_t i = 0; i < N; ++i) {
                temp[i][j] = rhs[i][0];
            }
            rval[j][0] = temp.determinant() / d;
        }
        return rval;
    }

    // A very slow method of solving a linear equation. It uses the determinant
    // to perform its calclations
    Vec<T, N> slow_solve(const Vec<T, N>& rhs) const{
        static_assert(M == N, "Only square matrices are solvable");
        Vec<T, N> rval;
        T d = determinant();
        for (size_t j = 0; j < N; ++j) {
            Matrix<T, N> temp = *this;
            for (size_t i = 0; i < N; ++i) {
                temp[i][j] = rhs[i];
            }
            rval[j] = temp.determinant() / d;
        }
        return rval;
    }
};

// This is really slow for large matrices (and would also be very prone
// to stack overflows)
template <typename T, size_t N>
T _determinant(const Matrix<T, N>& mat) {
    T rval(0);
    bool sign = true;
    for (size_t skip = 0; skip < N; ++skip) {
        Matrix<T, N - 1> sub_mat;
        for (size_t i = 1; i < N; ++i) {
            for (size_t j = 0; j < skip; ++j) {
                sub_mat[i - 1][j] = mat[i][j];
            }
        }
        for (size_t i = 1; i < N; ++i) {
            for (size_t j = skip + 1; j < N; ++j) {
                sub_mat[i - 1][j - 1] = mat[i][j];
            }
        }
        if (sign) {
            rval += mat[0][skip] * _determinant(sub_mat);
        } else {
            rval -= mat[0][skip] * _determinant(sub_mat);
        }
        sign = !sign;
    }
    return rval;
}

// Specialize the above function so that it terminates
template <typename T>
T _determinant(const Matrix<T, 1>& mat) {
    return mat[0][0];
}

template <typename T, typename U, size_t M, size_t N>
Matrix<T, M, N> operator*(U f, const Matrix<T, M, N>& rhs) {
    T factor = T(f);
    Matrix<T, M, N> new_mat;
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            new_mat[i][j] = factor * rhs[i][j];
        }
    }
    return new_mat;
}

template <typename T, size_t M, size_t N>
std::ostream& operator<<(std::ostream& os, const Matrix<T, M, N>& rhs) {
    std::array<size_t, N> v;
    for (size_t j = 0; j < N; ++j) {
        v[j] = 0;
        for (size_t i = 0; i < M; ++i) {
            std::ostringstream s;
            s << +rhs[i][j];
            std::string str = s.str();
            if (str.size() > v[j]) {
                v[j] = str.size();
            }
        }
    }
    for (size_t i = 0; i < M; ++i) {
        os << "[ ";
        for (size_t j = 0; j < N; ++j) {
            os << std::setw(v[j]) << +rhs[i][j] << " ";
        }
        os << "]";
        if (i != M - 1) {
            os << std::endl;
        }
    }
    return os;
}

}

#endif
