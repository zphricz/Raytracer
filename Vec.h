#ifndef VEC_H
#define VEC_H

#include <cmath>
#include <array>
#include <iomanip>
#include <stdint.h>
#include <sstream>

/*template <typename T>
constexpr pi<T> = T(M_PI);*/

namespace Linear {

template <typename T, size_t M, size_t N>
class Matrix;

// Default value for C is column vector
template <typename T, size_t N, bool C = true>
class Vec {
private:
    std::array<T, N> v;

public:
    Vec() {
        static_assert(N > 0, "Vector cannot have a zero dimension");
    }

    template <typename U>
    Vec(const Vec<U, N, C>& other) {
        static_assert(N > 0, "Vector cannot have a zero dimension");
        for (size_t i = 0; i < N; ++i) {
            v[i] = other[i];
        }
    }

    template <typename U>
    Vec(std::initializer_list<U> list) {
        static_assert(N > 0, "Vector cannot have a zero dimension");
        size_t i = 0;
        for (auto it = list.begin(); it != list.end() && i < N; ++it, ++i) {
            v[i] = *it;
        }
    }

    template <size_t M>
    Vec(const Matrix<T, 1, M>& other) {
        static_assert(M > 0, "Vector cannot have a zero dimension");
        static_assert(N == M, "Mismatched vector/matrix length");
        static_assert(C == false, "Mismatched vector/matrix dimensions");
        for (size_t i = 0; i < M; ++i) {
            v[i] = other[0][i];
        }
    }

    template <size_t M>
    Vec(const Matrix<T, M, 1>& other) {
        static_assert(M > 0, "Vector cannot have a zero dimension");
        static_assert(N == M, "Mismatched vector/matrix length");
        static_assert(C == true, "Mismatched vector/matrix dimensions");
        for (size_t i = 0; i < M; ++i) {
            v[i] = other[i][0];
        }
    }

    Vec(const Matrix<T, 1, 1>& other) {
        static_assert(N == 1, "Mismatched vector/matrix length");
        v[0] = other[0][0];
    }

    static Vec zeros() {
        static_assert(N > 0, "Vector cannot have a zero dimension");
        Vec<T, N, C> rval;
        rval.zero();
        return rval;
    }

    ~Vec() {
    }

    auto magnitude() -> decltype(sqrt(T())) const {
        return sqrt(magnitude_square());
    }

    T magnitude_square() const {
        T sum = T();
        for (size_t i = 0; i < N; ++i) {
            sum += v[i] * v[i];
        }
        return sum;
    }

    void zero() {
        for (size_t i = 0; i < N; ++i) {
            v[i] = T();
        }
    }

    void normalize() {
        *this /= magnitude();
    }

    Vec& operator+=(const Vec<T, N, C>& rhs) {
        for (size_t i = 0; i < N; ++i) {
            v[i] += rhs[i];
        }
        return *this;
    }

    Vec& operator-=(const Vec<T, N, C>& rhs) {
        for (size_t i = 0; i < N; ++i) {
            v[i] -= rhs[i];
        }
        return *this;
    }

    Vec& operator*=(T factor) {
        for (size_t i = 0; i < N; ++i) {
            v[i] *= factor;
        }
        return *this;
    }

    Vec& operator/=(T factor) {
        for (size_t i = 0; i < N; ++i) {
            v[i] /= factor;
        }
        return *this;
    }

    Vec operator+(const Vec<T, N, C>& rhs) const {
        Vec<T, N, C> new_vec;
        for (size_t i = 0; i < N; ++i) {
            new_vec[i] = v[i] + rhs[i];
        }
        return new_vec;
    }

    Vec operator-(const Vec<T, N, C>& rhs) const {
        Vec<T, N, C> new_vec;
        for (size_t i = 0; i < N; ++i) {
            new_vec[i] = v[i] - rhs[i];
        }
        return new_vec;
    }

    Vec operator*(T factor) const {
        Vec<T, N, C> new_vec;
        for (size_t i = 0; i < N; ++i) {
            new_vec[i] = v[i] * factor;
        }
        return new_vec;
    }

    Vec operator/(T factor) const {
        Vec<T, N, C> new_vec;
        for (size_t i = 0; i < N; ++i) {
            new_vec[i] = v[i] / factor;
        }
        return new_vec;
    }

    T& operator[](size_t index) {
        return v[index];
    }

    T operator[](size_t index) const {
        return v[index];
    }

    bool operator==(const Vec<T, N, C>& rhs) const {
        return v == rhs.v;
    }

    bool operator!=(const Vec<T, N, C>& rhs) const {
        return v != rhs.v;
    }

    T dot(const Vec<T, N, C>& other) const {
        T sum = T();
        for (size_t i = 0; i < N; ++i) {
            sum += v[i] * other[i];
        }
        return sum;
    }

    size_t size() const {
        return N;
    }
};

template<typename T, size_t N, bool C>
static Vec<T, N, false> to_row(const Vec<T, N, C>& other) {
    Vec<T, N, false> rval;
    static_assert(N > 0, "Vector cannot have a zero dimension");
    for (size_t i = 0; i < N; ++i) {
        rval[i] = other[i];
    }
    return rval;
}

template<typename T, size_t N, bool C>
static Vec<T, N> to_col(const Vec<T, N, C>& other) {
    Vec<T, N> rval;
    static_assert(N > 0, "Vector cannot have a zero dimension");
    for (size_t i = 0; i < N; ++i) {
        rval[i] = other[i];
    }
    return rval;
}

/* Prints out the vector as a column vector
 */
template <typename T, size_t N, bool C>
std::ostream& operator<<(std::ostream& os, const Vec<T, N, C>& rhs) {
    if (!C) {
        os << "[ ";
    }
    size_t max_size = 0;
    for (size_t i = 0; i < N; ++i) {
        std::ostringstream s;
        s << +rhs[i];
        std::string str = s.str();
        if (str.size() > max_size) {
            max_size = str.size();
        }
    }
    for (size_t i = 0; i < N; ++i) {
        if (C) {
            os << "[ " << std::setw(max_size) << +rhs[i] << " ]";
            if (i != N - 1) {
                os << std::endl;
            }
        } else {
            os << rhs[i] << " ";
        }
    }
    if (!C) {
        os << "]";
    }
    return os;
}

template <typename T, typename U, size_t N, bool C>
Vec<T, N, C> operator*(U f, const Vec<T, N, C>& rhs) {
    T factor = T(f);
    Vec<T, N, C> new_vec;
    for (size_t i = 0; i < N; ++i) {
        new_vec[i] = factor * rhs[i];
    }
    return new_vec;
}


template <typename T, bool C>
class Vec<T, 2, C> {
public:
    T x;
    T y;

    Vec() {
    }

    template <typename U>
    Vec(const Vec<U, 2, C>& other) :
        x(other.x),
        y(other.y) {
    }

    template <typename U>
    Vec(std::initializer_list<U> list) {
        size_t i = 0;
        for (auto it = list.begin(); it != list.end() && i < 2; ++it, ++i) {
            (*this)[i] = *it;
        }
    }

    Vec(const Matrix<T, 1, 2>& other) :
        x(other[0][0]),
        y(other[0][1]) {
        static_assert(C == false, "Mismatched vector/matrix dimensions");
    }

    Vec(const Matrix<T, 2, 1>& other) :
        x(other[0][0]),
        y(other[1][0]) {
        static_assert(C == true, "Mismatched vector/matrix dimensions");
    }

    Vec(T angle) :
        x(cos(angle)),
        y(sin(angle)) {
    }

    Vec(T x, T y) :
        x(x),
        y(y) {
    }

    static Vec zeros() {
        Vec<T, 2, C> rval;
        rval.zero();
        return rval;
    }

    ~Vec() {
    }

    auto theta() -> decltype(atan2(T(), T())) const {
        return atan2(y, x);
    }

    auto magnitude() -> decltype(sqrt(T())) const {
        return sqrt(magnitude_square());
    }

    T magnitude_square() const {
        return x * x + y * y;
    }

    void zero() {
        x = T();
        y = T();
    }

    void normalize() {
        *this /= magnitude();
    }

    Vec& operator+=(const Vec<T, 2, C>& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    Vec& operator-=(const Vec<T, 2, C>& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    Vec& operator*=(T factor) {
        x *= factor;
        y *= factor;
        return *this;
    }

    Vec& operator/=(T factor) {
        x /= factor;
        y /= factor;
        return *this;
    }

    Vec operator+(const Vec<T, 2, C>& rhs) const {
        return Vec<T, 2, C>(x + rhs.x, y + rhs.y);
    }

    Vec operator-(const Vec<T, 2, C>& rhs) const {
        return Vec<T, 2, C>(x - rhs.x, y - rhs.y);
    }

    Vec operator*(T factor) const {
        return Vec<T, 2, C>(x * factor, y * factor);
    }

    Vec operator/(T factor) const {
        return Vec<T, 2, C>(x / factor, y / factor);
    }

    // Prefer accessing x and y directly to this
    T& operator[](size_t index) {
        switch (index) {
            case 0:
                return x;
            case 1:
                return y;
            default:
                return *(T*)nullptr; // Not quite safe
        }
    }

    // Prefer accessing x and y directly to this
    T operator[](size_t index) const {
        switch (index) {
            case 0:
                return x;
            case 1:
                return y;
            default:
                return *(volatile T*)nullptr; // Not quite safe
        }
    }

    bool operator==(const Vec<T, 2, C>& rhs) const {
        return x == rhs.x && y == rhs.y;
    }

    bool operator!=(const Vec<T, 2, C>& rhs) const {
        return !(*this == rhs);
    }

    T dot(const Vec<T, 2, C>& other) const {
        return x * other.x + y * other.y;
    }

    size_t size() const {
        return 2;
    }
};

typedef Vec<float, 2> Vec2f;
typedef Vec<double, 2> Vec2d;
typedef Vec<int8_t, 2> Vec2i8;
typedef Vec<int32_t, 2> Vec2i32;
typedef Vec<int64_t, 2> Vec2i64;
typedef Vec<uint8_t, 2> Vec2u8;
typedef Vec<uint32_t, 2> Vec2u32;
typedef Vec<uint64_t, 2> Vec2u64;

template <typename T, bool C>
class Vec<T, 3, C> {
public:
    T x;
    T y;
    T z;

    Vec() {
    }

    template <typename U>
    Vec(const Vec<U, 3, C>& other) :
        x(other.x),
        y(other.y),
        z(other.z) {
    }

    template <typename U>
    Vec(std::initializer_list<U> list) {
        size_t i = 0;
        for (auto it = list.begin(); it != list.end() && i < 3; ++it, ++i) {
            (*this)[i] = *it;
        }
    }

    Vec(const Matrix<T, 1, 3>& other) :
        x(other[0][0]),
        y(other[0][1]),
        z(other[0][2]) {
        static_assert(C == false, "Mismatched vector/matrix dimensions");
    }

    Vec(const Matrix<T, 3, 1>& other) :
        x(other[0][0]),
        y(other[1][0]),
        z(other[2][0]) {
        static_assert(C == true, "Mismatched vector/matrix dimensions");
    }

    /*
     * Pitch needs to be bounded by [-PI/2, PI/2]
     */
    Vec(T _yaw, T _pitch) :
        x(cos(_pitch) * sin(_yaw)),
        y(sin(_pitch)),
        z(cos(_pitch) * cos(_yaw)) {
    }

    Vec(T x, T y, T z) :
        x(x),
        y(y),
        z(z) {
    }

    static Vec zeros() {
        Vec<T, 3, C> rval;
        rval.zero();
        return rval;
    }

    ~Vec() {
    }

    auto magnitude() -> decltype(sqrt(T())) const {
        return sqrt(magnitude_square());
    }

    T magnitude_square() const {
        return x * x + y * y + z * z;
    }

    void zero() {
        x = T();
        y = T();
        z = T();
    }

    void normalize() {
        *this /= magnitude();
    }

    Vec& operator+=(const Vec<T, 3, C>& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    Vec& operator-=(const Vec<T, 3, C>& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    Vec& operator*=(T factor) {
        x *= factor;
        y *= factor;
        z *= factor;
        return *this;
    }

    Vec& operator/=(T factor) {
        x /= factor;
        y /= factor;
        z /= factor;
        return *this;
    }

    Vec operator+(const Vec<T, 3, C>& rhs) const {
        return Vec<T, 3, C>(x + rhs.x, y + rhs.y, z + rhs.z);
    }

    Vec operator-(const Vec<T, 3, C>& rhs) const {
        return Vec<T, 3, C>(x - rhs.x, y - rhs.y, z - rhs.z);
    }

    Vec operator*(T factor) const {
        return Vec<T, 3, C>(x * factor, y * factor, z * factor);
    }

    Vec operator/(T factor) const {
        return Vec<T, 3, C>(x / factor, y / factor, z / factor);
    }

    // Prefer accessing x, y, and z directly to this
    T& operator[](size_t index) {
        switch (index) {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            default:
                return *(T*)nullptr; // Not quite safe
        }
    }

    // Prefer accessing x, y, and z directly to this
    T operator[](size_t index) const {
        switch (index) {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            default:
                return *(volatile T*)nullptr; // Not quite safe
        }
    }

    bool operator==(const Vec<T, 3, C>& rhs) const {
        return x == rhs.x && y == rhs.y && z == rhs.z;
    }

    bool operator!=(const Vec<T, 3, C>& rhs) const {
        return !(*this == rhs);
    }

    T dot(const Vec<T, 3, C>& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    size_t size() const {
        return 3;
    }


    /* ^
     * |
     * |
     * |
     * Y  (changes in height along Y)
     * |
     * |
     * |
     * v
     *  <--------X-------->
     */
    
    // A vector laying entirely within the xz plane has a pitch of 0.0
    // As the vector points towards the positive y axis, it's pitch increases
    // As the vector points towards the negative y axis, it's pitch decreases
    auto pitch() -> decltype(atan2(sqrt(T() * T() + T() * T()), T())) const {
        return atan2(y, sqrt(x * x + z * z));
    }

    // A vector pointing entirely along the positive z axis has a yaw of 0.0
    // As the vector points towards the positive x axis, it's yaw increases
    // As the vector points towards the negative x axis, it's yaw decreases
    auto yaw() -> decltype(atan2(T(), T())) const {
        return atan2(x, z);
    }

    // A positive delta will rotate the vector clockwise on the xz plane when
    // viewed from the positive y axis
    void rotate_yaw(T delta) {
#if 1
        auto temp = x * cos(delta) + z * sin(delta);
        z = z * cos(delta) - x * sin(delta);
        x = temp;
#else
        Matrix<T, 3, 3> lhs{{cos(delta),  T(), sin(delta)},
                            {T(),         T(1),       T()},
                            {-sin(delta), T(), cos(delta)}};

        // Figure out why the above doesn't work, but this does
        /* Matrix<T, 3, 3> lhs = {{T(cos(delta)),  T(), T(sin(delta))},
                            {T(),         T(1),       T()},
                            {T(-sin(delta)), T(), T(cos(delta))}};*/
        *this = lhs * *this;
#endif
    }

    // A positive delta will rotate the vector towards the positive y axis
    // A negative delta will rotate the vector towards the negative y axis
    // The result of this function should not make pitch leave the range
    //   of [-PI/2, PI/2]
    void rotate_pitch(T delta) {
        auto v = sqrt(x * x + z * z);
        auto temp = y * cos(delta) + v * sin(delta);
        auto new_v = v * cos(delta) - y * sin(delta);
        auto factor = new_v / v;
        x *= factor;
        z *= factor;
        y = temp;
    }

    void set_yaw(T _yaw) {
        rotate_yaw(_yaw - yaw());
    }

    // _pitch should be in the range of  [-PI/2, PI/2]
    void set_pitch(T _pitch) {
        rotate_pitch(_pitch - pitch());
    }

    Vec cross(const Vec<T, 3, C>& other) const {
        return Vec<T, 3, C>(y * other.z - z * other.y,
                         z * other.x - x * other.z,
                         x * other.y - y * other.x);
    }
};


typedef Vec<float, 3> Vec3f;
typedef Vec<double, 3> Vec3d;
typedef Vec<int8_t, 3> Vec3i8;
typedef Vec<int32_t, 3> Vec3i32;
typedef Vec<int64_t, 3> Vec3i64;
typedef Vec<uint8_t, 3> Vec3u8;
typedef Vec<uint32_t, 3> Vec3u32;
typedef Vec<uint64_t, 3> Vec3u64;

}

#endif

