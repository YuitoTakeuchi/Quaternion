#pragma once
#include <cmath>

template<typename T>
class Quaternion {
public:
    Quaternion() :  w(1), x(0), y(0), z(0) {}
    Quaternion(T w, T x, T y, T z) : w(w), x(x), y(y), z(z) {}
    Quaternion(const Quaternion& q) : w(q.w), x(q.x), y(q.y), z(q.z) {}
    template<typename Vector>
    Quaternion(const Vector& v) : w(v(0)), x(v(1)), y(v(2)), z(v(3)) {}

    template<typename U>
    Quaternion operator=(const Quaternion<U>& q) {
        w = q.w;
        x = q.x;
        y = q.y;
        z = q.z;
        return *this;
    }
    template<typename U>
    Quaternion operator=(const U& s) {
        w = s;
        x = 0;
        y = 0;
        z = 0;
        return *this;
    }

    template<typename U>
    T dot(const Quaternion<U>& q) const {
        return w * q.w + x * q.x + y * q.y + z * q.z;
    }
    Quaternion conj() const {
        return {w, -x, -y, -z};
    }
    Quaternion inverse() const {
        return conj() / dot(*this);
    }
    T norm() const {
        return sqrt(w * w + x * x + y * y + z * z);
    }
    Quaternion normalized() const {
        return *this / norm();
    }
    Quaternion normalize() {
        return *this /= norm();
    }

    template<typename U>
    Quaternion lerp(const Quaternion& q, U t) const {
        return *this + (q - *this) * t;
    }
    template<typename U>
    Quaternion slerp(const Quaternion& q, U t) const {
        T dot = this->dot(q);
        if (dot < 0.0) {
            dot = -dot;
            return lerp(-q, t);
        }
        if (dot > 0.95) {
            return lerp(q, t);
        }
        T angle = acos(dot);
        return (*this * sin(angle * (1 - t)) + q * sin(angle * t)) / sin(angle);
    }
    void to_euler(T& yaw, T& pitch, T& roll) const {
        roll  = atan2(2.0 * (y*z + w*x), 1.0-2.0*(x*x+y*y));
        pitch = asin (2.0 * (w*y - x*z));
        yaw   = atan2(2.0 * (x*y + w*z), 1.0-2.0*(y*y+z*z));
    }

    template<typename Mat>
    Mat rotation_matrix() const {
        Quaternion q = normalized();
        Mat m;
        m(0, 0) = q.w*q.w + q.x*q.x - q.y*q.y - q.z*q.z;
        m(0, 1) = 2.0*(q.x*q.y - q.w*q.z);
        m(0, 2) = 2.0*(q.x*q.z + q.w*q.y);
        m(1, 0) = 2.0*(q.x*q.y + q.w*q.z);
        m(1, 1) = q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z;
        m(1, 2) = 2.0*(q.y*q.z - q.w*q.x);
        m(2, 0) = 2.0*(q.x*q.z - q.w*q.y);
        m(2, 1) = 2.0*(q.y*q.z + q.w*q.x);
        m(2, 2) = q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z;
        return m;
    }
    template<typename Mat>
    void rotation_matrix(Mat& m) const {
        Quaternion q = normalized();
        m(0, 0) = q.w*q.w + q.x*q.x - q.y*q.y - q.z*q.z;
        m(0, 1) = 2.0*(q.x*q.y - q.w*q.z);
        m(0, 2) = 2.0*(q.x*q.z + q.w*q.y);
        m(1, 0) = 2.0*(q.x*q.y + q.w*q.z);
        m(1, 1) = q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z;
        m(1, 2) = 2.0*(q.y*q.z - q.w*q.x);
        m(2, 0) = 2.0*(q.x*q.z - q.w*q.y);
        m(2, 1) = 2.0*(q.y*q.z + q.w*q.x);
        m(2, 2) = q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z;
    }

    // w: real part, x, y, z: imaginary parts
    T w, x, y, z;

    Quaternion operator+() const {
        return *this;
    }
    Quaternion operator-() const {
        return {-w, -x, -y, -z};
    }
    Quaternion operator+=(const T& s) {
        w += s;
        return *this;
    }
    Quaternion operator -= (const T& s) {
        w -= s;
        return *this;
    }
    Quaternion operator *= (const T& s) {
        w *= s;
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }
    Quaternion operator /= (const T& s) {
        w /= s;
        x /= s;
        y /= s;
        z /= s;
        return *this;
    }
    template<typename U>
    Quaternion operator+=(const Quaternion<U>& q) {
        w += q.w;
        x += q.x;
        y += q.y;
        z += q.z;
        return *this;
    }
    template<typename U>
    Quaternion operator -= (const Quaternion<U>& q) {
        w -= q.w;
        x -= q.x;
        y -= q.y;
        z -= q.z;
        return *this;
    }
    template<typename U>
    Quaternion operator *= (const Quaternion<U>& q) {
        T wt = w * q.w - x * q.x - y * q.y - z * q.z;
        T xt = w * q.x + x * q.w + y * q.z - z * q.y;
        T yt = w * q.y - x * q.z + y * q.w + z * q.x;
        T zt = w * q.z + x * q.y - y * q.x + z * q.w;
        w = wt;
        x = xt;
        y = yt;
        z = zt;
        return *this;
    }
};

using Quaterniond = Quaternion<double>;
using Quaternionf = Quaternion<float>;

template<typename T, typename U>
inline Quaternion<T> operator+(const Quaternion<T>& q, const U& s) {
    return {q.w + s, q.x, q.y, q.z};
}

template<typename T, typename U>
inline Quaternion<T> operator-(const Quaternion<T>& q, const U& s) {
    return {q.w - s, q.x, q.y, q.z};
}

template<typename T, typename U>
inline Quaternion<T> operator*(const Quaternion<T>& q, const U& s) {
    return {q.w * s, q.x * s, q.y * s, q.z * s};
}

template<typename T, typename U>
inline Quaternion<T> operator/(const Quaternion<T>& q, const U& s) {
    return {q.w / s, q.x / s, q.y / s, q.z / s};
}

template<typename T, typename U>
inline Quaternion<T> operator+(const U& s, const Quaternion<T>& q) {
    return {q.w + s, q.x, q.y, q.z};
}

template<typename T, typename U>
inline Quaternion<T> operator-(const U& s, const Quaternion<T>& q) {
    return {q.w - s, q.x, q.y, q.z};
}

template<typename T, typename U>
inline Quaternion<T> operator*(const U& s, const Quaternion<T>& q) {
    return {q.w * s, q.x * s, q.y * s, q.z * s};
}

template<typename T, typename U>
inline Quaternion<T> operator/(const U& s, const Quaternion<T>& q) {
    return {q.w / s, q.x / s, q.y / s, q.z / s};
}

template<typename T, typename U>
inline Quaternion<T> operator+(const Quaternion<T>& q1, const Quaternion<U>& q2) {
    return {q1.w + q2.w, q1.x + q2.x, q1.y + q2.y, q1.z + q2.z};
}

template<typename T, typename U>
inline Quaternion<T> operator-(const Quaternion<T>& q1, const Quaternion<U>& q2) {
    return {q1.w - q2.w, q1.x - q2.x, q1.y - q2.y, q1.z - q2.z};
}

template<typename T, typename U>
inline Quaternion<T> operator*(const Quaternion<T>& q1, const Quaternion<U>& q2) {
    T w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
    T x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
    T y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
    T z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;
    return {w, x, y, z};
}

template<typename T, typename Stream>
inline Stream& operator<<(Stream& stream, const Quaternion<T>& q) {
    stream << q.w << " " << q.x << " " << q.y << " " << q.z;
    return stream;
}

template<typename T, typename U>
auto dot(const Quaternion<T>& q1, const Quaternion<U>& q2) -> decltype(q1.w * q2.w) {
    return q1.dot(q2);
}