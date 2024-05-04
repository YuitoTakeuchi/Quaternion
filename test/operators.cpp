#include "../Quaternion.hpp"
#include <cassert>
// Unary operators
void test_equals() {
    Quaternion p(1.0, 0.0, 0.0, 0.0);
    Quaternion q(0.0f, 1.0f, 0.0f, 0.0f);

    p = q;
}

void test_sums() {
    Quaterniond p(1.0, 0.0, 0.0, 0.0);
    Quaterniond q(1.0, 0.0, 0.0, 0.0);
    Quaternionf r(0.0f, 1.0f, 0.0f, 0.0f);

    Quaterniond s;
    // Quaternion同士
    s = p + q;
    assert(s.w == 2.0f && s.x == 0.0f && s.y == 0.0f && s.z == 0.0f);
    s = p + r;
    assert(s.w == 1.0f && s.x == 1.0f && s.y == 0.0f && s.z == 0.0f);

    s += p;
    assert(s.w == 2.0f && s.x == 1.0f && s.y == 0.0f && s.z == 0.0f);
    s += r;
    assert(s.w == 2.0f && s.x == 2.0f && s.y == 0.0f && s.z == 0.0f);

    // Quaternionとスカラー
    r = p + 1.0f;
    r = 1.0f + p;
}

void basic() {
    // test with wolfram alpha
    // (0+3i+4j+3k)*(4+3.9i-j-3k)
    Quaterniond p(0.0, 3.0, 4.0, 3.0);
    Quaterniond q(4.0f, 3.9f, -1.0f, -3.0f);

    auto r = p * q;

    // std::cout << r << std::endl;
    // std::cout << "Norm: " << r.norm() << std::endl;
    // std::cout << "Normalized: " << r.normalized() << std::endl;
    // std::cout << "Conjugate: " << r.conj() << std::endl;
    // std::cout << "Inverse: " << r.inverse() << std::endl;
    // Eigen::Matrix3d m = r.rotation_matrix<Eigen::Matrix3d>();
    // std::cout << "Rotation matrix:\n" << m << std::endl;
}

int main() {
    test_equals();
    test_sums();
    return 0;
}