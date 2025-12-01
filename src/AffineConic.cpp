#include "AffineConic.h"

Eigen::Matrix3f AffineConic::translate(const Eigen::Vector2f& t) {
    Eigen::Matrix3f T = Eigen::Matrix3f::Identity();
    T(0, 2) = t.x();
    T(1, 2) = t.y();
    return T;
}

Eigen::Matrix3f AffineConic::rotate(float a) {
    float c = std::cos(a);
    float s = std::sin(a);

    Eigen::Matrix3f R;
    R << c, -s, 0.0f,
        s, c, 0.0f,
        0.0f, 0.0f, 1.0f;

    return R;
}

Eigen::Matrix3f AffineConic::scale(const Eigen::Vector2f& s) {
    Eigen::Matrix3f S = Eigen::Matrix3f::Identity();
    S(0, 0) = s.x();
    S(1, 1) = s.y();
    return S;
}

Eigen::Matrix3f AffineConic::getEllipseMatrix(const Eigen::Vector2f& t, float a, const Eigen::Vector2f& s) {
    Eigen::Matrix3f T = translate(t) * rotate(a) * scale(s);
    Eigen::Matrix3f D = Eigen::Vector3f(1.0f, 1.0f, -1.0f).asDiagonal();
    return T * D * T.transpose();
}

Eigen::Matrix3f AffineConic::getParabolaMatrix(const Eigen::Vector2f& t, float a, const Eigen::Vector2f& s) {
    const float alpha = 1.0f / std::sqrt(2.0f);
    Eigen::Matrix3f M;
    M << 1.0f, 0.0f, 0.0f,
        0.0f, alpha, -alpha,
        0.0f, alpha, alpha;

    Eigen::Matrix3f T = translate(t) * rotate(a) * scale(s) * M;
    Eigen::Matrix3f D = Eigen::Vector3f(1.0f, -0.5f, 0.5f).asDiagonal();
    return T * D * T.transpose();
}

Eigen::Matrix3f AffineConic::getHyperbolaMatrix(const Eigen::Vector2f& t, float a, const Eigen::Vector2f& s) {
    Eigen::Matrix3f T = translate(t) * rotate(a) * scale(s);
    Eigen::Matrix3f D = Eigen::Vector3f(1.0f, -1.0f, -1.0f).asDiagonal();
    return T * D * T.transpose();
}
