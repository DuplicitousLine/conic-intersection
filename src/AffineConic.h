#ifndef AFFINE_CONIC_H
#define AFFINE_CONIC_H

#include <Eigen/Dense>

using Conic = Eigen::Matrix3f;

class AffineConic {
public:
    static Conic applyTransform(const Eigen::Matrix3f& T, const Eigen::Matrix3f& M);

    static Eigen::Matrix3f translate(const Eigen::Vector2f& t);
    static Eigen::Matrix3f rotate(float a);
    static Eigen::Matrix3f scale(const Eigen::Vector2f& s);

    static Conic empty();
    static Conic point(const Eigen::Vector2f& p);
    static Conic line(const Eigen::Vector2f& a, const Eigen::Vector2f& b);
    static Conic parallelLines(const Eigen::Vector2f& p, const Eigen::Vector2f& d, float s);
    static Conic intersectingLines(const Eigen::Vector2f& p1, const Eigen::Vector2f& q1, const Eigen::Vector2f& p2, const Eigen::Vector2f& q2);
    static Conic ellipse(const Eigen::Vector2f& center, float a, float b, float angle);
    static Conic hyperbola(const Eigen::Vector2f& center, float a, float b, float angle);
    static Conic parabola(const Eigen::Vector2f& vertex, float p, float angle);

    static Conic fromPoints(const std::vector<Eigen::Vector2f>& points);
    static Conic fromCoefficients(float A, float B, float C, float D, float E, float F);
};

#endif // AFFINE_CONIC_H
