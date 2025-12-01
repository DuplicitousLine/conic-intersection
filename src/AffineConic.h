#ifndef AFFINE_CONIC_H
#define AFFINE_CONIC_H

#include <Eigen/Core>
#include <Eigen/Geometry>

class AffineConic {
public:
    static Eigen::Matrix3f translate(const Eigen::Vector2f& t);
    static Eigen::Matrix3f rotate(float a);
    static Eigen::Matrix3f scale(const Eigen::Vector2f& s);
    static Eigen::Matrix3f getEllipseMatrix(const Eigen::Vector2f& t, float a, const Eigen::Vector2f& s);
    static Eigen::Matrix3f getParabolaMatrix(const Eigen::Vector2f& t, float a, const Eigen::Vector2f& s);
    static Eigen::Matrix3f getHyperbolaMatrix(const Eigen::Vector2f& t, float a, const Eigen::Vector2f& s);
};

#endif // AFFINE_CONIC_H
