#ifndef CONIC_INTERSECTION_H
#define CONIC_INTERSECTION_H

#include <vector>
#include <tuple>

#include <Eigen/Core>

using Point = Eigen::Vector3f;
using Line = Eigen::Matrix<float, 1, 3>;
using Conic = Eigen::Matrix3f;

class ConicIntersection {
public:
    ConicIntersection() = default;
    ~ConicIntersection() = default;

    std::vector<Point> intersect(const Conic& a, const Conic& b);

private:
    std::vector<Point> completeIntersection(const Conic& e1, const Conic& e2);
    std::vector<Point> degenerateIntersection(const Conic& fullE, const Conic& degenerateE);
    std::vector<Point> degenerateIntersection(const Conic& fullE, const Conic& degenerateE, int degenerateRank);
    std::vector<Point> intersectLine(const Conic& c, const Line& l);

    std::tuple<bool, Line, Line> decomposeDegenerateConic(const Conic& de);
    std::tuple<bool, Line, Line> decomposeDegenerateConic(const Conic& de, int rank);

    std::pair<Point, Point> pointsOnLine(const Line& l);

    Eigen::Matrix3f crossMatrix(const Point& p);
    Eigen::Matrix3f adjointSym(const Eigen::Matrix3f& M);
};

#endif // CONIC_INTERSECTION_H
