#ifndef CONIC_INTERSECTION_H
#define CONIC_INTERSECTION_H

#include <vector>
#include <optional>

#include <Eigen/Core>

using Conic = Eigen::Matrix3f;
using Point = Eigen::Vector3f;
using Line = Eigen::RowVector3f;

class ConicIntersection {
public:
    ConicIntersection() = default;
    ~ConicIntersection() = default;

    std::vector<Point> intersect(const Conic& E1, const Conic& E2);

private:
    std::vector<Point> completeIntersection(const Conic& E1, const Conic& E2);

    std::vector<Point> intersectFullWithDegenerate(
        const Conic& fullE, const Conic& degenerateE, int rankD);

    std::vector<Point> intersectBothDegenerate(
        const Conic& E1, int rank1,
        const Conic& E2, int rank2);

    std::vector<Line> decomposeDegenerateConic(
        const Conic& de, int dRank);

    std::vector<Point> intersectLine(const Conic& C, const Line& l);

    bool isSingleLineConic(const Conic& C);
    bool isPointConic(const Conic& C);
    Point recoverPoint(const Conic& C);

    std::pair<Point, Point> pointsOnLine(const Line& l);
    Eigen::Matrix3f crossMatrix(const Point& p);
    Eigen::Matrix3f adjointSym(const Eigen::Matrix3f& M);
};

#endif // CONIC_INTERSECTION_H
