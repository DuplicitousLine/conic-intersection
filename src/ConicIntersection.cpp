#include "ConicIntersection.h"

#include <unsupported/Eigen/Polynomials>
#include <cassert>

std::vector<Point> ConicIntersection::intersect(const Conic& E1, const Conic& E2) {
    assert(E2.cols() == 3 && E2.rows() == 3 && E1.cols() == 3 && E1.rows() == 3);
    assert(E1(0, 1) == E1(1, 0) && E1(0, 2) == E1(2, 0) && E1(1, 2) == E1(2, 1));
    assert(E2(0, 1) == E2(1, 0) && E2(0, 2) == E2(2, 0) && E2(1, 2) == E2(2, 1));

    Eigen::FullPivLU<Eigen::Matrix3f> lu1(E1);
    Eigen::FullPivLU<Eigen::Matrix3f> lu2(E2);

    auto r1 = lu1.rank();
    auto r2 = lu2.rank();

    if (r1 == 3 && r2 == 3) {
        return completeIntersection(E1, E2);
    }
    else if (r2 < 3) {
        return degenerateIntersection(E1, E2, r2);
    }
    else {
        return degenerateIntersection(E2, E1, r1);
    }
}

std::vector<Point> ConicIntersection::completeIntersection(const Conic& E1, const Conic& E2) {
    Eigen::Matrix3f EE = E1 * (-E2).inverse();
    auto detEE12 = EE(0, 0) * EE(1, 1) - EE(0, 1) * EE(1, 0);
    auto detEE23 = EE(1, 1) * EE(2, 2) - EE(1, 2) * EE(2, 1);
    auto detEE13 = EE(0, 0) * EE(2, 2) - EE(0, 2) * EE(2, 0);

    Eigen::Vector4f k;
    k << EE.determinant(), -(detEE12 + detEE23 + detEE13), EE.trace(), -1;

    Eigen::PolynomialSolver<float, 3> solver;
    solver.compute(k);
    bool realSolutions;
    auto largestRealRoots = solver.greatestRealRoot(realSolutions);

    std::vector<Point> intersections;
    if (!realSolutions) {
        return intersections;
    }
    else {
        auto E0 = E1 + largestRealRoots * E2;
        Conic E0Conic(E0);
        auto output = decomposeDegenerateConic(E0Conic);

        if (!std::get<0>(output)) {
            return intersections;
        }

        Line m = std::get<1>(output);
        Line l = std::get<2>(output);

        auto P1 = intersectLine(E1, m);
        auto P2 = intersectLine(E1, l);

        intersections.insert(intersections.end(), P1.begin(), P1.end());
        intersections.insert(intersections.end(), P2.begin(), P2.end());

        return intersections;
    }
}

std::vector<Point> ConicIntersection::degenerateIntersection(const Conic& fullE, const Conic& degenerateE) {
    Eigen::FullPivLU<Eigen::Matrix3f> lu(degenerateE);
    return degenerateIntersection(fullE, degenerateE, lu.rank());
}

std::vector<Point> ConicIntersection::degenerateIntersection(const Conic& fullE, const Conic& degenerateE, int rankD) {
    auto output = decomposeDegenerateConic(degenerateE, rankD);
    std::vector<Point> intersections;
    if (!std::get<0>(output)) {
        return intersections;
    }

    Line m = std::get<1>(output);
    Line l = std::get<2>(output);

    auto P1 = intersectLine(fullE, m);
    auto P2 = intersectLine(fullE, l);

    intersections.insert(intersections.end(), P1.begin(), P1.end());
    intersections.insert(intersections.end(), P2.begin(), P2.end());

    return intersections;
}

std::vector<Point> ConicIntersection::intersectLine(const Conic& C, const Line& l) {
    assert(C.cols() == 3 && C.rows() == 3 && l.cols() == 3 && l.rows() == 1);

    std::vector<Point> pts;
    Point p1, p2;
    float k1, k2;

    if (l(0) == 0 && l(1) == 0) return pts;

    std::tie(p1, p2) = pointsOnLine(l);

    float p1Cp1 = p1.transpose() * C * p1;
    float p2Cp2 = p2.transpose() * C * p2;
    float p1Cp2 = p1.transpose() * C * p2;

    if (p2Cp2 == 0) {
        k1 = -0.5f * p1Cp1 / p1Cp2;
        pts.push_back(p1 + k1 * p2);
    }
    else {
        float delta = p1Cp2 * p1Cp2 - p1Cp1 * p2Cp2;
        if (delta >= 0) {
            float deltaSqrt = std::sqrt(delta);
            k1 = (-p1Cp2 + deltaSqrt) / p2Cp2;
            k2 = (-p1Cp2 - deltaSqrt) / p2Cp2;

            pts.push_back(p1 + k1 * p2);
            pts.push_back(p1 + k2 * p2);
        }
    }
    return pts;
}

std::tuple<bool, Line, Line> ConicIntersection::decomposeDegenerateConic(const Conic& de) {
    Eigen::FullPivLU<Eigen::Matrix3f> lu(de);
    return decomposeDegenerateConic(de, lu.rank());
}

std::tuple<bool, Line, Line> ConicIntersection::decomposeDegenerateConic(const Conic& de, int dRank) {
    assert(std::abs(de.determinant()) < 0.0001f);

    Eigen::Matrix3f C;
    if (dRank == 1) {
        C = de;
    }
    else {
        auto B = -adjointSym(de);
        Eigen::Vector3f diagonal = B.diagonal().array().abs();
        float v = diagonal.maxCoeff();
        int idx;
        diagonal.maxCoeff(&idx);

        if (B(idx, idx) < 0) {
            return std::make_tuple(false, Line(), Line());
        }

        float b = std::sqrt(v);
        auto p = B.col(idx) / b;
        auto Mp = crossMatrix(p);
        C = de + Mp;
    }

    float currentMax = 0;
    int maxI = 0, maxJ = 0;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            if (currentMax < std::abs(C(i, j))) {
                maxI = i;
                maxJ = j;
                currentMax = std::abs(C(i, j));
            }
        }
    }

    Line l = C.row(maxI);
    Line m = C.col(maxJ).transpose();

    return std::make_tuple(true, l, m);
}

std::pair<Point, Point> ConicIntersection::pointsOnLine(const Line& l) {
    Point p1, p2;
    assert(l(0) != 0 || l(1) != 0);

    p2 << -l(1), l(0), 0;
    if (std::abs(l(0)) < std::abs(l(1))) {
        p1 << 0, -l(2), l(1);
    } else {
        p1 << -l(2), 0, l(0);
    }
    return std::make_pair(p1, p2);
}

Eigen::Matrix3f ConicIntersection::crossMatrix(const Point& p) {
    Eigen::Matrix3f Mp = Eigen::Matrix3f::Zero();
    Mp(0, 1) = p(2);
    Mp(0, 2) = -p(1);
    Mp(1, 0) = -p(2);
    Mp(1, 2) = p(0);
    Mp(2, 0) = p(1);
    Mp(2, 1) = -p(0);

    return Mp;
}

Eigen::Matrix3f ConicIntersection::adjointSym(const Eigen::Matrix3f& M) {
    Eigen::Matrix3f A;
    auto a = M(0, 0), b = M(0, 1), d = M(0, 2);
    auto c = M(1, 1), e = M(1, 2);
    auto f = M(2, 2);

    A(0, 0) = c * f - e * e;
    A(0, 1) = -b * f + e * d;
    A(0, 2) = b * e - c * d;

    A(1, 0) = A(0, 1);
    A(1, 1) = a * f - d * d;
    A(1, 2) = -a * e + b * d;

    A(2, 0) = A(0, 2);
    A(2, 1) = A(1, 2);
    A(2, 2) = a * c - b * b;

    return A;
}
