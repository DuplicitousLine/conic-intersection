#include "ConicIntersection.h"

#include <unsupported/Eigen/Polynomials>
#include <cassert>
#include <iostream>

constexpr float EPSILON = 1e-6f;

/*
 * Computes the intersection of two conics (E1 and E2).
 * It handles the various degenerate cases based on the rank of the conics.
 *
 * Note on rank interpretations:
 * - Rank 0: A degenerate conic representing an invalid or degenerate case (zero matrix).
 * - Rank 1: A degenerate conic representing a double line (l^T * l).
 * - Rank 2 with indefinite quadratic: Two distinct real lines.
 * - Rank 2 with definite quadratic: A single point (complex conjugate lines).
 * - Rank 3: A proper (non-degenerate) conic.
 */
std::vector<Point> ConicIntersection::intersect(const Conic& E1, const Conic& E2) {
    assert(E2.cols() == 3 && E2.rows() == 3 && E1.cols() == 3 && E1.rows() == 3);
    assert(E1(0, 1) == E1(1, 0) && E1(0, 2) == E1(2, 0) && E1(1, 2) == E1(2, 1));
    assert(E2(0, 1) == E2(1, 0) && E2(0, 2) == E2(2, 0) && E2(1, 2) == E2(2, 1));

    Eigen::FullPivLU<Eigen::Matrix3f> lu1(E1);
    Eigen::FullPivLU<Eigen::Matrix3f> lu2(E2);

    int r1 = lu1.rank();
    int r2 = lu2.rank();

    if (r1 == 3 && r2 == 3) {
        return completeIntersection(E1, E2);
    }
    else if (r1 == 3 && r2 < 3) {
        return intersectFullWithDegenerate(E1, E2, r2);
    }
    else if (r2 == 3 && r1 < 3) {
        return intersectFullWithDegenerate(E2, E1, r1);
    }

    return intersectBothDegenerate(E1, r1, E2, r2);
}

/*
 * Compute the intersection of two full-rank conics.
 * Solves the generalized eigenvalue problem to find the linear combination
 * that produces a degenerate conic, then intersects the full conic with it.
 */
std::vector<Point> ConicIntersection::completeIntersection(const Conic& E1, const Conic& E2) {
    Eigen::Matrix3f EE = E1 * (-E2).inverse();
    float detEE12 = EE(0, 0) * EE(1, 1) - EE(0, 1) * EE(1, 0);
    float detEE23 = EE(1, 1) * EE(2, 2) - EE(1, 2) * EE(2, 1);
    float detEE13 = EE(0, 0) * EE(2, 2) - EE(0, 2) * EE(2, 0);

    Eigen::Vector4f k;
    k << EE.determinant(), -(detEE12 + detEE23 + detEE13), EE.trace(), -1;

    Eigen::PolynomialSolver<float, 3> solver;
    solver.compute(k);
    bool realSolutions;
    float largestRealRoots = solver.greatestRealRoot(realSolutions);

    std::vector<Point> intersections;
    if (!realSolutions) {
        return intersections;
    }
    else {
        Conic E0 = E1 + largestRealRoots * E2;

        Eigen::FullPivLU<Eigen::Matrix3f> lu(E0);
        lu.setThreshold(EPSILON);

        int r0 = lu.rank();

        return intersectFullWithDegenerate(E1, E0, 2);
    }
}

/*
 * Intersect a full-rank conic with a degenerate conic.
 */
std::vector<Point> ConicIntersection::intersectFullWithDegenerate(
    const Conic& fullE, const Conic& degenerateE, int rankD) {

    assert(rankD < 3);

    std::vector<Point> intersections;

    if (rankD == 2 && isPointConic(degenerateE)) {
        Point p = recoverPoint(degenerateE);

        if (std::abs(p.transpose() * fullE * p) < EPSILON) {
            intersections.push_back(p);
        }
    }
    else {
        std::vector<Line> lines = decomposeDegenerateConic(degenerateE, rankD);

        for (const auto& line : lines) {
            std::vector<Point> P = intersectLine(fullE, line);
            intersections.insert(intersections.end(), P.begin(), P.end());
        }
    }

    return intersections;
}

/*
 * Intersect two degenerate conics (rank < 3).
 */
std::vector<Point> ConicIntersection::intersectBothDegenerate(
    const Conic& E1, int rank1,
    const Conic& E2, int rank2) {

    assert(rank1 < 3 && rank2 < 3);

    std::vector<Point> intersections;

    bool isPoint1 = (rank1 == 2) && isPointConic(E1);
    bool isPoint2 = (rank2 == 2) && isPointConic(E2);

    if (isPoint1 && isPoint2) {
        Point p1 = recoverPoint(E1);
        Point p2 = recoverPoint(E2);

        if (p1.isApprox(p2, EPSILON)) {
            intersections.push_back(p1);
        }
    }

    if (isPoint1 || isPoint2) {
        const Conic& pointConic = isPoint1 ? E1 : E2;
        const Conic& lineConic = isPoint1 ? E2 : E1;
        int lineRank = isPoint1 ? rank2 : rank1;

        Point p = recoverPoint(pointConic);

        std::vector<Line> lines = decomposeDegenerateConic(lineConic, lineRank);

        for (const Line& l : lines) {
            if (std::abs(l.dot(p)) < EPSILON) {
                intersections.push_back(p);
                break;
            }
        }
    }

    if (!isPoint1 && !isPoint2) {
        std::vector<Line> lines1 = decomposeDegenerateConic(E1, rank1);
        std::vector<Line> lines2 = decomposeDegenerateConic(E2, rank2);

        for (const Line& l1 : lines1) {
            for (const Line& l2 : lines2) {
                Point p;
                p << l1(1) * l2(2) - l1(2) * l2(1),
                     l1(2) * l2(0) - l1(0) * l2(2),
                     l1(0) * l2(1) - l1(1) * l2(0);

                if (std::abs(p(2)) > EPSILON) {
                    intersections.push_back(p);
                }
            }
        }
    }

    return intersections;
}

/*
 * Decompose a degenerate conic into one or two lines.
 *
 * - Rank 1: the conic represents a single line (l^T * l).
 * - Rank 2: the conic represents two lines (l^T * m + m^T * l)
 */
std::vector<Line> ConicIntersection::decomposeDegenerateConic(
    const Conic& de, int dRank) {

    assert(std::abs(de.determinant()) < EPSILON);

    Eigen::Matrix2f A = de.topLeftCorner<2, 2>();
    return A.isZero();

    bool singleLine = isSingleLineConic(de);

    Eigen::Matrix3f C;
    if (dRank == 1 && !singleLine) {
        C = de;
    }
    else {
        Point p = recoverPoint(de);
        C = de + crossMatrix(p);
    }

    float currentMax = 0;
    int maxI, maxJ;

    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
            if (currentMax < std::abs(C(i, j))) {
                maxI = i; maxJ = j;
                currentMax = std::abs(C(i, j));
            }
        }
    }

    std::vector<Line> lines;

    if (dRank == 2 || !singleLine) {
        lines.push_back(C.row(maxI));
    }

    if (dRank == 2 || singleLine) {
        lines.push_back(C.col(maxJ).transpose());
    }

    return lines;
}

/*
 * Intersect a full conic with a line.
 * Solves the quadratic equation along the line to find intersection points.
 */
std::vector<Point> ConicIntersection::intersectLine(const Conic& C, const Line& l) {
    assert(C.cols() == 3 && C.rows() == 3 && l.cols() == 3 && l.rows() == 1);

    std::vector<Point> points;
    Point p1, p2;
    float k1, k2;

    if (l(0) == 0 && l(1) == 0) return points;

    std::tie(p1, p2) = pointsOnLine(l);

    float p1Cp1 = p1.transpose() * C * p1;
    float p2Cp2 = p2.transpose() * C * p2;
    float p1Cp2 = p1.transpose() * C * p2;

    if (p2Cp2 == 0) {
        k1 = -0.5f * p1Cp1 / p1Cp2;
        points.push_back(p1 + k1 * p2);
    }
    else {
        float delta = p1Cp2 * p1Cp2 - p1Cp1 * p2Cp2;

        if (delta >= 0) {
            float deltaSqrt = std::sqrt(delta);
            k1 = (-p1Cp2 + deltaSqrt) / p2Cp2;
            k2 = (-p1Cp2 - deltaSqrt) / p2Cp2;

            points.push_back(p1 + k1 * p2);
            points.push_back(p1 + k2 * p2);
        }
    }

    return points;
}

/*
 * For a rank-2 conic to represent a point (rather than two real lines),
 * the quadratic part (the upper-left 2x2 block) must be definite.
 */
bool ConicIntersection::isPointConic(const Conic& C) {
    Eigen::Matrix2f A = C.topLeftCorner<2, 2>();

    float a11 = A(0, 0);
    float det = A.determinant();

    bool positiveDefinite = (a11 > EPSILON && det > EPSILON);
    bool negativeDefinite = (a11 < -EPSILON && det > EPSILON);

    return positiveDefinite || negativeDefinite;
}

bool ConicIntersection::isSingleLineConic(const Conic& C) const {
    bool quadraticZero = std::abs(C(0, 0)) < EPSILON &&
                         std::abs(C(1, 1)) < EPSILON &&
                         std::abs(C(0, 1)) < EPSILON;
    bool linearNonZero = std::abs(C(0, 2)) > EPSILON ||
                         std::abs(C(1, 2)) > EPSILON;
    return quadraticZero && linearNonZero;
}

/*
 * Recover the intersection point of a degenerate conic (two lines).
 *
 * To localize the intersection point p, we use the dual conic B associated with
 * the conic C. The dual conic represents the locus of all tangent lines to C,
 * and it can be expressed as a quadratic form. For a degenerate conic, the dual
 * conic consists of a double point, which represents the intersection of the
 * original lines l and m.
 *
 * The intersection point is recovered by:
 * 1. Computing the adjoint matrix B of the conic C.
 * 2. Finding the column B_j corresponding to the largest diagonal element B_ij.
 * 3. Normalizing the column to recover the intersection point p = B_j / sqrt(B_ij).
 */
Point ConicIntersection::recoverPoint(const Conic& C) {
    Eigen::Matrix3f B = adjointSym(C);

    Eigen::Vector3f diagonal = B.diagonal().array().abs();
    float v = 0;
    int idx = 0;
    for (int i = 0; i < 3; i++) {
        if (v < diagonal(i)) {
            v = diagonal(i);
            idx = i;
        }
    }

    assert(std::abs(B(idx, idx)) > 0);

    float b = std::sqrt(v);
    Eigen::Vector3f p = B.col(idx) / b;

    return p;
}

/*
 * Generate two points lying on a given homogeneous line `l`.
 *
 * - p1: a finite point on the line, computed by setting x=0 or y=0 and solving the line equation.
 * - p2: a point at infinity in the direction of the line, used to represent the line's direction in homogeneous coordinates.
 *
 * Together, these two points can be used to represent or parametrize the line.
 */
std::pair<Point, Point> ConicIntersection::pointsOnLine(const Line& l) {
    Point p1, p2;
    assert(l(0) != 0 || l(1) != 0);

    p2 << -l(1), l(0), 0;
    if (std::abs(l(0)) < std::abs(l(1))) {
        p1 << 0, -l(2), l(1);
    }
    else {
        p1 << -l(2), 0, l(0);
    }

    return std::make_pair(p1, p2);
}

/*
 * Return the cross-product matrix for a point.
 */
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

/*
 * Compute the adjoint of a symmetric 3x3 matrix.
 */
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
