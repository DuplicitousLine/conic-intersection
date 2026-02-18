#include "AffineConic.h"

Conic AffineConic::applyTransform(const Eigen::Matrix3f& T, const Eigen::Matrix3f& M) {
    Eigen::Matrix3f invT = T.inverse();
    Eigen::Matrix3f C = invT.transpose() * M * invT;
    return Conic(C);
}

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

Conic AffineConic::empty() {
    return Eigen::Matrix3f::Identity();
}

Conic AffineConic::point(const Eigen::Vector2f& p) {
    Eigen::Matrix3f Q;
    Q << 1.0f, 0.0f, -p.x(),
        0.0f, 1.0f, -p.y(),
        -p.x(), -p.y(), p.x() * p.x() + p.y() * p.y();
    return Conic(Q);
}

Conic AffineConic::line(const Eigen::Vector2f& a, const Eigen::Vector2f& b) {
    float A = b.y() - a.y();
    float B = a.x() - b.x();
    float C = b.x() * a.y() - a.x() * b.y();

    Eigen::RowVector3f l(A, B, C);
    return l.transpose() * l;
}

Conic AffineConic::parallelLines(
    const Eigen::Vector2f& p,
    const Eigen::Vector2f& d,
    float s
) {
    Eigen::Vector2f dir = d.normalized();
    Eigen::Vector2f n(-dir.y(), dir.x());

    Eigen::Vector2f p1 = p + 0.5f * s * n;
    Eigen::Vector2f p2 = p - 0.5f * s * n;

    Eigen::Vector2f q1 = p1 + dir;
    Eigen::Vector2f q2 = p2 + dir;

    Eigen::RowVector3f l1;
    l1 << (p1.y() - q1.y()),
        (q1.x() - p1.x()),
        (p1.x() * q1.y() - q1.x() * p1.y());

    Eigen::RowVector3f l2;
    l2 << (p2.y() - q2.y()),
        (q2.x() - p2.x()),
        (p2.x() * q2.y() - q2.x() * p2.y());

    return 0.5f * (l1.transpose() * l2 + l2.transpose() * l1);
}

Conic AffineConic::intersectingLines(const Eigen::Vector2f& p1, const Eigen::Vector2f& q1,
    const Eigen::Vector2f& p2, const Eigen::Vector2f& q2) {
    Eigen::RowVector3f l1;
    l1 << (p1.y() - q1.y()),
        (q1.x() - p1.x()),
        (p1.x() * q1.y() - q1.x() * p1.y());

    Eigen::RowVector3f l2;
    l2 << (p2.y() - q2.y()),
        (q2.x() - p2.x()),
        (p2.x() * q2.y() - q2.x() * p2.y());

    return 0.5f * (l1.transpose() * l2 + l2.transpose() * l1);
}

Conic AffineConic::ellipse(const Eigen::Vector2f& center, float a, float b, float angle) {
    Eigen::Matrix3f T = translate(center) * rotate(angle) * scale(Eigen::Vector2f(a, b));
    Eigen::Matrix3f M;
    M << 1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, -1.0f;

    return applyTransform(T, M);
}

Conic AffineConic::hyperbola(const Eigen::Vector2f& center, float a, float b, float angle) {
    Eigen::Matrix3f T = translate(center) * rotate(angle) * scale(Eigen::Vector2f(a, b));
    Eigen::Matrix3f M;
    M << 1.0f, 0.0f, 0.0f,
        0.0f, -1.0f, 0.0f,
        0.0f, 0.0f, -1.0f;

    return applyTransform(T, M);
}

Conic AffineConic::parabola(const Eigen::Vector2f& vertex, float p, float angle) {
    Eigen::Matrix3f T = translate(vertex) * rotate(angle);
    Eigen::Matrix3f M;
    M << 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, -2.0f * p,
        0.0f, -2.0f * p, 0.0f;

    return applyTransform(T, M);
}

Conic AffineConic::fromPoints(const std::vector<Eigen::Vector2f>& points) {
    if (points.size() != 5) {
        throw std::invalid_argument("Exactly 5 points required to define a conic");
    }

    Eigen::MatrixXf A(5, 6);
    for (int i = 0; i < 5; i++) {
        float x = points[i].x();
        float y = points[i].y();
        A.row(i) << x * x, x* y, y* y, x, y, 1.0f;
    }

    Eigen::JacobiSVD<Eigen::MatrixXf> svd(A, Eigen::ComputeFullV);
    Eigen::VectorXf c = svd.matrixV().col(5);

    return fromCoefficients(c(0), c(1), c(2), c(3), c(4), c(5));
}

Conic AffineConic::fromCoefficients(float A, float B, float C, float M, float E, float F) {
    Eigen::Matrix3f Q;
    Q << A, B / 2.0f, M / 2.0f,
        B / 2.0f, C, E / 2.0f,
        M / 2.0f, E / 2.0f, F;
    return Q;
}