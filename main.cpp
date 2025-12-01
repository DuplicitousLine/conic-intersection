#include <iostream>
#include <vector>
#include <cmath>

#include <Eigen/Dense>

#include "AffineConic.h"
#include "ConicIntersection.h"

static void printConicConstruction(const std::string& name, const Eigen::Vector2f& t, float a, const Eigen::Vector2f& s) {
    std::cout << name << "\n"
        << "  Translation (" << t.x() << ", " << t.y() << ")\n"
        << "  Rotation    " << a << " radians\n"
        << "  Scale       (" << s.x() << ", " << s.y() << ")\n\n";
}

static void printIntersection(const std::string& label, const std::vector<Eigen::Vector3f>& points) {
    std::cout << label << "\n";
    if (points.empty()) {
        std::cout << "  No intersection points found.\n";
    }
    else {
        for (const auto& p : points) {
            if (std::abs(p.z()) > 1e-6f) {
                std::cout << "  (" << p.x() / p.z() << ", " << p.y() / p.z() << ")\n";
            }
            else {
                std::cout << "  (" << p.x() << ", " << p.y() << ", " << p.z() << ")\n";
            }
        }
    }
    std::cout << "\n";
}

int main() {
    Eigen::Vector2f t1(0.0f, 0.0f);
    float a1 = 0.0f;
    Eigen::Vector2f s1(1.0f, 2.0f);

    Eigen::Vector2f t2(0.5f, 0.5f);
    float a2 = 0.5f;
    Eigen::Vector2f s2(1.5f, 1.0f);

    Eigen::Vector2f t3(1.0f, 0.0f);
    float a3 = 0.0f;
    Eigen::Vector2f s3(1.0f, 1.0f);

    Eigen::Vector2f t4(0.0f, 1.0f);
    float a4 = 0.0f;
    Eigen::Vector2f s4(1.0f, 1.0f);

    Conic ellipse1 = AffineConic::getEllipseMatrix(t1, a1, s1);
    Conic ellipse2 = AffineConic::getEllipseMatrix(t2, a2, s2);
    Conic parabola = AffineConic::getParabolaMatrix(t3, a3, s3);
    Conic hyperbola = AffineConic::getHyperbolaMatrix(t4, a4, s4);

    printConicConstruction("Ellipse 1", t1, a1, s1);
    printConicConstruction("Ellipse 2", t2, a2, s2);
    printConicConstruction("Parabola", t3, a3, s3);
    printConicConstruction("Hyperbola", t4, a4, s4);

    ConicIntersection engine;

    printIntersection("Intersection between Ellipse 1 and Ellipse 2:", engine.intersect(ellipse1, ellipse2));
    printIntersection("Intersection between Ellipse 1 and Parabola:", engine.intersect(ellipse1, parabola));
    printIntersection("Intersection between Ellipse 1 and Hyperbola:", engine.intersect(ellipse1, hyperbola));

    return 0;
}
