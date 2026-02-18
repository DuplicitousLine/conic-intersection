#include <iostream>
#include <Eigen/Core>
#include <vector>
#include <cassert>
#include "ConicIntersection.h"
#include "AffineConic.h"

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

static void printEllipseParams(const Eigen::Vector2f& center, float radiusX, float radiusY, float rotation) {
    std::cout << "Ellipse parameters:\n"
        << "  Center: (" << center(0) << ", " << center(1) << ")\n"
        << "  Radius X: " << radiusX << ", Radius Y: " << radiusY << "\n"
        << "  Rotation: " << rotation << "\n";
}

static void printLineParams(const Eigen::Vector2f& point1, const Eigen::Vector2f& point2) {
    std::cout << "Line parameters:\n"
        << "  Point 1: (" << point1(0) << ", " << point1(1) << ")\n"
        << "  Point 2: (" << point2(0) << ", " << point2(1) << ")\n";
}

static void printCircleParams(const Eigen::Vector2f& center, float radius) {
    std::cout << "Circle parameters:\n"
        << "  Center: (" << center(0) << ", " << center(1) << ")\n"
        << "  Radius: " << radius << "\n";
}

static void printHyperbolaParams(const Eigen::Vector2f& center, float a, float b) {
    std::cout << "Hyperbola parameters:\n"
        << "  Center: (" << center(0) << ", " << center(1) << ")\n"
        << "  a: " << a << ", b: " << b << "\n";
}

static void printParabolaParams(const Eigen::Vector2f& vertex, float p) {
    std::cout << "Parabola parameters:\n"
        << "  Vertex: (" << vertex(0) << ", " << vertex(1) << ")\n"
        << "  p: " << p << "\n";
}

static void printParallelLines(const Eigen::Vector2f& point, const Eigen::Vector2f& direction, float distance) {
    std::cout << "Parallel Line Parameters:\n";
    std::cout << "  Point: (" << point.x() << ", " << point.y() << ")\n";
    std::cout << "  Direction: (" << direction.x() << ", " << direction.y() << ")\n";
    std::cout << "  Distance between parallel lines: " << distance << "\n";
}

// Test case 1: Two ellipses that intersect at two points
void testEllipseIntersection(ConicIntersection & conicIntersection, AffineConic & affineConic) {
    std::cout << "\n==== Testing Ellipse Intersection ====\n";
    Eigen::Vector2f center1(0.0f, 0.0f);
    Eigen::Vector2f center2(1.0f, 1.0f);
    float radiusX1 = 3.0f;
    float radiusY1 = 2.0f;
    float rotation1 = 0.0f;

    float radiusX2 = 3.0f;
    float radiusY2 = 2.0f;
    float rotation2 = 0.0f;

    printEllipseParams(center1, radiusX1, radiusY1, rotation1);
    printEllipseParams(center2, radiusX2, radiusY2, rotation2);

    Conic ellipse1 = affineConic.ellipse(center1, radiusX1, radiusY1, rotation1);
    Conic ellipse2 = affineConic.ellipse(center2, radiusX2, radiusY2, rotation2);
    auto ellipseIntersections = conicIntersection.intersect(ellipse1, ellipse2);
    printIntersection("Ellipse Intersection", ellipseIntersections);
}

// Test case 2: A circle and a hyperbola that do not intersect
void testCircleAndHyperbolaIntersection(ConicIntersection & conicIntersection, AffineConic & affineConic) {
    std::cout << "\n==== Testing Circle and Hyperbola Intersection (No Intersection) ====\n";
    Eigen::Vector2f circleCenter(0.0f, 0.0f);
    float circleRadius = 5.0f;
    printCircleParams(circleCenter, circleRadius);
    Conic circle = affineConic.ellipse(circleCenter, circleRadius, circleRadius, 0.0f);

    Eigen::Vector2f hyperbolaCenter(0.0f, 0.0f);
    float hyperbolaA = 3.0f;
    float hyperbolaB = 2.0f;
    printHyperbolaParams(hyperbolaCenter, hyperbolaA, hyperbolaB);
    Conic hyperbola = affineConic.hyperbola(hyperbolaCenter, hyperbolaA, hyperbolaB, 0.0f);

    auto noIntersection = conicIntersection.intersect(circle, hyperbola);
    printIntersection("Circle and Hyperbola Intersection", noIntersection);
}

// Test case 3: Parallel lines (should be no intersection)
void testParallelLinesIntersection(ConicIntersection& conicIntersection, AffineConic& affineConic) {
    std::cout << "\n==== Testing Parallel Lines Intersection (No Intersection) ====\n";

    // Define points and direction for the parallel lines
    Eigen::Vector2f line1Point(0.0f, 1.0f);    // A point on line 1
    Eigen::Vector2f line1Direction(1.0f, 0.0f); // Direction of the lines (horizontal in this case)
    float distanceBetweenLines = 2.0f;          // Distance between the two parallel lines

    // Print the parameters of the parallel lines
    printParallelLines(line1Point, line1Direction, distanceBetweenLines);
    printParallelLines(line1Point + Eigen::Vector2f(0.0f, 1.0f), line1Direction, distanceBetweenLines);

    // Create the two parallel lines using the affineConic.parallelLines method
    Conic parallelLines = affineConic.parallelLines(line1Point, line1Direction, distanceBetweenLines);
    Conic parallelLines2 = affineConic.parallelLines(line1Point + Eigen::Vector2f(0.0f, 1.0f), line1Direction, distanceBetweenLines);

    // Test intersection: Parallel lines do not intersect, so the result should be empty
    auto parallelIntersections = conicIntersection.intersect(parallelLines, parallelLines2);
    printIntersection("Parallel Lines Intersection", parallelIntersections);
}

// Test case 4: Two lines that intersect at a point
void testLineLineIntersection(ConicIntersection & conicIntersection, AffineConic & affineConic) {
    std::cout << "\n==== Testing Line-Line Intersection ====\n";
    Eigen::Vector2f line3Point1(0.0f, 0.0f);
    Eigen::Vector2f line3Point2(1.0f, 1.0f);
    Eigen::Vector2f line4Point1(1.0f, 0.0f);
    Eigen::Vector2f line4Point2(0.0f, 1.0f);

    printLineParams(line3Point1, line3Point2);
    printLineParams(line4Point1, line4Point2);

    Conic line1 = affineConic.line(line3Point1, line3Point2);
    Conic line2 = affineConic.line(line4Point1, line4Point2);
    auto lineIntersections = conicIntersection.intersect(line1, line2);
    printIntersection("Line-Line Intersection", lineIntersections);
}

// Test case 5: Point on ellipse
void testPointEllipseIntersection(ConicIntersection& conicIntersection, AffineConic& affineConic) {
    std::cout << "\n==== Testing Point and Ellipse Intersection ====\n";

    // Define the point to test (a point on the ellipse)
    Eigen::Vector2f point(2.0f, 0.0f);  // Point (1, 1)

    // Create a degenerate conic (representing the point (1, 1))
    Conic pointConic = affineConic.point(point);
    std::cout << "Point Conic Matrix:\n" << pointConic << std::endl;

    // Create an ellipse centered at (0, 0) with radius 2 and 3, no rotation
    Eigen::Vector2f center(0.0f, 0.0f);
    float radiusX = 2.0f;
    float radiusY = 3.0f;
    float rotation = 0.0f;  // No rotation

    Conic ellipseConic = affineConic.ellipse(center, radiusX, radiusY, rotation);
    std::cout << "Ellipse Conic Matrix:\n" << ellipseConic << std::endl;

    // Test the intersection: Check if the point lies on the ellipse
    auto intersections = conicIntersection.intersect(pointConic, ellipseConic);

    // Print the intersection results
    printIntersection("Point and Ellipse Intersection", intersections);
}

int main() {
    ConicIntersection conicIntersection;
    AffineConic affineConic;
    
    testEllipseIntersection(conicIntersection, affineConic);
    testCircleAndHyperbolaIntersection(conicIntersection, affineConic);
    testParallelLinesIntersection(conicIntersection, affineConic);
    testLineLineIntersection(conicIntersection, affineConic);
    testPointEllipseIntersection(conicIntersection, affineConic);

    return 0;
}
