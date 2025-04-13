#include <cmath>
#include <iostream>
#include <vector>
struct Point {
  double x;
  double y;
  Point(double x, double y) : x(x), y(y) {}
  Point static middle(const Point &first, const Point &second) {
    return {(first.x + second.x) / 2, (first.y + second.y) / 2};
  }
  Point change_plus(const Point &p) const { return {x + p.x, y + p.y}; }
  Point change_multi(double factor) const { return {x * factor, y * factor}; }
  Point rotate(double angle) const {
    return {std::cos(angle) * x - std::sin(angle) * y,
            std::sin(angle) * x + std::cos(angle) * y};
  }
  Point rotate(const Point &point, double angle) const {
    return change_plus(point.change_multi(-1)).rotate(angle).change_plus(point);
  }
  Point reflect(const Point &point) {
    double delta_x = point.x - x;
    double delta_y = point.y - y;
    return change_plus({2 * delta_x, 2 * delta_y});
  }
  Point scale(const Point &point, double factor) {
    return Point{x - point.x, y - point.y}.change_multi(factor).change_plus(
        point);
  }
};
bool operator==(const Point &first, const Point &second) {
  return (std::abs(first.x - second.x) < 1e-9) &&
         (std::abs(first.y - second.y) < 1e-9);
}
bool operator!=(const Point &first, const Point &second) {
  return !(first == second);
}
namespace Utility {
struct GeometricVector {
  double x;
  double y;
  GeometricVector(Point first, Point second)
      : x(second.x - first.x), y(second.y - first.y) {}
  double static dot(GeometricVector first, GeometricVector second) {
    return first.x * second.x + first.y * second.y;
  }
  double static cross(GeometricVector first, GeometricVector second) {
    return first.x * second.y - first.y * second.x;
  }
  double length() { return std::sqrt(x * x + y * y); }
  double static angle(GeometricVector first, GeometricVector second) {
    return std::acos(dot(first, second) / (first.length() * second.length()));
  }
};
} // namespace Utility
class Line {
public:
  double a;
  double b;
  double c;
  // конструирование от канонического уравнения
  Line(double a, double b, double c) : a(a), b(b), c(c) {}
  Line(Point first, Point second) {
    a = second.y - first.y;
    b = first.x - second.x;
    c = second.x * first.y - first.x * second.y;
    toCanon();
  }
  // конструирование вида y = kx + b
  Line(double angle, double shift) {
    a = -angle;
    b = 1;
    c = -shift;
  }
  // конструирование от точки и угла
  Line(Point point, double angle) {
    a = -angle;
    b = 1;
    c = angle * point.x - point.y;
  }
  // деление канонического уравнения
  void division(double factor) {
    a /= factor;
    b /= factor;
    c /= factor;
  }
  // пересечение двух прямых
  Point intersection(const Line &line) {
    return {(b * line.c - line.b * c) / (a * line.b - line.a * b),
            (line.a * c - a * line.c) / (a * line.b - line.a * b)};
  }
  // проекция точки на прямую
  Point projection(const Point &point) const {
    return {(point.x - a * point.y - a * c) / (1 + a * a),
            (a * a * point.y - a * point.x - c) / (1 + a * a)};
  }
  // нормаль прямой к точке
  Line get_normal(const Point &point) {
    Point project = projection(point);
    Line normal(-b, a, (b - a) * project.x - (a + b) * project.y - c);
    normal.toCanon();
    return normal;
  }
  Point reflect(const Point &point) const {
    Point project = projection(point);
    Point to_add = {2 * (project.x - point.x), 2 * (project.y - point.y)};
    return point.change_plus(to_add);
  }

private:
  // приводим к нормальному виду
  void toCanon() {
    if (b == 0) {
      division(a);
    } else {
      division(b);
    }
  }
};
// операторы сравнения
bool operator==(const Line &first, const Line &second) {
  std::cerr << ((std::abs(first.a - second.b) < 1e-9) &&
                (std::abs(first.b - second.b) < 1e-9) &&
                (std::abs(first.c - second.c) < 1e-9))
            << '\n';
  return !((std::abs(first.a - second.b) < 1e-9) &&
           (std::abs(first.b - second.b) < 1e-9) &&
           (std::abs(first.c - second.c) < 1e-9));
}
bool operator!=(const Line &first, const Line &second) {
  return !(first == second);
}
// абстрактный класс
class Shape {
public:
  virtual ~Shape() = default;
  virtual double perimeter() const = 0;
  virtual double area() const = 0;
  virtual bool isCongruentTo(const Shape &another) = 0;
  virtual bool isSimilarTo(const Shape &another) = 0;
  virtual bool containsPoint(const Point &point) = 0;
  virtual void rotate(const Point &center, double angle) = 0;
  virtual void reflect(const Point &center) = 0;
  virtual void reflect(const Line &axis) = 0;
  virtual void scale(const Point &center, double coefficient) = 0;
  virtual bool equal(const Shape &another) const = 0;
};
bool operator==(const Shape &shape_left, const Shape &shape_right) {
  return shape_left.equal(shape_right);
}
bool operator!=(const Shape &shape_left, const Shape &shape_right) {
  return !(shape_left == shape_right);
}
class Polygon : public Shape {
protected:
  std::vector<Point> vertices;

public:
  Polygon() = default;
  Polygon(std::vector<Point> vertex) { vertices = vertex; }
  template <typename... Args> Polygon(Args... args) {
    (vertices.push_back(args), ...);
  }
  size_t verticesCount() const { return vertices.size(); }
  std::vector<Point> getVertices() const { return vertices; }
  bool isConvex() {
    using Utility::GeometricVector;
    size_t N = verticesCount();
    std::vector<Point> points = getVertices();
    double cp = 0;
    // double prev = 1;
    if (N <= 3) {
      return true;
    }
    bool sign = false;
    for (size_t i = 0; i < N; i++) {
      cp = GeometricVector::cross(
          GeometricVector(vertices[i], vertices[(i + 1) % N]),
          GeometricVector(vertices[(i + 1) % N], vertices[(i + 2) % N]));
      if (cp != 0) {
        if (!sign) {
          sign = (cp > 0);
        } else if ((cp > 0) != sign) {
          return false;
        }
      }
    }
    return true;
  }
  double perimeter() const override {
    size_t N = vertices.size();
    double ans = 0;
    for (size_t i = 0; i < N; i++) {
      ans +=
          Utility::GeometricVector(vertices[i], vertices[(i + 1) % N]).length();
    }
    return ans;
  }
  double area() const override {
    size_t N = vertices.size();
    double ans = 0;
    for (size_t i = 1; i < N - 1; i++) {
      ans += Utility::GeometricVector::cross(
                 Utility::GeometricVector(vertices[0], vertices[i]),
                 Utility::GeometricVector(vertices[0], vertices[i])) /
             2;
    }
    return std::abs(ans);
  }
  bool equal(const Shape &another) const override {
    const Polygon *main = dynamic_cast<const Polygon *>(&another);
    if (main) {
      size_t count_first = vertices.size();
      size_t count_second = main->verticesCount();
      const std::vector<Point> &main_vertex = main->getVertices();
      size_t i = 0;
      size_t j = 0;
      while (j < count_second && vertices[i] != main_vertex[j]) {
        ++j;
      }
      if (j == count_second) {
        return false;
      }
      if (count_first != count_second) {
        return false;
      }
      if (vertices[i + 1] == main_vertex[(j + 1) % count_second]) {
        ++i;
        size_t start = j;
        for (size_t k = (start + 1) % count_second; k != start;
             k = (k + 1) % count_second) {
          if (vertices[i++] != main_vertex[k]) {
            return false;
          }
        }
        return true;
      }
      ++i;
      size_t start = j;
      for (size_t k = (start - 1 + count_second) % count_second; k != start;
           k = (k - 1 + count_second) % count_second) {
        if (vertices[i++] != main_vertex[k]) {
          return false;
        }
      }
      return true;
    }
    return false;
  }
  static bool IsCongruent(const std::vector<Point> &vertices_first,
                          const std::vector<Point> &vertices_second) {
    using Utility::GeometricVector;
    size_t n = vertices_first.size();
    for (size_t i = 0; i < n; ++i) {
      bool is_congruent = true;
      for (size_t j = 0; j < n; ++j) {
        GeometricVector edge1(vertices_first[(i + j) % n],
                              vertices_first[(i + j + 1) % n]);
        GeometricVector edge2(vertices_second[j], vertices_second[(j + 1) % n]);
        if (std::abs(edge1.length() - edge2.length()) > 1e-9) {
          is_congruent = false;
          break;
        }
        double area_first = GeometricVector::cross(
            edge1, GeometricVector(vertices_first[(i + j) % n],
                                   vertices_first[(i + j - 1 + n) % n]));
        double area_second = GeometricVector::cross(
            edge2, GeometricVector(vertices_second[j],
                                   vertices_second[(j - 1 + n) % n]));
        if (std::abs(std::abs(area_first) - std::abs(area_second)) > 1e-9) {
          is_congruent = false;
          break;
        }
      }
      if (is_congruent) {
        return true;
      }
    }
    return false;
  }
  static bool IsSimilar(const std::vector<Point> &vertices_first,
                        const std::vector<Point> &vertices_second) {
    using Utility::GeometricVector;

    size_t n = vertices_first.size();
    for (size_t i = 0; i < n; ++i) {
      bool is_similar = true;
      double k = -1;
      for (size_t j = 0; j < n; ++j) {
        GeometricVector edge1(vertices_first[(i + j) % n],
                              vertices_first[(i + j + 1) % n]);
        GeometricVector edge2(vertices_second[j], vertices_second[(j + 1) % n]);
        double factor_similar = edge1.length() / edge2.length();
        if (k == -1) {
          k = factor_similar;
        } else if (std::abs(k - factor_similar) > 1e-9) {
          is_similar = false;
          break;
        }
        double area_first = GeometricVector::cross(
            edge1, GeometricVector(vertices_first[(i + j) % n],
                                   vertices_first[(i + j - 1 + n) % n]));
        double area_second = GeometricVector::cross(
            edge2, GeometricVector(vertices_second[j],
                                   vertices_second[(j - 1 + n) % n]));
        if (std::abs(std::abs(area_first) - k * k * std::abs(area_second)) >
            1e-9) {
          is_similar = false;
          break;
        }
      }
      if (is_similar) {
        return true;
      }
    }
    return false;
  }
  bool isCongruentTo(const Shape &another) override {
    const Polygon *polygon = dynamic_cast<const Polygon *>(&another);
    if (polygon) {
      if (vertices.size() != polygon->verticesCount()) {
        return false;
      }
      std::vector<Point> vertices_other = polygon->getVertices();
      if (IsCongruent(vertices, vertices_other)) {
        return true;
      }
      std::reverse(vertices_other.begin(), vertices_other.end());
      return IsCongruent(vertices, vertices_other);
    }
    return false;
  }

  bool isSimilarTo(const Shape &another) override {
    const Polygon *polygon = dynamic_cast<const Polygon *>(&another);
    if (polygon) {
      if (vertices.size() != polygon->verticesCount()) {
        return false;
      }
      std::vector<Point> vertices_other = polygon->getVertices();
      if (IsSimilar(vertices, vertices_other)) {
        return true;
      }
      std::reverse(vertices_other.begin(), vertices_other.end());
      return (IsSimilar(vertices, vertices_other));
    }
    return false;
  }
  bool containsPoint(const Point &point) override {
    using Utility::GeometricVector;
    size_t sz = vertices.size();
    double sign = 0;
    for (size_t i = 0; i < sz; ++i) {
      double cross = GeometricVector::cross(
          GeometricVector(vertices[i], point),
          GeometricVector(vertices[i], vertices[(i + 1) % sz]));
      if (cross == 0) {
        sign = 0;
      } else if (cross != 0 && cross * sign < 0) {
        return false;
      }
      sign = cross;
    }
    return true;
  }
  void rotate(const Point &center, double angle) override {
    for (Point &p : vertices) {
      p = p.rotate(center, angle);
    }
  }

  void reflect(const Point &center) override {
    for (Point &p : vertices) {
      p = p.reflect(center);
    }
  }

  void reflect(const Line &axis) override {
    for (Point &p : vertices) {
      p = axis.reflect(p);
    }
  }

  void scale(const Point &center, double factor) override {
    for (Point &p : vertices) {
      p = p.scale(center, factor);
    }
  }
};
class Ellipse : public Shape {
protected:
  std::pair<Point, Point> points_focuces = {{0, 0}, {0, 0}};
  double big_axis;

public:
  Ellipse(Point p1, Point p2, double length) {
    if (p1.x < p2.x || (std::abs(p1.x - p2.y) < 1e-9 && p1.y < p2.y)) {
      points_focuces = {p1, p2};
    } else {
      points_focuces = {p2, p1};
    }
    big_axis = length / 2;
  }
  std::pair<Point, Point> focuses() { return points_focuces; }
  std::pair<Line, Line> directrices() {
    Line axis(points_focuces.first, points_focuces.second);
    return {axis.get_normal(points_focuces.first),
            axis.get_normal(points_focuces.second)};
  }
  Point center() const {
    return {(points_focuces.first.x + points_focuces.second.x) / 2,
            (points_focuces.first.y + points_focuces.second.y) / 2};
  }
  // расстоние от центра до фокуса
  double focus_length() const {
    return Utility::GeometricVector(points_focuces.first, center()).length();
  }
  // большая полуось^2 = малая полуось^2 + фокальное расстоние^2
  double small_axis() const {
    return std::sqrt(big_axis * big_axis - focus_length() * focus_length());
  }
  Ellipse toCanon() const {
    return Ellipse({0, 0},
                   points_focuces.second.change_plus(
                       points_focuces.first.change_multi(-1)),
                   2 * big_axis);
  }
  // эксцентриситет
  double eccentricity() { return focus_length() / big_axis; }
  // периметр
  double perimeter() const override {
    double a = big_axis;
    double b = small_axis();
    double h = (std::pow(a - b, 2)) / (std::pow(a + b, 2));
    double perimeters =
        M_PI * (a + b) * (1 + (3 * h) / (10 + std::sqrt(4 - 3 * h)));
    return perimeters;
  }
  // площадь
  double area() const override { return M_PI * big_axis * small_axis(); }
  bool equal(const Shape &another) const override {
    const Ellipse *ellipse = dynamic_cast<const Ellipse *>(&another);
    if (ellipse) {
      return (points_focuces.first == ellipse->points_focuces.first &&
              points_focuces.second == ellipse->points_focuces.second &&
              (big_axis - ellipse->big_axis) < 1e-9);
    }
    return false;
  }
  // лежит ли точка
  bool containsPoint(const Point &point) override {
    Point centre = center();
    double x = point.x - centre.x;
    double y = point.y - centre.y;
    return (x * x / (big_axis * big_axis) +
            y * y / (small_axis() * small_axis()) - 1) < 1e-9;
  }
  // геом. равенсво
  bool isCongruentTo(const Shape &another) override {
    const Ellipse *ellipse = dynamic_cast<const Ellipse *>(&another);
    using Utility::GeometricVector;
    if (ellipse) {
      Ellipse canonic_first = toCanon();
      Ellipse canonic_second = ellipse->toCanon();
      double angle1 = GeometricVector::angle(
          GeometricVector({0, 0}, {1, 0}),
          GeometricVector(canonic_first.points_focuces.first,
                          canonic_first.points_focuces.second));
      double angle2 = GeometricVector::angle(
          GeometricVector({0, 0}, {1, 0}),
          GeometricVector(canonic_second.points_focuces.first,
                          canonic_second.points_focuces.second));
      canonic_first.rotate({0, 0}, -angle1);
      canonic_second.rotate({0, 0}, -angle2);
      return canonic_first == canonic_second;
    }
    return false;
  }
  // подобие равенство
  bool isSimilarTo(const Shape &another) override {
    const Ellipse *ellipse = dynamic_cast<const Ellipse *>(&another);
    using Utility::GeometricVector;
    if (ellipse) {
      Ellipse canonic_first = toCanon();
      Ellipse canonic_second = ellipse->toCanon();
      double k1 = 1 / GeometricVector(canonic_first.points_focuces.first,
                                      canonic_first.points_focuces.second)
                          .length();
      double k2 = 1 / GeometricVector(canonic_second.points_focuces.first,
                                      canonic_second.points_focuces.second)
                          .length();
      canonic_first.scale({0, 0}, k1);
      canonic_second.scale({0, 0}, k2);
      return canonic_first.isCongruentTo(canonic_second);
    }
    return false;
  }
  void rotate(const Point &center, double angle) override {
    points_focuces.first = points_focuces.first.rotate(center, angle);
    points_focuces.second = points_focuces.second.rotate(center, angle);
  }

  void reflect(const Point &center) override {
    points_focuces.first = points_focuces.first.reflect(center);
    points_focuces.second = points_focuces.second.reflect(center);
  }

  void reflect(const Line &axis) override {
    points_focuces.first = axis.reflect(points_focuces.first);
    points_focuces.second = axis.reflect(points_focuces.second);
  }

  void scale(const Point &center, double factor) override {
    points_focuces.first = points_focuces.first.scale(center, factor);
    points_focuces.second = points_focuces.second.scale(center, factor);
    big_axis *= factor;
  }
};
class Circle : public Ellipse {
public:
  Circle(Point center, double radius) : Ellipse(center, center, 2 * radius) {}
  double radius() { return big_axis; }
};
class Rectangle : public Polygon {
public:
  Rectangle() = default;
  Rectangle(const Point &point_one, const Point &point_three, double ratio)
      : Polygon() {
    Point center = Point::middle(point_one, point_three);
    if (ratio < 1) {
      ratio = 1 / ratio;
    }
    Point point_second = point_one.rotate(center, M_PI - 2 * std::atan(ratio));
    Point point_fourth = point_second.reflect(center);
    vertices = {point_one, point_second, point_three, point_fourth};
  }
  double perimeter() const override {
    using Utility::GeometricVector;
    return 2 * (GeometricVector(vertices[0], vertices[1]).length() +
                GeometricVector(vertices[1], vertices[2]).length());
  }

  double area() const override {
    using Utility::GeometricVector;
    return GeometricVector(vertices[0], vertices[1]).length() *
           GeometricVector(vertices[1], vertices[2]).length();
  }

  Point center() const { return Point::middle(vertices[0], vertices[2]); }
  std::pair<Line, Line> diagonals() {
    return {Line(vertices[0], vertices[2]), Line(vertices[1], vertices[3])};
  }
};
class Square : public Rectangle {
public:
  Square(const Point &point1, const Point &point3) : Rectangle() {
    Point middle = Point::middle(point1, point3);
    Point shift = point1.change_plus(middle.change_multi(-1));
    Point point2 = {-shift.y + middle.x, shift.x + middle.y};
    Point point4 = {shift.y + middle.x, -shift.x + middle.y};
    vertices = {point1, point2, point3, point4};
  }
  Circle circumscribedCircle() const {
    return Circle(center(),
                  Utility::GeometricVector(vertices[0], center()).length());
  }

  Circle inscribedCircle() const {
    return Circle(center(),
                  Utility::GeometricVector(vertices[0], vertices[1]).length() /
                      2);
  }
};
class Triangle : public Polygon {
public:
  Triangle(const Point &point1, const Point &point2, const Point &point3) {
    vertices = {point1, point2, point3};
  }
  Circle circumscribedCircle() const {
    Point middle_first = Point::middle(vertices[0], vertices[1]);
    Point middle_second = Point::middle(vertices[1], vertices[2]);
    Line bisector_first =
        Line(vertices[0], vertices[1]).get_normal(middle_first);
    Line bisector_second =
        Line(vertices[1], vertices[2]).get_normal(middle_second);
    Point inter = bisector_first.intersection(bisector_second);
    return Circle(inter, Utility::GeometricVector(inter, vertices[0]).length());
  }

  Circle inscribedCircle() const {
    using Utility::GeometricVector;
    GeometricVector side_first(vertices[0], vertices[1]);
    GeometricVector side_second(vertices[1], vertices[2]);
    GeometricVector side_third(vertices[2], vertices[0]);
    double k1 = side_first.length() / side_second.length();
    double k2 = side_second.length() / side_third.length();
    Point point_first = {(vertices[0].x + k1 * vertices[2].x) / (1 + k1),
                         (vertices[0].y + k1 * vertices[2].y) / (1 + k1)};
    Point point_second = {(vertices[1].x + k2 * vertices[0].x) / (1 + k2),
                          (vertices[1].y + k2 * vertices[0].y) / (1 + k2)};
    Line median_first = Line(vertices[1], point_first);
    Line median_second = Line(vertices[2], point_second);
    Point internal = median_first.intersection(median_second);
    Point projection = Line(vertices[0], vertices[1]).projection(internal);
    return Circle(internal, GeometricVector(internal, projection).length());
  }
  Point centroid() const {
    double x = (vertices[0].x + vertices[1].x + vertices[2].x) / 3;
    double y = (vertices[0].y + vertices[1].y + vertices[2].y) / 3;
    return {x, y};
  }

  Point orthocenter() const {
    Line side_first = Line(vertices[0], vertices[1]);
    Line side_second = Line(vertices[1], vertices[2]);
    return side_first.get_normal(vertices[2])
        .intersection(side_second.get_normal(vertices[0]));
  }

  Line EulerLine() const { return Line(centroid(), orthocenter()); }

  Circle ninePointsCircle() const {
    Point m1 = Point::middle(vertices[0], vertices[1]);
    Point m2 = Point::middle(vertices[1], vertices[2]);
    Point m3 = Point::middle(vertices[0], vertices[2]);
    return Triangle(m1, m2, m3).circumscribedCircle();
  }
};
