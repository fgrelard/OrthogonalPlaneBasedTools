#ifndef POLYGON_H
#define POLYGON_H

#include <initializer_list>
#include <vector>
#include <cstddef>
#include <geometry/path/BresenhamAlgorithm.h>

template<typename Point>
class Polygon {

public:
    Polygon() {}

    template <typename Iterator>
    Polygon(Iterator it, Iterator ite) {
        myPoints = std::vector<Point>(it, ite);
    }

    Polygon(std::initializer_list<Point> l) {
        myPoints = l;
    }

    std::vector<Point> getPolygon() const { return myPoints; }

    bool isInside(const Point &p) const;

private:
    std::vector<Point> myPoints;
};

template <typename Point>
bool Polygon<Point>::
isInside(const Point &p) const {
    bool inside = false;
    for (size_t c= 0, d = myPoints.size()-1; c < myPoints.size(); d = c++)
    {
        Point current = myPoints[c];
        Point other = myPoints[d];
        Point vector = other - current;
        int indexMaj = -1, indexMin = -1;
        if (vector[0] != 0) {
            indexMaj = 1;
            indexMin = 0;
        } else if (vector[1] != 0) {
            indexMaj = 0;
            indexMin = 1;
        }
        if (indexMaj != -1) {
            double slope = vector[indexMaj] / vector[indexMin];
            double b = current[indexMaj] - slope * current[indexMin];
            double value = slope * p[indexMin] + b;
            if (value == p[indexMaj]) return true;
        }


        if (((current[1] >= p[1]) != (other[1] >= p[1])) &&
            (p[0] <= ((vector[0]) * (p[1] - current[1]) /
                      (vector[1]) + current[0])))
            inside = !inside;
    }
    return inside;
}


#endif
