#ifndef DMD_H
#define DMD_H
#include "cpplap/cpplap.h"
#include "ray.h"
#include <limits>
#include <tuple>
#include <vector>

struct DmdParams {
    int    nx;          // Number mirrors in x
    int    ny;          // Number mirrors in y
    double angle;       // Flip angle
    double mirror_size; // Length of mirror edge
    double mirror_dist; // Distance between mirrors
};

std::ostream &operator<<(std::ostream &stream, const DmdParams &a)
{
    stream << "nx " << a.nx << std::endl
           << "ny " << a.ny << std::endl
           << "angle " << a.angle << std::endl
           << "mirror_size " << a.mirror_size << std::endl
           << "mirror_dist " << a.mirror_dist << std::endl;
    return stream;
}

class Dmd {

    friend std::ostream &operator<<(std::ostream &stream, const Dmd &a);

    std::vector<char> _flipped;
    // Dmd parameters
    DmdParams _p;
    int       _nmirrors;

    double _eff_width; // mirror_dist + mirror_size
    double _size_x;    // Total dmd size in x direction
    double _size_y;    // Total dmd size in y direction

    // cpplap::Vect<float>or along edge of mirror in x direction
    cpplap::Vect<float> _dir_x_on;
    cpplap::Vect<float> _dir_x_off;

    // cpplap::Vect<float>or along edge of mirror in y direction
    cpplap::Vect<float> _dir_y_on;
    cpplap::Vect<float> _dir_y_off;

    // Normal cpplap::Vect<float>or of mirror
    cpplap::Vect<float> _dir_n_on;
    cpplap::Vect<float> _dir_n_off;

    cpplap::Vect<float> _mirror_axis = cpplap::Vect<float>(1, 1, 0);

    cpplap::Vect<float> planeLineIntersec(const cpplap::Line<float> &line) const
    {
        return cpplap::HessePlane<float>(cpplap::Vect<float>(0, 0, 0), cpplap::Vect<float>(0, 0, 1))
            .intersectionPointWith(line);
    }

public:
    Dmd(DmdParams params)
        :

          _p(params), _nmirrors(_p.nx * _p.ny), _eff_width(_p.mirror_size + _p.mirror_dist),
          _size_x(_p.nx * _eff_width), _size_y(_p.ny * _eff_width)

    {

        _flipped.resize(_nmirrors, 0);

        // stripe pattern
        // for (int i = 0; i < _p.nx; i++){
        //    for (int j = 0; j < _p.ny; j++){
        //        _flipped[index(i,j)] = (j/3 % 2 == 0);
        //    }
        //}

        _dir_x_on  = (_p.mirror_size / 2 * cpplap::Vect<float>(1, 0, 0)).rotate(_mirror_axis, _p.angle);
        _dir_x_off = (_p.mirror_size / 2 * cpplap::Vect<float>(1, 0, 0)).rotate(_mirror_axis, -_p.angle);

        _dir_y_on  = (_p.mirror_size / 2 * cpplap::Vect<float>(0, 1, 0)).rotate(_mirror_axis, _p.angle);
        _dir_y_off = (_p.mirror_size / 2 * cpplap::Vect<float>(0, 1, 0)).rotate(_mirror_axis, -_p.angle);

        _dir_n_on  = (cpplap::Vect<float>(0, 0, 1)).rotate(_mirror_axis, _p.angle);
        _dir_n_off = (cpplap::Vect<float>(0, 0, 1)).rotate(_mirror_axis, -_p.angle);
    }

    inline int index(const int x, const int y) const
    {
        if (x >= _p.nx)
            throw std::out_of_range("Dmd: Index x out of range");
        if (y >= _p.ny)
            throw std::out_of_range("Dmd: Index out of range");
        return _p.nx * y + x;
    }

    inline int index(const cpplap::Vect<float> &vect) const
    {
        auto [x, y] = findIndXY(vect);
        return index(int(x), int(y));
    }

    inline std::tuple<int, int> findIndXY(const cpplap::Vect<float> &vect) const
    {

        if (vect[0] < -_size_x / 2 || vect[0] >= _size_y / 2) {
            throw std::domain_error("Dmd: cpplap::Vect<float>or index x out of range");
        }

        if (vect[1] < -_size_y / 2 || vect[1] >= _size_y / 2) {
            throw std::domain_error("Dmd: cpplap::Vect<float>or index y out of range");
        }

        double xd = vect[0] + _size_x / 2.;
        double yd = vect[1] + _size_y / 2.;
        int    x  = xd / _eff_width;
        int    y  = yd / _eff_width;
        return std::tuple(x, y);
    }

    std::tuple<int, int> de_index(const int i) const
    {
        int y = i / _p.nx;
        int x = i % _p.nx;
        return std::tuple(x, y);
    }

    inline int up_x(const int i, const int n) const
    {
        const auto [x, y] = de_index(i);
        return index(x + n, y);
    }

    inline int up_y(const int i, const int n) const
    {
        const auto [x, y] = de_index(i);
        return index(x, y + n);
    }

    Mirror getMirror(const int i) const
    {
        const auto [x, y] = de_index(i);
        return getMirror(x, y);
    }

    std::tuple<double, double> getMirrorCoordinates(const int i) const
    {
        const auto [ix, iy] = de_index(i);
        return getMirrorCoordinates(ix, iy);
    }

    std::tuple<char, double, double> getMirrorStatus(const int i) const
    {
        const auto [ix, iy] = de_index(i);
        const auto [x, y]   = getMirrorCoordinates(ix, iy);
        return std::tuple(_flipped[i], x, y);
    }

    std::tuple<double, double> getMirrorCoordinates(const int i, const int j) const
    {
        const double x = i * _eff_width - _size_x / 2 + _eff_width / 2;
        const double y = j * _eff_width - _size_y / 2 + _eff_width / 2;
        return std::tuple(x, y);
    }

    Mirror getMirror(const int x, const int y) const
    {

        if (x >= _p.nx)
            throw std::out_of_range("Dmd::getMirror: Index x out of range");
        if (y >= _p.ny)
            throw std::out_of_range("Dmd::getMirror: Index out of range");

        const auto [xd, yd] = getMirrorCoordinates(x, y);
        cpplap::Vect<float> base_vect(xd, yd, 0);

        if (_flipped[index(x, y)]) {
            return Mirror(base_vect, _dir_x_on, _dir_y_on);
        }
        else {
            return Mirror(base_vect, _dir_x_off, _dir_y_off);
        }
    }

    /* Find mirror that is first hit by a given ray.
     * First find at which point we hit the dmd plane. Then search for
     * intersection points with mirrors surrounding this hit point.
     */
    DiffractionPoint findDiffPoint(const Ray &ray) const
    {
        // Intersection with dmd plane (i.e. x-y-plane)
        cpplap::Vect<float> dmd_intersec = planeLineIntersec(ray);

        int x, y;
        try {
            // find the corresponding mirror coordinates
            std::tie(x, y) = findIndXY(dmd_intersec);
        }
        catch (std::domain_error &e) {
            // If there is no mirror return directly with amplitude 0
            return DiffractionPoint(dmd_intersec, 0, ray.dist(dmd_intersec));
        }

        // Diffraction point with current mirror
        DiffractionPoint diff_point;

        // First mirror up the ray. Initialize with dmd_intersec.
        // If we do not hit a mirror, this will be our return value with amplitude 0
        DiffractionPoint first_diff_point(dmd_intersec, 0, ray.dist(dmd_intersec));

        // Distance of the ray base vector and the diffraction point
        double last_dist = std::numeric_limits<double>::max();

        // Test mirrors surrounding dmd_intersec
        for (int i = -1 * (x > 0); i <= (x < _p.nx - 1); i++) {
            for (int j = -1 * (y > 0); j <= (y < _p.nx - 1); j++) {

                Mirror mirror = getMirror(x + i, y + j);

                // Test if we hit that mirror
                if (mirror.rayInterSec(ray, diff_point)) {

                    // Calculate the distance within ray
                    double dist = diff_point.getPhaseShift();

                    // If distance is smaller than before, select this as the new
                    // diffraction point
                    if (dist < last_dist) {
                        first_diff_point = diff_point;
                        last_dist        = dist;
                    }
                }
            }
        }

        return first_diff_point;
    }

    double getMirrorSize() const { return _p.mirror_size; }

    double getMirrorAngle() const { return _p.angle; }

    int nMirrors() const { return _nmirrors; }
};

std::ostream &operator<<(std::ostream &stream, const Dmd &a)
{
    stream << a._p << "Mirror Axis " << std::endl
           << a._mirror_axis << "Mirror direction x OFF " << std::endl
           << a._dir_x_off << std::endl
           << "Mirror direction x ON " << std::endl
           << a._dir_x_on << std::endl
           << "Mirror direction y OFF " << std::endl
           << a._dir_y_off << std::endl
           << "Mirror direction y ON " << std::endl
           << a._dir_x_on << std::endl;

    return stream;
}

#endif
