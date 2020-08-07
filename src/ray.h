#ifndef RAY_H
#define RAY_H
#include "cpplap.h"

class Ray;
class DiffractionPoint;
class Mirror;

class DiffractionPoint : public cpplap::Vect<float> {
    friend class Mirror;

    friend bool          operator==(const DiffractionPoint &, const DiffractionPoint &);
    friend std::ostream &operator<<(std::ostream &stream, const DiffractionPoint &a);
    friend bool          cmp_rel(const DiffractionPoint &a, const DiffractionPoint &b);

    float _amplitude;
    float _phase_shift;

public:
    DiffractionPoint(const DiffractionPoint &) = default;

    DiffractionPoint(){};

    DiffractionPoint(const float x0, const float x1, const float x2, const float amplitude = 1,
                     const float phase_shift = 0)
        : cpplap::Vect<float>(x0, x1, x2), _amplitude(amplitude), _phase_shift(phase_shift)
    {
    }

    DiffractionPoint(const cpplap::Vect<float> vect, float amplitude = 1, float phase_shift = 0)
        : cpplap::Vect<float>(vect), _amplitude(amplitude), _phase_shift(phase_shift)
    {
    }

    float &getPhaseShift() { return _phase_shift; }

    float &getAmplitude() { return _amplitude; }
};

bool cmp_rel(const DiffractionPoint &a, const DiffractionPoint &b)
{
    return cmp_rel(static_cast<const cpplap::Vect<float> &>(a), static_cast<const cpplap::Vect<float> &>(b)) &&
           cpplap::cmp_rel(a._amplitude, b._amplitude) && cpplap::cmp_rel(a._phase_shift, b._phase_shift);
}

std::ostream &operator<<(std::ostream &stream, const DiffractionPoint &a)
{
    stream << static_cast<const cpplap::Vect<float> &>(a) << "amplitude = " << a._amplitude << std::endl
           << "phase_shift = " << a._phase_shift << std::endl;
    return stream;
}

bool operator==(const DiffractionPoint &a, const DiffractionPoint &b)
{
    return (static_cast<const cpplap::Vect<float> &>(a) == static_cast<const cpplap::Vect<float> &>(b)) &&
           (a._amplitude == b._amplitude) && (a._phase_shift == b._phase_shift);
}

class Ray : public cpplap::Line<float> {
    friend class Mirror;

    float _amplitude;

public:
    Ray(const cpplap::Vect<float> r, const cpplap::Vect<float> v, const float amplitude = 1)
        : cpplap::Line<float>(r, v), _amplitude(amplitude)
    {
    }

    Ray(const cpplap::Line<float> line, const float amplitude = 1) : cpplap::Line<float>(line), _amplitude(amplitude) {}
};

class Mirror : public cpplap::ParametricPlane<float> {
    float _amplitude;

public:
    Mirror(const cpplap::Vect<float> r, const cpplap::Vect<float> u, const cpplap::Vect<float> v,
           const float amplitude = 1)
        : cpplap::ParametricPlane<float>(r, u, v), _amplitude(amplitude)
    {
    }

    bool rayInterSec(const Ray ray, DiffractionPoint &diff_point)
    {
        diff_point              = intersectionPointWith(ray);
        diff_point._amplitude   = ray._amplitude;
        diff_point._phase_shift = ray.distFromRTo(diff_point);
        return this->checkBounds(diff_point);
    }
};

#endif
