#ifndef RAYTRACING_H
#define RAYTRACING_H
#include "dmd.h"
#include "measure_time.h"
#include <algorithm>
#include <complex>
#include <random>
#include <thread>

#define BLOCK_SIZE 10

template <class TDist = float> TDist gauss(TDist s, TDist m, TDist x)
{
    return (1 / (s * sqrt(2 * M_PI))) * exp(-0.5 * pow((x - m) / s, 2.0));
}

class RayTracer {

    Dmd _dmd;

    // Direction of input beam
    cpplap::Vect<float> _input_dir;

    // Angles of input beam
    float _input_phi;
    float _input_theta;

    float _wave_length;

    // Standard deviation for beam intensity measured by distance from center of
    // beam
    float _sig;

    // Spherical coordinates along X:
    // Number of pixels of image in phi direction
    int _Nphi = 100;

    // Number of pixels of image in theta direction
    int _Ntheta = 100;

    float _phi_min   = -M_PI / 4;
    float _phi_max   = M_PI / 4;
    float _theta_min = -M_PI / 4;
    float _theta_max = M_PI / 4;
    float _dphi      = (_phi_max - _phi_min) / _Nphi;
    float _dtheta    = (_theta_max - _theta_min) / _Ntheta;

    // cpplap::Vect<float>or containing all diffraction points on the dmd
    std::vector<DiffractionPoint> _diff_points;
    std::vector<float>            _image;

    std::vector<std::complex<float>> _ref_image_on;
    std::vector<std::complex<float>> _ref_image_off;

    // Random generator
    std::default_random_engine _generator;

    // Normal distribution with _sig as standard deviation
    std::normal_distribution<float> _normal_distribution;

    // Uniform distribution between 0 and 2*pi
    std::uniform_real_distribution<float> _uniform_distribution;

    friend std::ostream &operator<<(std::ostream &stream, const RayTracer &a);

    std::complex<float> computeOneMirrorPixel(const cpplap::Vect<float> out, const bool flipped)
    {
        float m  = _dmd.getMirrorSize();
        float ax = _input_dir[0];
        float ay = _input_dir[1];
        float az = _input_dir[2];
        float bx = out[0];
        float by = out[1];
        float bz = out[2];

        float s2 = sqrt(2);
        float ca = cos((flipped * 2 - 1) * _dmd.getMirrorAngle());
        float sa = sin((flipped * 2 - 1) * _dmd.getMirrorAngle());

        float r0 = _wave_length * _wave_length;
        float r1 = M_PI * M_PI;
        float r2 = ax + ay - bx - by + (ax - ay - bx + by) * ca - s2 * (az - bz) * sa;
        float r3 = -ax - ay + bx + by + (ax - ay - bx + by) * ca - s2 * (az - bz) * sa;
        float r  = r0 / r1 / r2 / r3;

        float argFactor = M_PI / _wave_length;
        float arg0      = 0;
        float arg1      = (2 * ax * m + 2 * ay * m - 2 * bx * m - 2 * by * m) * argFactor;
        float arg2 =
            (ax * m + ay * m - bx * m - by * m + (ax - ay - bx + by) * m * ca - s2 * (az - bz) * m * sa) * argFactor;
        float arg3 =
            (ax * m + ay * m - bx * m - by * m - (ax - ay - bx + by) * m * ca + s2 * (az - bz) * m * sa) * argFactor;

        float re0 = cos(arg0);
        float im0 = sin(arg0);
        float re1 = cos(arg1);
        float im1 = sin(arg1);
        float re2 = cos(arg2);
        float im2 = sin(arg2);
        float re3 = cos(arg3);
        float im3 = sin(arg3);

        float               nx = 1 / s2 * sa;
        float               ny = -nx;
        float               nz = sqrt(1 - sa * sa);
        cpplap::Vect<float> n(nx, ny, nz);
        float               intesityFactor = std::abs(_input_dir * n);

        float re = intesityFactor * r * (re0 + re1 - re2 - re3);
        float im = intesityFactor * r * (im0 + im1 - im2 - im3);

        return std::complex<float>(re, im);
    }

    // Assuming dir is normalized!
    float phaseShiftToReference(const cpplap::Vect<float> vect, const cpplap::Vect<float> dir) const
    {
        return vect * dir;
    }

public:
    RayTracer(const DmdParams params, const float input_phi, const float input_theta, const float wave_length,
              const float sig)
        : RayTracer(params, -1.0 * cpplap::Vect<float>::xyAngularCoords(input_phi, input_theta, 1), wave_length, sig)
    {
    }

    RayTracer(const DmdParams params, const cpplap::Vect<float> input_dir, const float wave_length, const float sig)
        : _dmd(params), _input_dir(input_dir), _input_phi(_input_dir.phi()), _input_theta(_input_dir.theta()),
          _wave_length(wave_length), _sig(sig), _normal_distribution(0, _sig), _uniform_distribution(0, 2 * M_PI)
    {
        _input_dir.normalize();
        _generator.seed(time(nullptr));
        _image.resize(_Nphi * _Ntheta);
    }

    cpplap::Vect<float> getRandomNormalVect()
    {
        float r   = _normal_distribution(_generator);
        float phi = _uniform_distribution(_generator);

        cpplap::Vect<float> ret = cpplap::Vect<float>::SphericalCoords(phi, M_PI / 2, r);

        ret.rotate(cpplap::Vect<float>(0, 1, 0), _input_theta);
        ret.rotate(cpplap::Vect<float>(0, 0, 1), _input_phi);

        return ret;
    }

    Ray getRandomNormalRay(const float amplitude = 1)
    {
        cpplap::Vect<float> ret = getRandomNormalVect();
        return Ray(ret, _input_dir, amplitude);
    }

    void calculateNormalDiffractionPoints(const int numb_points)
    {
        for (int i = 0; i < numb_points; i++) {
            _diff_points.push_back(_dmd.findDiffPoint(getRandomNormalRay()));
        }
    }

    int indexPixel(const int nphi, const int ntheta) const { return (ntheta + nphi * _Ntheta); }

    void setImageParams(const int Nphi = -500, const int Ntheta = -500, const float phi_min = -500,
                        const float phi_max = -500, const float theta_min = -500, const float theta_max = -500)
    {

        if (_ref_image_on.size() + _ref_image_off.size() > 1) {
            throw std::domain_error("Reference image has already been computed!");
        }

        // Get default value from default class variables
        if (Nphi != -500) {
            _Nphi = Nphi;
        }

        if (Ntheta != -500) {
            _Ntheta = Ntheta;
        }

        if (phi_min != -500) {
            _phi_min = phi_min;
        }

        if (phi_max != -500) {
            _phi_max = phi_max;
        }

        if (theta_min != -500) {
            _theta_min = theta_min;
        }

        if (theta_max != -500) {
            _theta_max = theta_max;
        }

        _dphi   = (_phi_max - _phi_min) / _Nphi;
        _dtheta = (_theta_max - _theta_min) / _Ntheta;

        _image.resize(_Nphi * _Ntheta);
    }

    void computeOneMirrorImages()
    {

        _ref_image_off.resize(_Nphi * _Ntheta);
        _ref_image_on.resize(_Nphi * _Ntheta);
        for (int i = 0; i < _Nphi; i++) {
            for (int j = 0; j < _Ntheta; j++) {
                float               phi          = _phi_min + i * _dphi;
                float               theta        = _theta_min + j * _dtheta;
                cpplap::Vect<float> output_dir   = cpplap::Vect<float>::xyAngularCoords(phi, theta, 1);
                _ref_image_off[indexPixel(i, j)] = computeOneMirrorPixel(output_dir, false);
                _ref_image_on[indexPixel(i, j)]  = computeOneMirrorPixel(output_dir, true);
            }
        }
    }

    void runAnalyticSimulation()
    {

        if (_ref_image_on.size() + _ref_image_off.size() < (size_t)2 * _Nphi * _Ntheta) {
            throw std::domain_error("Compute reference Image first");
        }

        const int Nthreads = std::thread::hardware_concurrency();

        std::thread threads[Nthreads];
        int         threadSize = _Nphi / (Nthreads - (_Nphi % Nthreads != 0));

        for (int i = 0; i < Nthreads; i++) {
            int start  = i * threadSize;
            int stop   = std::min((i + 1) * threadSize, _Nphi);
            threads[i] = std::thread([this, start, stop] { runPartialAnalyticSimulation(start, stop); });
        }
        for (std::thread &t : threads) {
            t.join();
        }
    }

    void runPartialAnalyticSimulation(int start, int stop)
    {
        cpplap::HessePlane<float> beam_plane(cpplap::Vect<float>(0, 0, 0), _input_dir);

        for (int i = start; i < stop; i++) {

            for (int j = 0; j < _Ntheta; j++) {

                float phi   = _phi_min + i * _dphi;
                float theta = _theta_min + j * _dtheta;

                cpplap::Vect<float> output_dir    = cpplap::Vect<float>::xyAngularCoords(phi, theta, 1);
                std::complex<float> tot_amplitude = 0;

                for (int k = 0; k < _dmd.nMirrors(); k++) {

                    const auto [flipped, x, y]      = _dmd.getMirrorStatus(k);
                    cpplap::Vect<float> mirror_vect = cpplap::Vect<float>(x, y, 0);

                    float tot_phase_shift = phaseShiftToReference(mirror_vect, _input_dir);

                    tot_phase_shift -= phaseShiftToReference(mirror_vect, output_dir);

                    cpplap::Line<float> beam_line(mirror_vect, _input_dir);
                    cpplap::Vect<float> beam_origin = beam_plane.intersectionPointWith(beam_line);

                    float amplitude = gauss(_sig, 0.0f, beam_origin.norm());

                    // Uncomment for non gaussian beam
                    // float amplitude = 1;

                    if (flipped) {
                        tot_amplitude +=
                            _ref_image_on[indexPixel(i, j)] *
                            std::polar(amplitude, static_cast<float>(2 * M_PI / _wave_length * tot_phase_shift));
                    }
                    else {
                        tot_amplitude +=
                            _ref_image_off[indexPixel(i, j)] *
                            std::polar(amplitude, static_cast<float>(2 * M_PI / _wave_length * tot_phase_shift));
                    }
                }
                _image[indexPixel(i, j)] = std::abs(tot_amplitude * tot_amplitude);
            }
        }
    }

    void runRayTracing()
    {

        _image.resize(_Nphi * _Ntheta);

        for (int i = 0; i < _Nphi; i++) {
            std::cout << i << " of " << _Nphi << std::endl;

#pragma omp parallel for
            for (int j = 0; j < _Ntheta; j++) {

                float               phi           = _phi_min + i * _dphi;
                float               theta         = _theta_min + j * _dtheta;
                std::complex<float> tot_amplitude = 0;
                cpplap::Vect<float> output_dir    = cpplap::Vect<float>::xyAngularCoords(phi, theta, 1);

                cpplap::HessePlane<float> image_plane(cpplap::Vect<float>(0, 0, 0), output_dir);

                for (DiffractionPoint diff_point : _diff_points) {

                    float tot_phase_shift = -image_plane.scalarDistanceTo(diff_point) + diff_point.getPhaseShift();

                    float amplitude = diff_point.getAmplitude();

                    tot_amplitude +=
                        std::polar(amplitude, static_cast<float>(2 * M_PI / _wave_length * tot_phase_shift));
                }
                _image[indexPixel(i, j)] = std::abs(tot_amplitude);
            }
        }
    }

    void printDiffPoints(std::ostream &stream) const
    {
        for (DiffractionPoint diff_point : _diff_points) {
            stream << diff_point[0] << " " << diff_point[1] << " " << diff_point[2] << " " << diff_point.getAmplitude()
                   << " " << diff_point.getPhaseShift() << std::endl;
        }
    }

    void normalizeImage()
    {
        float max = *std::max_element(_image.begin(), _image.end());
        for (auto &pixel : _image) {
            pixel /= max;
        }
    }

    void printImage(std::ostream &stream) const
    {

        for (int i = 0; i < _Nphi; i++) {
            for (int j = 0; j < _Ntheta; j++) {
                float phi   = _phi_min + i * _dphi;
                float theta = _theta_min + j * _dtheta;
                stream << phi << " " << theta << " " << _image[indexPixel(i, j)] << std::endl;
            }
        }
    }

    void printOneMirrorImage(std::ostream &stream, const bool flipped)
    {

        std::vector<std::complex<float>> &image = flipped ? _ref_image_on : _ref_image_off;

        if (image.size() < (size_t)_Nphi * _Ntheta) {
            throw std::domain_error("Compute reference Image first");
        }
        for (int i = 0; i < _Nphi; i++) {
            for (int j = 0; j < _Ntheta; j++) {
                float phi   = _phi_min + i * _dphi;
                float theta = _theta_min + j * _dtheta;
                stream << phi << " " << theta << " " << std::abs(image[indexPixel(i, j)]) << std::endl;
            }
        }
    }
};

std::ostream &operator<<(std::ostream &stream, const RayTracer &a)
{

    stream << "Dmd parameters: " << std::endl
           << a._dmd << "Beam input direction:" << std::endl
           << a._input_dir << "Wavelength " << a._wave_length << std::endl
           << "Beam sigma " << a._sig << std::endl
           << "Nphi " << a._Nphi << std::endl
           << "Ntheta " << a._Ntheta << std::endl
           << "phi_min " << a._phi_min << std::endl
           << "phi_max " << a._phi_max << std::endl
           << "theta_min " << a._theta_min << std::endl
           << "theta_max " << a._theta_max << std::endl;

    return stream;
}
#endif
