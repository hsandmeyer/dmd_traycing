#include "dmd.h"
#include "test.h"

int main()
{
    DmdParams dparams;
    dparams.nx          = 10;
    dparams.ny          = 10;
    dparams.angle       = M_PI / 4;
    dparams.mirror_size = 1;
    dparams.mirror_dist = 0;

    Dmd dmd(dparams);

    cpplap::test_exact(dmd.index(3, 1), 13, "Index");
    cpplap::test_exact(dmd.up_x(13, 2), dmd.index(5, 1), "Index up x");
    cpplap::test_exact(dmd.up_y(13, -1), dmd.index(3, 0), "Index down y");

    cpplap::test_exact(dmd.index(cpplap::Vect<float>(-5, -5, 3)), dmd.index(0, 0), "Index of vector 1");
    cpplap::test_exact(dmd.index(cpplap::Vect<float>(0, 0, 3)), dmd.index(5, 5), "Index of vector 2");
    cpplap::test_exact(dmd.index(cpplap::Vect<float>(-5, 4.9, 3)), dmd.index(0, 9), "Index of vector 3");
    cpplap::test_exact(dmd.index(cpplap::Vect<float>(4.9, 4.9, 3)), dmd.index(9, 9), "Index of vector 4");

    dparams.mirror_size = 2;

    dmd = Dmd(dparams);
    cpplap::test_exact(dmd.index(cpplap::Vect<float>(9.9, 9.9, 0)), dmd.index(9, 9), "Index of vector 5");

    cpplap::test_exact(dmd.getMirror(0, 0).getR(), cpplap::Vect<float>(-9, -9, 0), "Get Mirror 1");
    cpplap::test_exact(dmd.getMirror(9, 0).getR(), cpplap::Vect<float>(9, -9, 0), "Get Mirror 2");
    cpplap::test_exact(dmd.getMirror(5, 5).getR(), cpplap::Vect<float>(1, 1, 0), "Get Mirror 3");

    Mirror mirror(cpplap::Vect<float>(-3, 1, 1), cpplap::Vect<float>(1, -2, -1), cpplap::Vect<float>(0, -1, 2));
    cpplap::Line<float> line = cpplap::Line<float>(cpplap::Vect<float>(2, -3, 2), cpplap::Vect<float>(1, -1, 3));
    cpplap::test_exact(mirror.intersectionPointWith(line), cpplap::Vect<float>(-1, 0, -7), "Mirror line intersec");

    // Vertical ray directly on mirror diagonal
    Ray ray(cpplap::Vect<float>(1, 1, 0), cpplap::Vect<float>(0, 0, 1), 1);

    cpplap::test_exact(dmd.findDiffPoint(ray), DiffractionPoint(1, 1, 0, 1), "Mirror Diffraction Point 1");

    // Vertical ray missing diagonal
    ray = Ray(cpplap::Vect<float>(0.5, 3.5, 0), cpplap::Vect<float>(0, 0, 1), 1);

    cpplap::test_rel(dmd.findDiffPoint(ray), DiffractionPoint(0.5, 3.5, -sqrt(0.5), 1, -sqrt(0.5)),
                     "Mirror Diffraction Point 2");

    // Ray in direction of normal vector of mirror plane pointing to mirror
    // diagonal
    ray = Ray(cpplap::Vect<float>(1, 1, 0), cpplap::Vect<float>(1, -1, 1), 1);

    cpplap::test_exact(dmd.findDiffPoint(ray), DiffractionPoint(1, 1, 0, 1), "Mirror Diffraction Point 3");

    // Ray in direction of negative normal vector of mirror plane pointing to
    // mirror diagonal
    ray = Ray(cpplap::Vect<float>(1, -3, 0), cpplap::Vect<float>(1, -1, 1), 1);

    cpplap::test_exact(dmd.findDiffPoint(ray), DiffractionPoint(1, -3, 0, 1), "Mirror Diffraction Point 4");

    // Ray in random direction missing mirror diagonal. Base vector of ray lies in
    // mirror.
    ray = Ray(cpplap::Vect<float>(3.5, -1.5, sqrt(0.5)), cpplap::Vect<float>(1, 2, 4), 2);

    cpplap::test_rel(dmd.findDiffPoint(ray), DiffractionPoint(3.5, -1.5, sqrt(0.5), 2), "Mirror Diffraction Point 5");

    // Ray that goes through three mirrors. We want to select the first one
    ray =
        Ray(cpplap::Vect<float>(1, 1, 0), cpplap::Vect<float>(2.5, -0.5, -sqrt(0.5)) - cpplap::Vect<float>(1, 1, 0), 2);

    cpplap::test_rel(dmd.findDiffPoint(ray),
                     DiffractionPoint(-0.5, 2.5, sqrt(0.5), 2, -sqrt(1.5 * 1.5 + 1.5 * 1.5 + 0.5)),
                     "Mirror Diffraction Point 6");

    dparams.mirror_dist = 2;
    dmd                 = Dmd(dparams);

    ray = Ray(cpplap::Vect<float>(3.3, 3.3, 5), cpplap::Vect<float>(0, 0, 1), 2);

    // Ray hitting area between mirrors.
    cpplap::test_rel(dmd.findDiffPoint(ray), DiffractionPoint(3.3, 3.3, 0, 0, -5), "Mirror Diffraction Point 7");
    return 0;
}
