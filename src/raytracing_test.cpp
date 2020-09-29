#include "raytracing.h"
#include "test.h"
#include <fstream>

int main()
{

    DmdParams dparams;
    dparams.nx          = 10;
    dparams.ny          = 10;
    dparams.angle       = M_PI * 2 / 3.;
    dparams.mirror_size = 1;
    dparams.mirror_dist = 0.2;

    cpplap::Vect<float> beam_direction = cpplap::Vect<float>::SphericalCoords(3, M_PI / 1.5, 1);
    RayTracer           traycer(dparams, beam_direction, 1, 0.5);

    cpplap::HessePlane<float> beam_plane(cpplap::Vect<float>(0, 0, 0), beam_direction);

    cpplap::test_exact(beam_plane.isInPlane(traycer.getRandomNormalVect()), true, "Normal random ray in plane");
}
