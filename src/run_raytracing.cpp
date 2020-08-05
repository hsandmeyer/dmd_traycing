#include "raytracing.h"
#include "test.h"
#include "units.h"
#include <cstdio>
#include <fstream>

using namespace cpplap::units;

int main()
{

    DmdParams dparams;

    /*
    Set up dmd parameters
    */

    dparams.nx          = 20;        // Number of mirrors in x-direction
    dparams.ny          = 20;        // Number of mirrors in y-direction
    dparams.mirror_size = 7250.0_nm; // Length of one edge of one mirror
    dparams.mirror_dist = 310.0_nm;  // Distance between mirrors
    dparams.angle       = 12_deg;    // Tilting angle of mirrors

    float     input_phi   = -21.0_deg; // Input angle in phi-direction
    float     input_theta = 21.0_deg;  // Input angle in theta-direction
    float     wave_length = 532.0_nm;  // Wavelength of the beam
    float     sig_spot    = 50.0_um;   // The intensity of the beam is modelled by a guass. This is sigma of this gauss
    RayTracer traycer     = RayTracer(dparams, input_phi, input_theta, wave_length, sig_spot);

    /*Print out parameters*/
    std::cout << traycer << std::endl;

    /*
    First, wie compute the diffraction image by raytraycing. To do so, we need to
    compute the diffraction points of the rays on the individual mirrors, first
    */
    traycer.calculateNormalDiffractionPoints(100000);

    /*
    Write diffraction points to file for late 3D plotting
    */
    std::ofstream out;
    out.open("diff_points.txt");

    if (out.is_open()) {
        traycer.printDiffPoints(out);
        out.close();
    }

    /*
    Setup image parameters
    (Number pixels in phi direction, number of pixels in theta direction,
    minimal outputh phi, maximal output phi, minimal output theta, maximal outputh theta)
    */
    traycer.setImageParams(100, 100, -10_deg, 10_deg, -10_deg, 10_deg);

    /*
    Run the raytracing. This will take some time
    */
    measure_time(std::cout, [&] { traycer.runRayTracing(); });
    /*
    Normalize the image
    */
    traycer.normalizeImage();

    /*
    Print out the calculated image for later plotting
    */
    out.open("ray_tracing.txt");
    if (out.is_open()) {
        traycer.printImage(out);
        out.close();
    }

    /*
    Now we simulate the same, using the analytic solution. First we need to compute
    the diffraction image of one mirror
    */
    traycer.computeOneMirrorImages();

    /*
    Print out the diffraction image for one mirror
    */
    out.open("one_mirror_ref_image.txt");
    if (out.is_open()) {
        traycer.printOneMirrorImage(out, false);
        out.close();
    }

    /*
    Run the analytic solution by phase shifting the image of each individual mirror
    */
    measure_time(std::cout, [&] { traycer.runAnalyticSimulation(); });

    /*
    Normalize the image
    */
    traycer.normalizeImage();

    /*
    Print out the analytic solution
    */
    out.open("analytic.txt");
    if (out.is_open()) {
        traycer.printImage(out);
        out.close();
    }
}
