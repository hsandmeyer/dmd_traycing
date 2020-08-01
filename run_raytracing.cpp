#include "raytracing.h"
//#include "raytracing_omp.h"
#include "test.h"
#include "units.h"
#include <cstdio>
#include <fstream>

int main()
{

    DmdParams dparams;

    dparams.nx          = 20;
    dparams.ny          = 20;
    dparams.mirror_size = 7250.0_nm;
    dparams.mirror_dist = 310.0_nm;
    dparams.angle       = 12_deg;
    // dparams.angle = 0;

    // for 2 x 2
    // RayTracer traycer = RayTracer(dparams, -21_deg, 21_deg, 532.0_nm, 3.0_um);

    // for 10 x 10
    // RayTracer traycer = RayTracer(dparams, -21_deg, 21_deg, 532.0_nm, 16.0_um);

    // for 100 x 100
    // RayTracer traycer = RayTracer(dparams, 0, 0, 5.32e-7, 20e-5);
    RayTracer traycer = RayTracer(dparams, -21_deg, 21_deg, 532.0_nm, 50.0_um);

    std::cout << traycer << std::endl;

    traycer.calculateNormalDiffractionPoints(100000);

    std::ofstream out;
    out.open("diff_points.txt");

    if (out.is_open()) {
        traycer.printDiffPoints(out);
        out.close();
    }

    traycer.setImageParams(100, 100, -10_deg, 10_deg, -10_deg, 10_deg);

    traycer.runRayTracing();
    traycer.normalizeImage();

    out.open("image.txt");
    if (out.is_open()) {
        traycer.printImage(out);
        out.close();
    }

    traycer.setImageParams(100, 100);
    traycer.computeReferenceImages();

    out.open("one_mirror_ref_image.txt");
    if (out.is_open()) {
        traycer.printReferenceImage(out, false);
        out.close();
    }

    // measure_time(std::cout, printf, "Hallo %d", 5);
    measure_time(std::cout, [&] { traycer.runReferenceTraycing(); });

    traycer.normalizeImage();

    out.open("ref_image.txt");
    if (out.is_open()) {
        traycer.printImage(out);
        out.close();
    }
}
