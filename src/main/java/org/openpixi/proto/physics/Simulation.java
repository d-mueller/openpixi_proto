package org.openpixi.proto.physics;

import org.openpixi.proto.initial.Pulse;
import org.openpixi.proto.observables.ProjectedEnergyDensity;
import org.openpixi.proto.util.Logger;
import org.openpixi.proto.util.MemoryInfo;
import org.openpixi.proto.util.PerformanceTimer;
import sun.rmi.runtime.Log;

public class Simulation {
    public int NL, NT;
    public double aL, aT, dt;
    public int iterations;

    public Grid grid;

    PerformanceTimer timer;

    public Simulation() {
        // Configuration
        NL = 512;
        NT = 1;
        aL = 0.1;
        aT = 1.0;
        dt = aL / 4.0;
        iterations = 100000;


        timer = new PerformanceTimer();
        timer.active = true;

        Logger.activated = false;
    }

    public void initialize() {
        timer.reset();

        MemoryInfo.writeInfo();

        // Initialization
        int[] numCells = new int[] {NL, NT, NT};
        grid = new Grid(numCells);
        grid.a = new double[] {aL, aT, aT};
        grid.dt   = dt;
        timer.lap("Init      ");

        MemoryInfo.writeInfo();

        timer.reset();
        // Initial conditions
        Pulse pulse = new Pulse(NL / 2, 1, 16.0 * aL, 1.0);
        pulse.initializeFields(this);
        timer.lap("Initial C ");

        MemoryInfo.writeInfo();
    }

    public void run() {
        ProjectedEnergyDensity diagnostic = new ProjectedEnergyDensity("pe1.dat", NL * aL);
        diagnostic.initialize(this);

        // Simulation loop
        for (int step = 0; step < iterations; step++) {
            Logger.log("SIM STEP   "+(step+1)+"/"+iterations);


            // Switch U0 and U1
            timer.reset();
            grid.switchU();
            timer.lap("Switch    ");

            timer.reset();
            // Evolve fields
            grid.evolve();
            timer.lap("Evolve    ");


            // Compute observables
            if(step % (aL / dt * NL) == 0) {
                timer.reset();
                diagnostic.evaluate(grid);
                timer.lap("Observe 1 ");
                diagnostic.writeToFile(this);
                timer.lap("Observe 2 ");

                MemoryInfo.writeInfo();
                Logger.logIgnore(step +" / " + iterations);
            }
        }
    }
}
