package org.openpixi.proto.physics;

import org.openpixi.proto.initial.Pulse;
import org.openpixi.proto.observables.ProjectedEnergyDensity;
import org.openpixi.proto.util.Logger;
import org.openpixi.proto.util.MemoryInfo;
import org.openpixi.proto.util.PerformanceTimer;

public class Simulation {
    public int NL, NT;
    public double aL, aT, dt;
    public int iterations;
    public int observablesSteps;
    public int numberOfThreads;

    public Grid grid;

    PerformanceTimer timer;

    public Simulation() {
        // Configuration
        NL = 512;
        NT = 1;
        aL = 0.5;
        aT = 1.0;
        dt = aL / 4.0;
        iterations = 10000;
        observablesSteps = 10;
        numberOfThreads = 1;


        timer = new PerformanceTimer();
        timer.active = true;

        Logger.activated = true;
    }

    public void initialize() {
        timer.reset();

        MemoryInfo.writeInfo();

        // Initialization
        int[] numCells = new int[] {NL, NT, NT};
        grid = new Grid(numCells);
        grid.a = new double[] {aL, aT, aT};
        grid.dt   = dt;
        grid.numberOfThreads = numberOfThreads;
        timer.lap("Init      ");

        MemoryInfo.writeInfo();

        timer.reset();
        // Initial conditions
        Pulse pulse;
        pulse = new Pulse((int) (NL * 0.25), 1, 4 * aL, 1.0, 0.0, 0.0);
        pulse.initializeFields(this);
        pulse = new Pulse((int) (NL * 0.75), -1, 4 * aL, 0.0, 1.0, 0.0);
        pulse.initializeFields(this);
        timer.lap("Initial C ");

        MemoryInfo.writeInfo();
    }

    public void run() {
        ProjectedEnergyDensity diagnostic = new ProjectedEnergyDensity("pe1.dat", observablesSteps * dt);
        diagnostic.initialize(this);

        // Simulation loop
        for (int step = 0; step < iterations; step++) {
            Logger.log("SIM STEP   "+(step)+"/"+iterations);


            // Switch U0 and U1
            timer.reset();
            grid.switchU();
            timer.lap("Switch    ");

            timer.reset();
            // Evolve fields
            grid.evolve();
            timer.lap("Evolve    ");


            // Compute observables
            if(step % observablesSteps == 0) {
                timer.reset();
                diagnostic.evaluate(grid);
                timer.lap("Observe 1 ");
                diagnostic.writeToFile(this);
                timer.lap("Observe 2 ");

                MemoryInfo.writeInfo();
            }
        }
    }
}
