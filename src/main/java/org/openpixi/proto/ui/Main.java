package org.openpixi.proto.ui;

import org.openpixi.proto.physics.Simulation;

public class Main {
    public static void main(String[] args) {
        Simulation simulation = new Simulation();
        simulation.initialize();
        simulation.run();
    }
}
