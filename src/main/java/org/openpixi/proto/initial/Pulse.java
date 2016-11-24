package org.openpixi.proto.initial;

import org.openpixi.proto.math.SU2AlgebraElement;
import org.openpixi.proto.math.SU2GroupElement;
import org.openpixi.proto.physics.Simulation;

public class Pulse {
    int start, orientation;
    double width, a0, a1, a2;

    public Pulse(int start, int orientation, double width, double a0, double a1, double a2) {
        this.start = start;
        this.orientation = orientation;
        this.width = width;
        this.a0 = a0;
        this.a1 = a1;
        this.a2 = a2;
    }

    public void initializeFields(Simulation simulation) {
        for (int i = 0; i < simulation.grid.totalNumberOfCells; i++) {
            int[] pos = simulation.grid.getCellPos(i);
            double x = (pos[0] - start) * simulation.aL;
            double f0 = 1.0 / (Math.sqrt(2.0 * Math.PI) * width) * Math.exp( - 0.5 * Math.pow(x / width, 2.0));
            double f1 = 1.0 / (Math.sqrt(2.0 * Math.PI) * width) * Math.exp( - 0.5 * Math.pow((x - simulation.dt * orientation) / width, 2.0));
            SU2AlgebraElement A0 = new SU2AlgebraElement(a0 * f0, a1 * f0, a2 * f0);
            SU2AlgebraElement A1 = new SU2AlgebraElement(a0 * f1, a1 * f1, a2 * f1);
            SU2GroupElement U0 = simulation.grid.cells[i].U0[1];
            SU2GroupElement U1 = simulation.grid.cells[i].U1[1];
            U0.multAssign(A0.getLink());
            U1.multAssign(A1.getLink());
        }

        for (int i = 0; i < simulation.grid.totalNumberOfCells; i++) {
            for (int j = 0; j < 3; j++) {
                SU2AlgebraElement E = simulation.grid.cells[i].E[j];
                E.set(simulation.grid.getTP(i, j).proj());
                E.multAssign(-1.0 / simulation.dt);
            }
        }
    }
}
