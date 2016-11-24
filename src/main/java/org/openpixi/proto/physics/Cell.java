package org.openpixi.proto.physics;

import org.openpixi.proto.math.SU2AlgebraElement;
import org.openpixi.proto.math.SU2GroupElement;

public class Cell {
    int index;

    public SU2GroupElement[] U0 = new SU2GroupElement[3];
    public SU2GroupElement[] U1 = new SU2GroupElement[3];
    public SU2AlgebraElement[] E = new SU2AlgebraElement[3];

    public Cell(int index) {
        this.index = index;
        for (int i = 0; i < 3; i++) {
            U0[i] = new SU2GroupElement();
            U1[i] = new SU2GroupElement();
            E[i] = new SU2AlgebraElement();
        }
    }

    public void switchU() {
        /*
        for (int i = 0; i < 3; i++) {
            SU2GroupElement tmp = U0[i];
            U0[i] = U1[i];
            U1[i] = tmp;
        }
        */
        SU2GroupElement[] tmp = U0;
        U0 = U1;
        U1 = tmp;
    }
}
