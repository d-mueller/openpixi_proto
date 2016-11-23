package org.openpixi.proto.physics;

import org.openpixi.proto.math.SU2AlgebraElement;
import org.openpixi.proto.math.SU2GroupElement;

public class Grid {
    public int totalNumberOfCells;
    public int[] numCells;
    public Cell[] cells;
    public int acummulatedCellCount[];
    public double[] a;
    public double dt;

    public Grid(int[] numCells) {
        this.numCells = numCells;
        totalNumberOfCells = numCells[0] * numCells[1] * numCells[2];
        cells = new Cell[totalNumberOfCells];
        for (int i = 0; i < totalNumberOfCells; i++) {
            cells[i] = new Cell(i);
        }

        acummulatedCellCount = new int[3 + 1];
        acummulatedCellCount[3] = 1;
        for (int i = 3 - 1; i >= 0; i--) {
            acummulatedCellCount[i] = acummulatedCellCount[i + 1] * numCells[i];
        }
    }

    public int shift(int index, int direction, int orientation)
    {
        int result = index;
        int directionIndex = index / acummulatedCellCount[direction + 1];
        int withinDirectionIndex = directionIndex % numCells[direction];
        if (orientation > 0) {
            if (withinDirectionIndex == numCells[direction] - 1) {
                // wrap around along positive direction
                result -= acummulatedCellCount[direction];
            }
            result += acummulatedCellCount[direction + 1];
        } else if (orientation < 0) {
            if (withinDirectionIndex == 0) {
                // wrap around along negative direction
                result += acummulatedCellCount[direction];
            }
            result -= acummulatedCellCount[direction + 1];
        }
        return result;
    }

    public int[] getCellPos(int index)
    {
        int[] pos = new int[3];

        for(int i = 3-1; i >= 0; i--)
        {
            pos[i] = index % this.numCells[i];
            index -= pos[i];
            index /= this.numCells[i];
        }

        return pos;
    }

    public SU2GroupElement getStapleSum(int index, int d) {
        SU2GroupElement S = new SU2GroupElement(0, 0, 0, 0);
        int ci1 = shift(index, d, 1);
        int ci2, ci3, ci4;
        for (int i = 0; i < 3; i++) {
            if(i != d) {
                ci2 = shift(index, i, 1);
                ci3 = shift(ci1, i, -1);
                ci4 = shift(index, i, -1);
                //SU2GroupElement U1 = getU(ci1, i).mult(getU(ci2, d).adj());
                SU2GroupElement U1 = cells[ci1].U0[i].mult(cells[ci2].U0[d].adj());
                //U1.multAssign(getU(index, i).adj());
                U1.multAssign(cells[index].U0[i].adj());
                //SU2GroupElement U2 = getU(ci4, d).mult(getU(ci3, i));
                SU2GroupElement U2 = cells[ci4].U0[d].mult(cells[ci3].U0[i]);
                U2.adjAssign();
                //U2.multAssign(getU(ci4, i));
                U2.multAssign(cells[ci4].U0[i]);
                double areaFactor = 1.0 / Math.pow(a[i], 2);
                U1.addAssign(U2);
                U1.multAssign(areaFactor);
                S.addAssign(U1);
            }
        }
        return S;
    }

    public void evolve() {
        SU2AlgebraElement E;
        SU2GroupElement tempGroup = new SU2GroupElement();
        for (int i = 0; i < totalNumberOfCells; i++) {
            for (int d = 0; d < 3; d++) {
                E = cells[i].E[d];

                // Rot B
                tempGroup.set(cells[i].U0[d]);
                tempGroup.multAssign(getStapleSum(i, d));
                tempGroup.multAssign(dt);

                // Evolve E
                E.addAssign(tempGroup.proj());

                // Evolve U
                tempGroup.set(E.mult(-dt).getLink());
                tempGroup.multAssign(cells[i].U0[d]);
                cells[i].U1[d].set(tempGroup);
            }
        }
    }

    public void switchU() {
        for (int index = 0; index < totalNumberOfCells; index++) {
            cells[index].switchU();
        }
    }

    public SU2AlgebraElement getAvgB(int c0, int i) {
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;

        int c1 = shift(c0, j, 1);
        int c2 = shift(c0, k, 1);

        Cell cell0 = cells[c0];
        Cell cell1 = cells[c1];
        Cell cell2 = cells[c2];

        SU2GroupElement P1 = cell0.U0[j].copy();
        P1.multAssign(       cell1.U0[k]);
        P1.multAssignAdj(    cell2.U0[j]);
        P1.multAssignAdj(    cell2.U0[k]);
        SU2GroupElement P2 = cell0.U1[j].copy();
        P2.multAssign(       cell1.U1[k]);
        P2.multAssignAdj(    cell2.U1[j]);
        P2.multAssignAdj(    cell2.U1[k]);

        P1.addAssign(P2);
        P1.multAssign(0.5 * a[i] / (a[j] * a[k]));

        return P1.proj();
    }

    public SU2GroupElement getTP(int index, int direction) {
        SU2GroupElement tempGroup = cells[index].U1[direction].copy();
        tempGroup.multAssignAdj(cells[index].U0[direction]);
        return tempGroup;
    }
}
