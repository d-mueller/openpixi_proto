package org.openpixi.proto.observables;

import org.openpixi.proto.math.SU2AlgebraElement;
import org.openpixi.proto.physics.Grid;
import org.openpixi.proto.physics.Simulation;
import org.openpixi.proto.util.FileFunctions;

import java.io.*;

public class ProjectedEnergyDensity {
    private String path;
    private double timeInterval;
    private int stepInterval;
    private int maxWrites;
    private int writes;

    private int numberOfCells;
    private int longitudinalCells;
    private int[] numCells;
    private int[] longitudinalIndexArray;

    private double[] ET;
    private double[] BT;
    private double[] EL;
    private double[] BL;
    private double[] SL;
    private double[] JE;

    public ProjectedEnergyDensity(String path, double timeInterval) {
        this.path = path;
        this.timeInterval = timeInterval;
    }

    public void initialize(Simulation s) {
        this.stepInterval = (int) Math.max(Math.round(timeInterval / s.dt), 1);
        maxWrites = s.iterations / stepInterval;
        writes = 0;

        numCells = s.grid.numCells;
        numberOfCells = s.grid.totalNumberOfCells;
        longitudinalCells = numCells[0];

        initializeIndexArray(s.grid);

        ET = new double[longitudinalCells];
        BT = new double[longitudinalCells];
        EL = new double[longitudinalCells];
        BL = new double[longitudinalCells];
        SL = new double[longitudinalCells];
        JE = new double[longitudinalCells];

        reset();

        FileFunctions.clearFile(path);
        File file = FileFunctions.getFile(path);
        writeBinaryHeader(file, s);
    }

    public void writeToFile(Simulation simulation) {
        if(writes < maxWrites) {
            finalizeArrays(simulation.grid);

            // Write to file
            File file = FileFunctions.getFile(path);

            try {
                DataOutputStream stream = null;
                try {
                    stream = getStream(file);
                    double time = 0.0;
                    stream.writeDouble(time);
                    writeBinaryDoubleArray(stream, ET);
                    writeBinaryDoubleArray(stream, BT);
                    writeBinaryDoubleArray(stream, EL);
                    writeBinaryDoubleArray(stream, BL);
                    writeBinaryDoubleArray(stream, SL);
                    writeBinaryDoubleArray(stream, JE);
                } finally {
                    stream.flush();
                    stream.close();
                    writes++;
                }
            } catch (IOException ex) {
                System.out.println("ProjectedEnergyDensity2: Error writing to file.");
            }
        }
    }

    public void reset() {
        for (int i = 0; i < longitudinalCells; i++) {
            ET[i] = 0.0;
            BT[i] = 0.0;
            EL[i] = 0.0;
            BL[i] = 0.0;
            SL[i] = 0.0;
            JE[i] = 0.0;
        }
    }

    public void finalizeArrays(Grid grid) {
        multiplyArraysWithFactors(grid);
        shiftFields();
    }

    public void evaluate(Grid grid) {
        for (int index = 0; index < grid.totalNumberOfCells; index++) {

            // Use index array to get longitudinal index.
            int lindex = longitudinalIndexArray[index];

            // Field components
            int iShiftX = grid.shift(index, 0, 1);

            SU2AlgebraElement Ex =  grid.cells[index  ].E[0];
            SU2AlgebraElement Ey0 = grid.cells[index  ].E[1];
            SU2AlgebraElement Ey1 = grid.cells[iShiftX].E[1];
            SU2AlgebraElement Ez0 = grid.cells[index  ].E[2];
            SU2AlgebraElement Ez1 = grid.cells[iShiftX].E[2];

            // Temporally averaged magnetic fields
            SU2AlgebraElement By = grid.getAvgB(index, 1);
            SU2AlgebraElement Bz = grid.getAvgB(index, 2);

            // Squared field components components
            double ExSq = Ex.square();
            double EySq = Ey0.square();
            double EzSq = Ez0.square();

            double BxSq = grid.getAvgB(index, 0).square();
            double BySq = By.square();
            double BzSq = Bz.square();

            // Poynting vector

            // Spatially averaged electric fields
            SU2AlgebraElement Ey = Ey0.add(Ey1);
            SU2AlgebraElement Ez = Ez0.add(Ez1);


            // Two parts of the Poynting vector SL1 and SL2.
            // These are not defined at the same positions in the transverse plane, but since we average over
            // the transverse coordinates it doesn't matter.
            double SL1 = - Ez.mult(By);
            double SL2 = Ey.mult(Bz);

            // Power input (not correctly time-averaged, this is problematic. How to fix this without buffering?)
            //AlgebraElement jx = grid.getJ(index, 0);
            //double jInE = jx.mult(Ex);
            double jInE = 0.0;

            // Synchronized write to arrays
            synchronized (this) {
                ET[lindex] += EySq + EzSq;
                BT[lindex] += BySq + BzSq;
                EL[lindex] += ExSq;
                BL[lindex] += BxSq;
                SL[lindex] += SL1 + SL2;
                JE[lindex] += jInE;
            }
        }
    }

    private void initializeIndexArray(Grid grid) {
        longitudinalIndexArray = new int[numberOfCells];

        for (int i = 0; i < grid.totalNumberOfCells; i++) {
            int[] gridPos = grid.getCellPos(i);
            longitudinalIndexArray[i] = gridPos[0];
        }
    }

    private void multiplyArraysWithFactors(Grid grid) {
        double cellFactor = grid.totalNumberOfCells / ((double) grid.numCells[0]);
        double unitFactorLSquared = Math.pow(grid.a[0], -2) / cellFactor;
        double unitFactorTSquared = Math.pow(grid.a[1], -2) / cellFactor;
        double unitFactorL = 1.0 / (cellFactor * grid.a[0]);

        for (int i = 0; i < longitudinalCells; i++) {
            ET[i] *= 0.5 * unitFactorTSquared;
            BT[i] *= 0.5 * unitFactorTSquared;
            EL[i] *= 0.5 * unitFactorLSquared;
            BL[i] *= 0.5 * unitFactorLSquared;
            SL[i] *= 0.5 * unitFactorTSquared;
            JE[i] *= unitFactorL;
        }
    }

    private void shiftFields() {
        /*
        Before shifting:
        JE, SL, EL, BT are half-shifted
        ET, BL are unshifted

        After shifting:
        All quantities are defined at unshifted lattice points
        */
        shiftArray(EL);
        shiftArray(BT);
        shiftArray(SL);
        shiftArray(JE);
    }

    private void shiftArray(double[] array) {
        double[] source = array.clone();
        int n = array.length;

        for (int i = 1; i < n; i++) {
            array[i] = (source[i] + source[i - 1]) * 0.5;
        }

        array[0] = (source[n-1] + source[0]) * 0.5;
    }

    private void writeBinaryHeader(File file, Simulation s) {
        int longitudinalCells = s.grid.numCells[0];
        try {
            DataOutputStream stream = null;
            try {
                stream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file, true)));
                stream.writeInt(longitudinalCells);
                stream.writeInt(maxWrites);
            } finally {
                stream.flush();
                stream.close();
            }
        } catch (IOException ex) {
            System.out.println("ProjectedEnergyDensity: Error writing to file.");
        }
    }

    private DataOutputStream getStream(File file) throws IOException {
        return new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file, true)));
    }

    private void writeBinaryDoubleArray(DataOutputStream stream, double[] array) throws IOException {
        for (int i = 0; i < array.length; i++) {
            stream.writeDouble(array[i]);
        }
    }
}
