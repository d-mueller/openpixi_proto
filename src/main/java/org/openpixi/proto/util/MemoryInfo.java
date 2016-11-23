package org.openpixi.proto.util;

public class MemoryInfo {
    private static final long GIGABYTE = 1024L * 1024L * 1024L;

    public static void writeInfo() {
        Runtime runtime = Runtime.getRuntime();
        double memory = ((int) (100 * runtime.totalMemory() / GIGABYTE)) / 100.0;
        int mempercent = (int) (100 * runtime.totalMemory() / runtime.maxMemory());

        Logger.log("Memory    " + memory + "gb (" + mempercent + "%)");
    }
}
