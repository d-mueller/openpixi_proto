package org.openpixi.proto.util;

public class Logger {
    public static boolean activated = true;

    public static void log(String s) {
        if(activated) System.out.println(s);
    }

    public static void logIgnore(String s) {
        System.out.println(s);
    }
}
