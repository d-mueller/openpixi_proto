package org.openpixi.proto.util;

public class PerformanceTimer {

	private long time = 0;
	public boolean active = false;

	public void reset() {
		time = System.nanoTime();
	}

	public void lap(String s) {
		if(active) {
			long delta = (System.nanoTime() - time) / (1000 * 1000);
			Logger.log(s.substring(0, 10) + " " + delta + "ms");
			reset();
		}
	}

}
