package org.panda.misc.teaching;

import org.panda.utility.ArrayUtil;
import org.panda.utility.statistics.Histogram;

public class ABETHistogram
{
	static String[] gradeStr = ("18\n" +
		"15\n" +
		"0\n" +
		"14\n" +
		"15\n" +
		"20\n" +
		"5\n" +
		"15\n" +
		"15\n" +
		"15\n" +
		"8\n" +
		"5\n" +
		"5\n" +
		"0\n" +
		"15\n" +
		"18\n" +
		"5\n" +
		"14\n" +
		"20\n" +
		"20\n" +
		"18\n" +
		"0").split("\n");

	static double range = 5;

	public static void main(String[] args)
	{
		double[] grades = ArrayUtil.readToDoubleArrayPreserveMissing(gradeStr);
		Histogram h = new Histogram(range);
		h.setBorderAtZero(true);
		h.countAll(grades);

		h.print();
	}
}
