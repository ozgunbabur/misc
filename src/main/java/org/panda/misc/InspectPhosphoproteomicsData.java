package org.panda.misc;

import org.panda.utility.ArrayUtil;
import org.panda.utility.statistics.*;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.DoubleStream;

/**
 * @author Ozgun Babur
 */
public class InspectPhosphoproteomicsData
{
	static final String DIR = "/home/babur/Documents/RPPA/TCGA/PNNL/";
	static final String DATA_FILE = DIR + "PNLL-causality-formatted.txt";
	static final String GROUPS_FILE = DIR + "subtypes/Subtype-Immunoreactive/parameters.txt";
	static final String PLOT = DIR + "plot.txt";
	static String[] control;
	static String[] test;


	static Map<String, Map<String, Double>> data;

	public static void main(String[] args) throws IOException
	{
		loadData();
//		replaceNansWithLowValues();
		loadGroups();
		writeDistribution(PLOT);
	}

	static void writeDistribution(String file) throws IOException
	{
		List<Double> pvals = new ArrayList<>();
		Map<String, Double> pMap = new HashMap<>();
		Map<String, Double> spMap = new HashMap<>();

		for (String id : data.keySet())
		{
			double[] c = getValues(id, control);
			double[] t = getValues(id, test);

			double p = getTTestSignedPValue(c, t);

			if (!Double.isNaN(p))
			{
				spMap.put(id, p);
				pvals.add(Math.abs(p));
				pMap.put(id, Math.abs(p));
			}
		}

		UniformityChecker.plot(pvals, file);

		Map<String, Double> qVals = FDR.getQVals(pMap, null);
		qVals.keySet().stream().filter(s -> qVals.get(s) <= 0.5)
			.sorted((s1, s2) -> qVals.get(s1).compareTo(qVals.get(s2))).forEach(id ->
			System.out.println(id + "\t" + spMap.get(id) + "\t" + qVals.get(id)));
	}

	static double logTransform(double p)
	{
		int sign = (int) Math.signum(p);
		p = Math.abs(p);
		return -Math.log(p) * sign;
	}

	static double getTTestSignedPValue(double[] control, double[] test)
	{
		double[] c = ArrayUtil.trimNaNs(control);
		double[] t = ArrayUtil.trimNaNs(test);

//		if (c.length < 13 || t.length < 13) return Double.NaN;

		double p = TTest.getPValOfMeanDifference(c, t);
		if (!Double.isNaN(p))
		{
			double dif = Summary.mean(t) - Summary.mean(c);
			if (dif < 0) p = -p;
		}
		return p;
	}

	static double getNanRatioDiff(double[] control, double[] test)
	{
		int tN = ArrayUtil.countNaNs(test);
		int cN = ArrayUtil.countNaNs(control);

		return (cN / (double) control.length) - (tN / (double) test.length);
	}

	static double getNanImbalanceSignedPValue(double[] control, double[] test)
	{
		int size = control.length + test.length;
		int tN = ArrayUtil.countNaNs(test);
		int cN = ArrayUtil.countNaNs(control);
		int featured = tN + cN;
		double p = FishersExactTest.calcImprovishmentPval(size, featured, test.length, tN);
		int sign = p > 0.5 ? -1 : 1;
		if (p > 0.5) p = 1-p;
		p *= 2 * sign;
		return p;
	}



	static double[] getValues(String id, String[] samples)
	{
		double[] v = new double[samples.length];
		for (int i = 0; i < v.length; i++)
		{
			v[i] = data.get(id).get(samples[i]);
		}
		return v;
	}

	static void loadGroups() throws IOException
	{
		Set<String> cSet = new HashSet<>();
		Set<String> tSet = new HashSet<>();
		Files.lines(Paths.get(GROUPS_FILE)).map(l -> l.split(" = ")).forEach(t ->
		{
			if (t[0].equals("control-value-column")) cSet.add(t[1]);
			if (t[0].equals("test-value-column")) tSet.add(t[1]);
		});

		Set<String> samples = getSamplesWithData();
		cSet.retainAll(samples);
		tSet.retainAll(samples);

		control = new ArrayList<>(cSet).toArray(new String[cSet.size()]);
		test = new ArrayList<>(tSet).toArray(new String[tSet.size()]);
	}

	static void loadData() throws IOException
	{
		data = new HashMap<>();
		Set<String> samples = getSamplesWithData();

		assert samples.contains("TCGA-13-1484");
		assert !samples.contains("TCGA-13-1483");
		assert samples.size() == 69;

		String[] header = readHeader();

		Files.lines(Paths.get(DATA_FILE)).skip(1).map(l -> l.split("\t")).filter(t -> !t[2].isEmpty()).forEach(t ->
		{
			String id = t[0];
			data.put(id, new HashMap<>());
			Map<String, Double> map = data.get(id);

			for (int i = 4; i < t.length; i++)
			{
				if (samples.contains(header[i]))
				{
					map.put(header[i], Double.valueOf(t[i]));
				}
			}
		});
	}
	
	static Set<String> getSamplesWithData() throws IOException
	{
		String[] header = readHeader();
		int[] cnt = new int[header.length];

		Files.lines(Paths.get(DATA_FILE)).skip(1).map(l -> l.split("\t")).filter(t -> !t[2].isEmpty()).forEach(t ->
		{
			for (int i = 4; i < cnt.length; i++)
			{
				if (t[i].equals("NaN")) cnt[i]++;
			}
		});

		int max = cnt[4];
		
		Set<String> samples = new HashSet<>();

		for (int i = 5; i < cnt.length; i++)
		{
			if (cnt[i] < max) samples.add(header[i]);
		}
		return samples;
	}

	private static String[] readHeader() throws IOException
	{
		return Files.lines(Paths.get(DATA_FILE)).findFirst().get().split("\t");
	}

	private static void replaceNansWithLowValues()
	{
		Random rand = new Random();

		for (String id : data.keySet())
		{
			Map<String, Double> map = data.get(id);

			double[] vals = map.values().stream().filter(d -> !d.isNaN()).flatMapToDouble(DoubleStream::of).toArray();
			double min = Summary.min(vals);
			double sd = Summary.stdev(vals);

			map.keySet().stream().filter(s -> map.get(s).isNaN()).forEach(s ->
				map.put(s, min + (rand.nextGaussian() * sd)));
		}
	}
}
