package org.panda.misc.analyses;

import org.panda.utility.FileUtil;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.Correlation;
import org.panda.utility.statistics.Summary;
import org.panda.utility.statistics.trendline.PolyTrendLine;
import org.panda.utility.statistics.trendline.TrendLine;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class SMMARTRPPANormalizer
{
	public static Map<String, Map<String, Double>> convertToZScores(Map<String, Map<String, Double>> data)
	{
		Map<String, Map<String, Double>> zData = new HashMap<>();

		for (String ab : data.keySet())
		{
			double[] x = new double[data.get(ab).keySet().size()];
			int i = 0;
			for (String sample : data.get(ab).keySet())
			{
				x[i++] = data.get(ab).get(sample);
			}

			double mean = Summary.mean(x);
			double stdev = Summary.stdev(x);

			zData.put(ab, new HashMap<>());

			for (String sample : data.get(ab).keySet())
			{
				zData.get(ab).put(sample, (data.get(ab).get(sample) - mean) / stdev);
			}
		}
		return zData;
	}

	public static void writeData(Map<String, Map<String, Double>> data, String filename) throws IOException
	{
		List<String> samples = data.keySet().stream().map(data::get).map(Map::keySet).flatMap(Collection::stream)
			.distinct().sorted().collect(Collectors.toList());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));
		writer.write("ID");
		samples.forEach(s -> FileUtil.tab_write(s, writer));

		data.keySet().stream().sorted().forEach(ab ->
		{
			FileUtil.lnwrite(ab, writer);

			samples.forEach(sample ->
				FileUtil.tab_write(data.get(ab).containsKey(sample) ? data.get(ab).get(sample).toString() : "NaN",
					writer));
		});

		writer.close();
	}

	/**
	 * The first dataset is modified and returned.
	 */
	public static Map<String, Map<String, Double>> normalizeAndUnite(Map<String, Map<String, Double>>... maps)
	{
		Map<String, Map<String, Double>> base = maps[0];

		for (int i = 1; i < maps.length; i++)
		{
			Map<String, Map<String, Double>> addMap = maps[i];
			normalizeAndAddOne(base, addMap);
		}

		return base;
	}

	private static void normalizeAndAddOne(Map<String, Map<String, Double>> base, Map<String, Map<String, Double>> add)
	{
		for (String ab : add.keySet())
		{
			if (base.containsKey(ab))
			{
				// focus on the data of the id
				Map<String, Double> map0 = base.get(ab);
				Map<String, Double> map1 = add.get(ab);

				// find intersection
				Set<String> common = new HashSet<>(map0.keySet());
				common.retainAll(map1.keySet());

				// find the scaling parameters
				double[] x = new double[common.size()];
				double[] y = new double[common.size()];
				int i = 0;
				for (String sample : common)
				{
					x[i] = map1.get(sample);
					y[i] = map0.get(sample);
					i++;
				}
				TrendLine trendLine = getTheTrendLine(x, y);

				// add the new set to the base
				Set<String> diff1 = new HashSet<>(map1.keySet());
				diff1.removeAll(common);
				diff1.forEach(sample -> map0.put(sample, trendLine.predict(map1.get(sample))));
			}
			else
			{
				base.put(ab, add.get(ab));
			}
		}
	}


	static TrendLine getTheTrendLine(double[] x, double[] y)
	{
		Tuple corr = Correlation.pearson(x, y);

		if (corr.v < 0.5 || corr.p > 0.05)
		{
			double mX = Summary.mean(x);
			double mY = Summary.mean(y);

			return new TrendLine()
			{
				@Override
				public void setValues(double[] y, double[] x)
				{

				}

				@Override
				public double predict(double x)
				{
					return x + (mY - mX);
				}
			};
		}
		else
		{
			PolyTrendLine ptl = new PolyTrendLine(1);
			ptl.setValues(y, x);
			return ptl;
		}
	}
}
