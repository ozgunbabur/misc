package org.panda.misc;

import org.panda.resource.tcga.ExpressionReader;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.ZScore;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * This class takes an RNA-seq and a TCGA study to convert the first to Z-scores, using the distribution in the second.
 *
 * @author Ozgun Babur
 */
public class ConvertRNAseqToZScores
{
	static final String TCGA_FILE = "/home/babur/Documents/TCGA/BRCA/expression.txt";
	static final String IN_FILE = "/home/babur/Documents/Analyses/SMMART/Patient1/SMMART-101-RNA-seq-rawcounts.txt";
	static final String OUT_FILE = "/home/babur/Documents/Analyses/SMMART/Patient1/SMMART-101-RNA-seq-Z-scores.txt";
	static final double LOG2 = Math.log(2D);
	static final int SYM_INDEX = 1;
	static final int[] VAL_INDEX = new int[]{2, 3};

	public static void main(String[] args) throws IOException
	{
		convert();
	}

	public static void convert() throws IOException
	{
		ExpressionReader er = new ExpressionReader(TCGA_FILE);
		String[] samples = er.getSamples().toArray(new String[0]);

		Map<String, double[]> distMap = new HashMap<>();
		List<Map<String, Double>> valMaps = new ArrayList<>();
		for (int i : VAL_INDEX)
		{
			valMaps.add(new HashMap<>());
		}

		Files.lines(Paths.get(IN_FILE)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String sym = t[SYM_INDEX];

			double[] dist = er.getGeneAlterationArray(sym, samples);
			if (dist != null) dist = removeZeros(dist);

			if (dist != null)
			{
				distMap.put(sym, dist);

				for (int i = 0; i < VAL_INDEX.length; i++)
				{
					double val = log(Double.valueOf(t[VAL_INDEX[i]]) + 1);
					valMaps.get(i).put(sym, val);
				}
			}
		});

		List<Map<String, Double>> zscores = new ArrayList<>();
		for (Map<String, Double> valMap : valMaps)
		{
			zscores.add(ZScore.get(distMap, valMap, null));
		}

		String[] header = Files.lines(Paths.get(IN_FILE)).findFirst().get().split("\t");
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(OUT_FILE));
		writer.write(header[SYM_INDEX]);
		for (int ind : VAL_INDEX)
		{
			writer.write("\t" + header[ind]);
		}

		distMap.keySet().stream().sorted().forEach(gene ->
		{
			FileUtil.lnwrite(gene, writer);

			for (Map<String, Double> map : zscores)
			{
				FileUtil.write("\t" + map.get(gene), writer);
			}
		});

		writer.close();
	}

	private static double log(double x)
	{
		return Math.log1p(x) / LOG2;
	}

	private static double[] removeZeros(double[] arr)
	{
		List<Double> list = new ArrayList<>();
		for (double v : arr)
		{
			if (v > 0) list.add(v);
		}
		if (list.size() < 10) return null;

		double[] v = new double[list.size()];
		for (int i = 0; i < list.size(); i++)
		{
			v[i] = list.get(i);
		}
		return v;
	}
}
