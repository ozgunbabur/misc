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
	static final double LOG2 = Math.log(2D);

	public static void main(String[] args) throws IOException
	{
		convertExternalFileComparingWithTCGA("/home/babur/Documents/TCGA/BRCA/expression.txt",
			new int[]{2, 3}, "/home/babur/Documents/Analyses/SMMART/Patient1/SMMART-101-RNA-seq-rawcounts.txt",
			1, "/home/babur/Documents/Analyses/SMMART/Patient1/SMMART-101-RNA-seq-Z-scores.txt");
	}

	public static void convertExternalFileComparingWithTCGA(String tcgaFile, int[] valIndex, String inFile,
		int symIndex, String outFile) throws IOException
	{
		ExpressionReader er = new ExpressionReader(tcgaFile);
		String[] samples = er.getSamples().toArray(new String[0]);

		Map<String, double[]> distMap = new HashMap<>();
		List<Map<String, Double>> valMaps = new ArrayList<>();
		for (int i : valIndex)
		{
			valMaps.add(new HashMap<>());
		}

		Files.lines(Paths.get(inFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String sym = t[symIndex];

			double[] dist = er.getGeneAlterationArray(sym, samples);
			if (dist != null) dist = removeZeros(dist);

			if (dist != null)
			{
				distMap.put(sym, dist);

				for (int i = 0; i < valIndex.length; i++)
				{
					double val = log(Double.valueOf(t[valIndex[i]]) + 1);
					valMaps.get(i).put(sym, val);
				}
			}
		});

		List<Map<String, Double>> zscores = new ArrayList<>();
		for (Map<String, Double> valMap : valMaps)
		{
			zscores.add(ZScore.get(distMap, valMap, null));
		}

		String[] header = Files.lines(Paths.get(inFile)).findFirst().get().split("\t");
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write(header[symIndex]);
		for (int ind : valIndex)
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

//	public static void convertTCGAFile(String tcgaFile, String outFile) throws IOException
//	{
//		ExpressionReader er = new ExpressionReader(tcgaFile);
//		String[] samples = er.getSamples().toArray(new String[0]);
//
//		Map<String, double[]> distMap = new HashMap<>();
//		List<Map<String, Double>> valMaps = new ArrayList<>();
//		for (int i : valIndex)
//		{
//			valMaps.add(new HashMap<>());
//		}
//
//		Files.lines(Paths.get(inFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
//		{
//			String sym = t[symIndex];
//
//			double[] dist = er.getGeneAlterationArray(sym, samples);
//			if (dist != null) dist = removeZeros(dist);
//
//			if (dist != null)
//			{
//				distMap.put(sym, dist);
//
//				for (int i = 0; i < valIndex.length; i++)
//				{
//					double val = log(Double.valueOf(t[valIndex[i]]) + 1);
//					valMaps.get(i).put(sym, val);
//				}
//			}
//		});
//
//		List<Map<String, Double>> zscores = new ArrayList<>();
//		for (Map<String, Double> valMap : valMaps)
//		{
//			zscores.add(ZScore.get(distMap, valMap, null));
//		}
//
//		String[] header = Files.lines(Paths.get(inFile)).findFirst().get().split("\t");
//		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
//		writer.write(header[symIndex]);
//		for (int ind : valIndex)
//		{
//			writer.write("\t" + header[ind]);
//		}
//
//		distMap.keySet().stream().sorted().forEach(gene ->
//		{
//			FileUtil.lnwrite(gene, writer);
//
//			for (Map<String, Double> map : zscores)
//			{
//				FileUtil.write("\t" + map.get(gene), writer);
//			}
//		});
//
//		writer.close();
//	}

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
