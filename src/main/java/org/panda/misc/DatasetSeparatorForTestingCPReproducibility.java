package org.panda.misc;

import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.FormatUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class DatasetSeparatorForTestingCPReproducibility
{
	public static void main(String[] args) throws IOException
	{
		String dir = "/home/babur/Documents/Analyses/CPTACBreastCancer77/correlation-based-phospho-0.1";
//		subset(dir, 0.5, 100);
		plotReproducibility(dir);
	}
	public static void subset(String inputDir, double ratio, int number) throws IOException
	{
		String[] t = inputDir.split("/");
		String d = t[t.length - 1];
		String s = inputDir.substring(0, inputDir.lastIndexOf("/"));
		String outDir = s + "/" + d + "-reproducibility";

		List<String> commonLines = Files.lines(Paths.get(inputDir + "/parameters.txt")).filter(l -> !l.startsWith("value-column = ")).collect(Collectors.toList());
		List<String> sampleLines = Files.lines(Paths.get(inputDir + "/parameters.txt")).filter(l ->  l.startsWith("value-column = ")).collect(Collectors.toList());

		int size = sampleLines.size();

		for (int i = 0; i < number; i++)
		{
			String dir = outDir + "/Subset-" + FormatUtil.posIntToString(i + 1, number);
			new File(dir).mkdirs();

			BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + "/parameters.txt"));
			commonLines.forEach(l -> FileUtil.writeln(l.replaceAll("\\.\\./", "../../").replaceAll("calculate-network-significance = true","calculate-network-significance = false"), writer));

			Collections.shuffle(sampleLines);

			for (int j = 0; j < size * ratio; j++)
			{
				writer.write(sampleLines.get(j) + "\n");
			}

			writer.close();
		}
	}

	public static void plotReproducibility(String dir) throws IOException
	{
		Set<String> original = Files.lines(Paths.get(dir + "/causative.sif")).map(l -> l.split("\t"))
			.filter(t -> t.length > 2).map(t -> ArrayUtil.getString("\t", t[0], t[1], t[2])).collect(Collectors.toSet());

		Map<String, Integer> relToCnt = new HashMap<>();
		original.forEach(rel -> relToCnt.put(rel, 0));

		File reproDir = new File(dir + "-reproducibility");
		for (File subDir : reproDir.listFiles())
		{
			if (!subDir.isDirectory()) continue;

			Set<String> subRels = Files.lines(Paths.get(subDir + "/causative.sif")).map(l -> l.split("\t"))
				.filter(t -> t.length > 2).map(t -> ArrayUtil.getString("\t", t[0], t[1], t[2])).collect(Collectors.toSet());

			subRels.forEach(rel ->
			{
				if (relToCnt.containsKey(rel)) relToCnt.put(rel, relToCnt.get(rel) + 1);
				else relToCnt.put(rel, 1);
			});
		}

		Map<Integer, Integer> origHist = new HashMap<>();
		Map<Integer, Integer> newHist = new HashMap<>();

		relToCnt.keySet().forEach(rel ->
		{
			Map<Integer, Integer> hist = original.contains(rel) ? origHist : newHist;
			int score = relToCnt.get(rel);
			if (hist.containsKey(score)) hist.put(score, hist.get(score) + 1);
			else hist.put(score, 1);
		});

		int max = relToCnt.values().stream().max(Integer::compare).get();
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(reproDir.getPath() + "/reproducibility-chart.txt"));
		writer.write("Reproduction count\tOriginal relation frequency\tNew relation frequency");
		for (int i = 0; i <= max; i++)
		{
			writer.write("\n" + i + "\t" + (origHist.containsKey(i) ? origHist.get(i) : 0) + "\t" + (newHist.containsKey(i) ? newHist.get(i) : 0));
		}

		writer.close();
	}
}
