package org.panda.misc.analyses;

import org.panda.causalpath.run.CausalityAnalysisSingleMethodInterface;
import org.panda.resource.HGNC;
import org.panda.resource.PCPathway;
import org.panda.resource.network.PathwayCommons;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.Summary;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * For the preliminary data from Karin.
 */
public class PlatformComparison
{
	public static final String DIR = "/home/babur/Documents/Analyses/PlatformComparison/";
	final static String FORMATTED_DATA = DIR + "Data-for-causality.txt";
	final static String NORMALIZED_DATA = DIR + "Data-for-causality-normalized.txt";


	public static void main(String[] args) throws IOException
	{
//		printAverages();
//		writeFormattedData();
		doPathwayEnrichment("OCvsPT", 1);
		doPathwayEnrichment("PDXvsPT", 1);
		doPathwayEnrichment("OCvsPDX", 1);
	}

	private static Map<String, String> readIDMapping() throws IOException
	{
		return Files.lines(Paths.get(DIR + "synergizer.tsv")).skip(4).map(l -> l.split("\t")).filter(t -> t.length > 1)
			.collect(Collectors.toMap(t -> t[0], t -> HGNC.get().getSymbol(t[1])));
	}

	private static void writeFormattedData() throws IOException
	{
		Map<String, String> mapping = readIDMapping();
		String inFile = DIR + "human_crossTab.txt";
		String header = Files.lines(Paths.get(inFile)).findFirst().get();
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(FORMATTED_DATA));

		writer.write("ID\tSymbols\tSites\tEffect" + header);

		Files.lines(Paths.get(inFile)).map(l -> l.split("\t")).filter(t -> mapping.get(t[0]) != null)
			.forEach(t -> FileUtil.lnwrite(mapping.get(t[0]) + "\t" + mapping.get(t[0]) + "\t\t\t" +
				t[1] + "\t" + t[2] + "\t" + t[3] + "\t" + t[4] + "\t" + t[5] + "\t" + t[6], writer));

		writer.close();
	}

	private static void doPathwayEnrichment(String dir, double thr) throws IOException
	{
		Map<String, Double> vals = Files.lines(
			Paths.get(DIR + File.separator + dir + File.separator + "value-changes.txt")).skip(3)
			.map(l -> l.split("\t")).collect(Collectors.toMap(t -> t[0], t -> Double.valueOf(t[1])));

		Set<String> background = new HashSet<>(vals.keySet());

		Set<String> selected = vals.keySet().stream().filter(g -> Math.abs(vals.get(g)) >= thr)
			.collect(Collectors.toSet());

		PCPathway.get().writeEnrichmentResults(selected, background, 3, 300, DIR + File.separator + dir + "-pathway-enrichment.txt");
	}

	private static void printAverages() throws IOException
	{
		List<Double>[] lists = new List[6];
		for (int i = 0; i < 6; i++)
		{
			lists[i] = new ArrayList<>();
		}

		Files.lines(Paths.get(FORMATTED_DATA)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			for (int i = 4; i < 10; i++)
			{
				lists[i - 4].add(Double.parseDouble(t[i]));
			}
		});

		Map<List<Double>, Double> avgMap = new HashMap<>();
		for (int i = 0; i < 6; i++)
		{
			avgMap.put(lists[i], Summary.meanOfDoubles(lists[i]));
		}

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(NORMALIZED_DATA));

		writer.write(Files.lines(Paths.get(FORMATTED_DATA)).findFirst().get());

		Iterator<Double>[] iters = new Iterator[lists.length];
		for (int i = 0; i < 6; i++)
		{
			iters[i] = lists[i].iterator();
		}
		Files.lines(Paths.get(FORMATTED_DATA)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			FileUtil.lnwrite(ArrayUtil.getString("\t", t[0], t[1], t[2], t[3]), writer);
			for (int i = 4; i < 10; i++)
			{
				FileUtil.write("\t" + (iters[i-4].next() - avgMap.get(lists[i-4])), writer);
			}
		});

		writer.close();
	}
}
