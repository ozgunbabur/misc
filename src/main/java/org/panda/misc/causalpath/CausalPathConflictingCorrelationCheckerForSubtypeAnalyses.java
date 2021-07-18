package org.panda.misc.causalpath;

import org.panda.utility.Tuple;
import org.panda.utility.statistics.Correlation;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class CausalPathConflictingCorrelationCheckerForSubtypeAnalyses
{
	public static void main(String[] args) throws IOException
	{
		String base = "/Users/ozgun/Documents/Analyses/CausalPath-data/CPTAC-BRCA/";
		String dir = base + "subtypes/Basal-vs-others/";

		List<String[]> rels = readRelations(dir);
		Set<String> ids = new HashSet<>();
		rels.stream().forEach(s ->
		{
			ids.add(s[0]);
			ids.add(s[1]);
		});

		Set<String> cols = readTestCols(dir);
		Map<String, double[]> data = readData(base + "CPTAC-TCGA-BRCA-data-77.txt", cols, ids);

		Map<String, Double> pvals = new HashMap<>();
		for (String[] rel : rels)
		{
			double[] v0 = data.get(rel[0]);
			double[] v1 = data.get(rel[1]);

			Tuple tup = Correlation.pearson(v0, v1);
			System.out.println(Arrays.toString(rel) + "\t" + tup);
		}
	}

	private static List<String[]> readRelations(String dir) throws IOException
	{
		List<String[]> rows = new ArrayList<>();
		Files.lines(Paths.get(dir + "results.txt")).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String source = t[4];
			String target = t[7];

			String sign = Double.valueOf(t[5]) * Double.valueOf(t[8]) < 0 ? "-" : "+";

			rows.add(new String[]{source, target, sign});
		});

		return rows;
	}

	private static Set<String> readTestCols(String dir) throws IOException
	{
		return Files.lines(Paths.get(dir + "parameters.txt")).filter(l -> l.startsWith("test-value-column = "))
			.map(l -> l.substring(l.indexOf("= ") + 2)).collect(Collectors.toSet());
	}

	private static Map<String, double[]> readData(String file, Set<String> cols, Set<String> ids) throws IOException
	{
		Map<String, double[]> res = new HashMap<>();

		String[] header = Files.lines(Paths.get(file)).findFirst().get().split("\t");

		Files.lines(Paths.get(file)).skip(1).filter(l -> ids.contains(l.substring(0, l.indexOf("\t")))).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[0];
			int k = 0;
			double[] v = new double[cols.size()];

			for (int i = 0; i < header.length; i++)
			{
				if (cols.contains(header[i]))
				{
					v[k++] = t[i].equals("NA") ? Double.NaN : Double.valueOf(t[i]);
				}
			}

			res.put(id, v);
		});

		return res;
	}
}
