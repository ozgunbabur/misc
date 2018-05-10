package org.panda.misc.analyses;

import org.panda.utility.FileUtil;
import org.panda.utility.statistics.Summary;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

/**
 * @author Ozgun Babur
 */
public class SMMARTExpressionZScoreConvertor
{
	static final String TCGA_FILE = "/home/babur/Documents/TCGA/PanCan/tcga_Kallisto_tpm.gz";
	static final String TCGA_BRCA_FILE = "/home/babur/Documents/TCGA/BRCA/expression-tpm.gz";
	static final String S_101_M1 = "/home/babur/Documents/Analyses/SMMART/Patient1-revisit/data/met1-transcriptome.tsv";
	static final String S_101_M2 = "/home/babur/Documents/Analyses/SMMART/Patient1-revisit/data/met2-transcriptome.tsv";
	static final double LOG2 = Math.log(2);

	public static void main(String[] args) throws IOException
	{
		Map<String, double[]> dist = readDistributions(TCGA_BRCA_FILE, loadENSTToGeneMap(S_101_M2));

		String dir = "/home/babur/Documents/Analyses/SMMART/Patient1-revisit/data/";
		writeAsZScores(dist, dir + "rna-exp-log2.txt", new String[]{"Met1", "Met2"},
			readPatientExpr(dir + "met1-transcriptome.tsv"), readPatientExpr(dir + "met2-transcriptome.tsv"));
	}

	static Map<String, double[]> readDistributions(String filename, Map<String, String> enstToGene) throws IOException
	{
		BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(
			new FileInputStream(filename))));

		String line = reader.readLine(); // skip header
//		String[] header = line.split("\t");

		Map<String, double[]> geneSumMap = new HashMap<>();
		line = reader.readLine();
		while (line != null)
		{
			String[] t = line.split("\t");
			if (enstToGene.containsKey(t[0]))
			{
				String gene = enstToGene.get(t[0]);
				if (!geneSumMap.containsKey(gene))
				{
					geneSumMap.put(gene, new double[t.length-1]);
				}

				double[] cnt = geneSumMap.get(gene);
				for (int i = 1; i < t.length; i++)
				{
					double d = Double.valueOf(t[1]);
					d = Math.pow(2, d) - 0.001;
					cnt[i-1] += d;
				}
			}

			line = reader.readLine();
		}
		reader.close();

		Map<String, double[]> statMap = new HashMap<>();
		for (String gene : geneSumMap.keySet())
		{
			double[] sums = geneSumMap.get(gene);

			for (int i = 0; i < sums.length; i++)
			{
				sums[i] = Math.log(sums[i] + 0.001) / LOG2;
			}

			statMap.put(gene, new double[]{Summary.mean(sums), Summary.stdev(sums)});
		}
		return statMap;
	}

	static Map<String, String> loadENSTToGeneMap(String filename) throws IOException
	{
		return Files.lines(Paths.get(filename)).skip(1).map(l -> l.split("\t")[0].split("\\|")).collect(
			Collectors.toMap(t -> t[0], t -> t[5]));
	}

	static Map<String, Double> readPatientExpr(String filename) throws IOException
	{
		Map<String, Double> map = new HashMap<>();

		Files.lines(Paths.get(filename)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String gene = t[0].split("\\|")[5];
			Double v = Double.valueOf(t[t.length - 1]);
			if (map.containsKey(gene)) map.put(gene, map.get(gene) + v);
			else map.put(gene, v);
		});

		Map<String, Double> logged = new HashMap<>();
		map.keySet().forEach(gene -> logged.put(gene, Math.log(map.get(gene)) / LOG2));

		return logged;
	}

	static void writeAsZScores(Map<String, double[]> stats, String filename, String[] dataNames,
		Map<String, Double>... data) throws IOException
	{
		if (dataNames.length != data.length) throw new IllegalArgumentException(
			"data name length = " + dataNames.length + ", data length = " + data.length);

		Set<String> genes = Arrays.stream(data).map(Map::keySet).flatMap(Collection::stream).collect(Collectors.toSet());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));
		writer.write("Gene");
		Arrays.stream(dataNames).forEach(name -> FileUtil.tab_write(name, writer));

		genes.stream().sorted().forEach(gene ->
		{
			if (stats.containsKey(gene))
			{
				double[] s = stats.get(gene);

				if (s[1] == 0) return;

				FileUtil.lnwrite(gene, writer);

				Arrays.stream(data).forEach(map ->
				{
					FileUtil.tab_write("", writer);
					if (map.containsKey(gene)) FileUtil.write(Double.toString((map.get(gene) - s[0]) / s[1]), writer);
				});
			}
			else System.out.println(gene + " is not in TCGA.");
		});

		writer.close();
	}
}
