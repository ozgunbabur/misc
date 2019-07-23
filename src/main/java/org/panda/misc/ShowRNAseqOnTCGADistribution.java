package org.panda.misc;

import org.panda.resource.tcga.ExpressionReader;
import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.Summary;
import org.panda.utility.statistics.ZScore;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class ShowRNAseqOnTCGADistribution
{
	static final String TCGA_FILE = "/home/ozgun/Data/TCGA/BRCA/expression.txt";
	static final String IN_FILE = "/home/babur/Documents/Analyses/SMMART/Patient1/SMMART-101-RNA-seq-rawcounts.txt";
	static final double LOG2 = Math.log(2D);
	static final int SYM_INDEX = 1;
	static final int[] VAL_INDEX = new int[]{2, 3, 4};
	static String[] genes = new String[]{"ERBB3", "CELSR1", "PGR", "ICAM1", "ESR1", "ITGB6", "IL6ST", "AR", "GFRA3"};

	public static void main(String[] args) throws IOException
	{
		plotTCGADistribution();
//		plot();

//		CustomExpressionReader reader = loadTCGABRCATPM();
//		Set<String> samples = reader.getSamples();
//		System.out.println("samples.size() = " + samples.size());

//		check();
//		writePatientZScores();
	}

	public static void plotTCGADistribution() throws FileNotFoundException
	{
		HashSet<String> geneSet = new HashSet<>(Arrays.asList(genes));
		ExpressionReader er = new ExpressionReader(TCGA_FILE, geneSet);
		String[] samples = er.getSamples().toArray(new String[0]);

		for (String gene : genes)
		{
			System.out.println("\ngene = " + gene);
			Histogram h = new Histogram(1);
			h.setBorderAtZero(true);
			h.countAll(er.getGeneAlterationArray(gene, samples));
			h.print();
		}
	}

	public static void plot() throws IOException
	{
		HashSet<String> geneSet = new HashSet<>(Arrays.asList(genes));
		ExpressionReader er = new ExpressionReader(TCGA_FILE, geneSet);
		String[] samples = er.getSamples().toArray(new String[0]);

		List<Map<String, Double>> valMaps = new ArrayList<>();
		for (int i : VAL_INDEX)
		{
			valMaps.add(new HashMap<>());
		}

		Files.lines(Paths.get(IN_FILE)).skip(1).map(l -> l.split("\t")).filter(t -> geneSet.contains(t[SYM_INDEX]))
			.forEach(t ->
		{
			String sym = t[SYM_INDEX];

			for (int i = 0; i < VAL_INDEX.length; i++)
			{
				double val = log(Double.valueOf(t[VAL_INDEX[i]]) + 1);
				valMaps.get(i).put(sym, val);
			}
		});

		for (String gene : genes)
		{
			System.out.println("\n-------\n\ngene = " + gene);
			double[] dist = er.getGeneAlterationArray(gene, samples);

			Histogram h = new Histogram(1, dist);
			h.print();


			for (Map<String, Double> map : valMaps)
			{
				System.out.println("\n" + map.get(gene) + "\t" + 0);
				System.out.println(map.get(gene) + "\t" + 0.5);
			}
		}
	}

	private static double log(double x)
	{
		return Math.log1p(x) / LOG2;
	}



	//--Section: Converting patient RNAseq into Z-scores using TCGA RNAseq in TPM---------------------------------------

	private static void check() throws IOException
	{
		Map<String, Double> map = Files.lines(Paths.get("/home/babur/Documents/Analyses/SMMART/Patient115/RNAseq-zscores.txt"))
			.skip(1).map(l -> l.split("\t")).collect(Collectors.toMap(t -> t[0], t -> Double.valueOf(t[1])));

		Histogram h = new Histogram(1);
		for (String gene : map.keySet())
		{
			if (map.get(gene) < 100) h.count(map.get(gene));
		}
		h.print();
	}

	private static void writePatientZScores() throws IOException
	{
		Map<String, Double> patientT = readPatientTranscriptLevelRNASeq();
		Map<String, Set<String>> gToT = readPatientGeneToTMap();
		Map<String, double[]> tcgaT = readTCGABRCATranscripts(patientT.keySet());
		Map<String, double[]> tcgaG = sumUpTCGAToGeneLevel(tcgaT, gToT);
		Map<String, Double> patientG = sumUpPatientTToGeneLevel(patientT, gToT, tcgaT.keySet());
		Map<String, Double> zscores = convertPatientToZScore(patientG, tcgaG);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(
			"/home/babur/Documents/Analyses/SMMART/Patient115/RNAseq-zscores.txt"));

		writer.write("Gene\tZScore");

		for (String gene : zscores.keySet())
		{
			writer.write("\n" + gene + "\t" + zscores.get(gene));
		}

		writer.close();
	}

	private static Map<String, Double> readPatientTranscriptLevelRNASeq() throws IOException
	{
		return Files.lines(Paths.get("/home/babur/Documents/Analyses/SMMART/Patient115/abundance.tsv"))
			.skip(1).map(l -> l.split("\t")).filter(t -> t[0].endsWith("|protein_coding|"))
			.collect(Collectors.toMap(t -> t[0].split("\\|")[0], t -> Math.log(Double.valueOf(t[4]) + 1E-3) / LOG2));
	}

	private static Map<String, Set<String>> readPatientGeneToTMap() throws IOException
	{
		Map<String, Set<String>> map = new HashMap<>();
		Files.lines(Paths.get("/home/babur/Documents/Analyses/SMMART/Patient115/abundance.tsv"))
			.skip(1).map(l -> l.split("\t")[0])
			.filter(t -> t.endsWith("|protein_coding|"))
			.map(l -> l.split("\\|"))
			.forEach(t ->
			{
				String sym = t[5];
				String tid = t[0];

				if (!map.containsKey(sym)) map.put(sym, new HashSet<>());
				map.get(sym).add(tid);
			});

		return map;
	}

	private static Map<String, double[]> readTCGABRCATranscripts(Set<String> tids) throws IOException
	{
		return Files.lines(Paths.get("/home/babur/Documents/TCGA/BRCA/expression-tpm")).skip(1)
			.filter(l -> tids.contains(l.substring(0, l.indexOf("\t"))))
			.map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], ShowRNAseqOnTCGADistribution::readVals));
	}

	private static double[] readVals(String[] t)
	{
		double[] d = new double[t.length - 1];
		for (int i = 1; i < t.length; i++)
		{
			d[i - 1] = Double.valueOf(t[i]);
		}
		return d;
	}

	private static Map<String, double[]> sumUpTCGAToGeneLevel(Map<String, double[]> tcgaT,
		Map<String, Set<String>> gToTMap)
	{
		Map<String, double[]> map = new HashMap<>();

		for (String sym : gToTMap.keySet())
		{
			Set<String> tids = new HashSet<>(gToTMap.get(sym));
			tids.retainAll(tcgaT.keySet());

			if (tids.isEmpty()) continue;

			double[] v = new double[tcgaT.get(tids.iterator().next()).length];

			for (int i = 0; i < v.length; i++)
			{
				for (String tid : tids)
				{
					v[i] += tcgaT.get(tid)[i];
				}
			}

			map.put(sym, v);
		}

		return map;
	}

	private static Map<String, Double> sumUpPatientTToGeneLevel(Map<String, Double> patientT,
		Map<String, Set<String>> gToTMap, Set<String> tcgaTIDs)
	{
		Map<String, Double> map = new HashMap<>();

		for (String sym : gToTMap.keySet())
		{
			Set<String> tids = new HashSet<>(gToTMap.get(sym));
			tids.retainAll(tcgaTIDs);

			if (tids.isEmpty()) continue;

			double tot = 0;

			for (String tid : tids)
			{
				tot += patientT.get(tid);
			}

			map.put(sym, tot);
		}

		return map;
	}

	private static Map<String, Double> convertPatientToZScore(Map<String, Double> patientG, Map<String, double[]> tcgaG)
	{
		return ZScore.get(tcgaG, patientG, null);
	}
}
