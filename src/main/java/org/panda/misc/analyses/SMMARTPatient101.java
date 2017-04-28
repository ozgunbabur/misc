package org.panda.misc.analyses;

import org.panda.utility.FileUtil;
import org.panda.utility.statistics.BoxPlot;
import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.Summary;
import org.panda.utility.statistics.UniformityChecker;

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
public class SMMARTPatient101
{
	public static final String DIR = "/home/babur/Documents/Analyses/BRCA-RPPA-from-MA/";
	public static final String RUNS_DIR = DIR + "sample-runs/";
	public static final String PATIENT_DIR = "/home/babur/Documents/Analyses/SMMART/Patient1/";

	public static void main(String[] args) throws IOException
	{
//		printUniformityOfGraphSizePvals(RUNS_DIR);

//		mapPlatforms();
		writeSample1InZScores();

//		prepareCausalPathDirs();
//		Map<String, List<Double>> map = loadDistributions();

//		int limit = 224;
//		BoxPlot.write(DIR + "plot.txt",
//			map.keySet().stream().limit(limit).collect(Collectors.toList()).toArray(new String[limit]),
//			map.keySet().stream().limit(limit).map(map::get).collect(Collectors.toList()).toArray(new List[limit]));

//		plotHistogram(map, "AKT_pT308", "EGFR_pY1068");
	}

	static void prepareCausalPathDirs() throws IOException
	{
		Map<String, List<Double>> map = loadDistributions();
		Map<String, Map<String, Double>> samples = loadSampleMaps();

		for (String sample : samples.keySet())
		{
			Files.createDirectory(Paths.get(RUNS_DIR + sample));
			writeParametersFile(RUNS_DIR + sample + "/parameters.txt");
			writeSampleData(getZScores(map, samples.get(sample)), RUNS_DIR + sample + "/values.txt");
		}
	}

	static Map<String, List<Double>> loadDistributions() throws IOException
	{
		String[] header = Files.lines(Paths.get(DIR + "TCGA-BRCA-L4.csv")).findFirst().get().split(",");

		Map<String, List<Double>> map = new HashMap<>();

		for (int i = 4; i < header.length; i++)
		{
			map.put(header[i], new ArrayList<>());
		}


		Files.lines(Paths.get(DIR + "TCGA-BRCA-L4.csv")).skip(1).map(l -> l.split(",")).forEach(t ->
		{
			for (int i = 4; i < t.length; i++)
			{
				if (!t[i].equals("NA")) map.get(header[i]).add(Double.valueOf(t[i]));
			}
		});

		map.values().forEach(Collections::sort);

		return map;
	}

	static Map<String, Map<String, Double>> loadSampleMaps() throws IOException
	{
		String[] header = Files.lines(Paths.get(DIR + "TCGA-BRCA-L4.csv")).findFirst().get().split(",");

		Map<String, Map<String, Double>> map = new HashMap<>();


		Files.lines(Paths.get(DIR + "TCGA-BRCA-L4.csv")).skip(1).map(l -> l.split(",")).forEach(t ->
		{
			String sampleID = t[0];
			Map<String, Double> sample = new HashMap<>();

			for (int i = 4; i < t.length; i++)
			{
				if (!t[i].equals("NA")) sample.put(header[i], Double.valueOf(t[i]));
			}

			map.put(sampleID, sample);
		});

		return map;
	}

	static void plotHistogram(Map<String, List<Double>> map, String... genes)
	{
		for (String gene : genes)
		{
			Histogram h = new Histogram(0.2);
			map.get(gene).stream().forEach(h::count);
			System.out.println("\ngene = " + gene);
			h.print();
		}
	}

	static Map<String, Double> getZScores(Map<String, List<Double>> distMap, Map<String, Double> sample)
	{
		Map<String, Double> z = new HashMap<>();

		for (String id : sample.keySet())
		{
			List<Double> dist = distMap.get(id);
			if (dist == null) continue;

			double mean = Summary.meanOfDoubles(dist);
			double sd = Summary.stdev(dist.toArray(new Double[dist.size()]));

			Double actual = sample.get(id);

			z.put(id, (actual - mean) / sd);
		}

		return z;
	}

	static void writeSampleData(Map<String, Double> sample, String filename) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));
		writer.write("ID\tValue");
		sample.keySet().stream().forEach(id -> FileUtil.lnwrite(id + "\t" + sample.get(id), writer));
		writer.close();
	}

	static void writeParametersFile(String filename) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));
		writer.write("proteomics-platform-file = ../../tcga-platform.txt");
		writer.write("\nproteomics-values-file = values.txt");
		writer.write("\nid-column = ID\n" +
			"symbols-column = Symbols\n" +
			"sites-column = Sites\n" +
			"effect-column = Effect\n");
		writer.write("\nthreshold-for-data-significance = 1\n" +
			"value-transformation = arithmetic-mean\n");

		writer.write("\ncalculate-network-significance = true\n" +
			"permutations-for-significance = 1000\n" +
			"color-saturation-value = 3\n" +
			"do-site-matching = true\n" +
			"value-column = Value");

		writer.close();
	}


	static void printUniformityOfGraphSizePvals(String dir) throws IOException
	{
		List<Double> vals = new ArrayList<>();
		for (File subDir : new File(dir).listFiles())
		{
			String filename = subDir.getPath() + "/significance-pvals.txt";
			if (new File(filename).exists())
			{
				Double pv = Double.valueOf(Files.lines(Paths.get(filename)).findFirst().get().split("=")[1].trim());
				vals.add(pv);
			}
		}
		UniformityChecker.plot(vals);
	}

	static void mapPlatforms() throws IOException
	{
		String idLine = Files.lines(Paths.get(PATIENT_DIR + "RPPA/Sample1/01_Gray_Johnson_Set131.txt")).skip(1).findFirst().get().trim();
		String geneLine = Files.lines(Paths.get(PATIENT_DIR + "RPPA/Sample1/01_Gray_Johnson_Set131.txt")).skip(3).findFirst().get().trim();

		String[] ids = idLine.substring(idLine.indexOf("\t") + 1).split("\t");
		String[] genes = geneLine.substring(geneLine.indexOf("\t") + 1).split("\t");

		for (int i = 0; i < genes.length; i++)
		{
			genes[i] = genes[i].split("/| ")[0];
		}

		System.out.println("ids.length = " + ids.length);

		Map<String, String> tcgaLinesMap = new HashMap<>();
		Map<String, Set<String>> tcgaGenesMap = new HashMap<>();

		Files.lines(Paths.get(PATIENT_DIR + "RPPA/tcga-platform.txt")).skip(1).forEach(l ->
		{
			String[] t = l.split("\t");
			tcgaLinesMap.put(t[0], l.substring(l.indexOf("\t") + 1));
			tcgaGenesMap.put(t[0], new HashSet<>(Arrays.asList(t[1].split(" "))));
		});

		int mapped = 0;

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(PATIENT_DIR + "RPPA/Sample1/platform.txt"));
		writer.write("ID-tcga\tID\tSymbols\tSites\tEffect");

		for (int i = 0; i < ids.length; i++)
		{
			String id = ids[i].toUpperCase().replaceAll("-", "").replaceAll("_P", "_p");
			if (Character.isDigit(id.charAt(0))) id = "X" + id;

			String[] t = id.split("_");
			if (t.length > 2)
			{
				id = t[0] + "_" + t[1];
				for (int j = 2; j < t.length; j++)
				{
					id += t[j];
				}
			}

			if (tcgaLinesMap.containsKey(id))
			{
				mapped++;
				writer.write("\n" + id + "\t" + ids[i] + "\t" + tcgaLinesMap.get(id));
			}
			else
			{
				Set<String> pids = new HashSet<>();
				for (String pid : tcgaGenesMap.keySet())
				{
					Set<String> set = tcgaGenesMap.get(pid);

					if (set.contains(genes[i]) &&
						((pid.contains("_") && id.contains("_")) || (!pid.contains("_") && !id.contains("_"))))
					{
						pids.add(pid);
					}
				}

				if (pids.size() == 1)
				{
					mapped++;
					writer.write("\n" + pids.iterator().next() + "\t" + ids[i] + "\t" + tcgaLinesMap.get(pids.iterator().next()));
				}
				else if (!pids.isEmpty())
				{
					writer.write("\n" + ids[i] + "\t" + ids[i] + "\t" + genes[i] + "\t\t" + pids);
				}
				else
				{
					writer.write("\n" + ids[i] + "\t" + ids[i] + "\t" + genes[i] + "\t\t\tnew!");
				}
			}
		}

		writer.close();
		System.out.println("mapped = " + mapped);
	}

	static void writeSample1InZScores() throws IOException
	{
		Map<String, List<Double>> distMap = loadDistributions();

		String idLine = Files.lines(Paths.get(PATIENT_DIR + "RPPA/Sample1/01_Gray_Johnson_Set131.txt")).skip(1).findFirst().get().trim();
		String valLine = Files.lines(Paths.get(PATIENT_DIR + "RPPA/Sample1/01_Gray_Johnson_Set131.txt")).skip(11).findFirst().get().trim();

		idLine = idLine.substring(idLine.indexOf("\t") + 1);
		valLine = valLine.substring(valLine.indexOf("Human Tissue") + 13);

		String[] ids = idLine.substring(idLine.indexOf("\t") + 1).split("\t");
		Double[] values = Arrays.stream(valLine.substring(valLine.indexOf("\t") + 1).split("\t")).map(s -> Double.valueOf(s)).collect(Collectors.toList()).toArray(new Double[0]);

		Map<String, String> id2tcga = new HashMap<>();
		Files.lines(Paths.get(PATIENT_DIR + "RPPA/Sample1/platform.txt")).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			id2tcga.put(t[1], t[0]);
		});

		Map<String, Double> vals = new HashMap<>();
		for (int i = 0; i < ids.length; i++)
		{
			vals.put(id2tcga.get(ids[i]), values[i]);
		}

		Map<String, Double> zScores = getZScores(distMap, vals);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(PATIENT_DIR + "RPPA/Sample1/values.txt"));
		writer.write("ID-tcga\tValue");
		for (String id : zScores.keySet())
		{
			writer.write("\n" + id + "\t" + zScores.get(id));
		}
		writer.close();
	}
}
