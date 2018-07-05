package org.panda.misc.analyses;

import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.*;
import org.panda.utility.statistics.trendline.PolyTrendLine;
import org.panda.utility.statistics.trendline.TrendLine;

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
	public static final String RPPAFile = "RPPA/Sample1/01_Gray_Johnson_Set131.txt";
	public static final String RPPAFile2 = "RPPA/Sample2/15_Joe_Gray__Brett_Johnson_Set140_111517.csv";
	public static final String RNA_FILE = "SMMART-101-RNA-seq-rawcounts.txt";
	public static final String GeneTrails_CNA_FILE = "irb16113_reported_copy_number_by_mrn_Subject_101_deID-CompBio.csv";

	public static void main(String[] args) throws IOException
	{
//		assessTumorSource();
//		printUniformityOfGraphSizePvals(RUNS_DIR);

		mapPlatforms();
		writeSamplesInZScores();

//		extractTCGABRCASubtype();

//		prepareCausalPathDirs();
//		Map<String, List<Double>> map = loadDistributions(DIR + "TCGA-BRCA-L4.csv");

//		int limit = 224;
//		BoxPlot.write(DIR + "plot.txt",
//			map.keySet().stream().limit(limit).collect(Collectors.toList()).toArray(new String[limit]),
//			map.keySet().stream().limit(limit).map(map::get).collect(Collectors.toList()).toArray(new List[limit]));

//		plotHistogram(map, "AKT_pT308", "EGFR_pY1068");
	}

	static void prepareCausalPathDirs() throws IOException
	{
		Map<String, List<Double>> map = loadDistributions(DIR + "TCGA-BRCA-L4.csv");
		Map<String, Map<String, Double>> samples = loadSampleMaps();

		for (String sample : samples.keySet())
		{
			Files.createDirectory(Paths.get(RUNS_DIR + sample));
			writeParametersFile(RUNS_DIR + sample + "/parameters.txt");
			writeSampleData(ZScore.get(map, samples.get(sample)), RUNS_DIR + sample + "/values.txt");
		}
	}

	static Map<String, List<Double>> loadDistributions(String filename) throws IOException
	{
		String[] header = Files.lines(Paths.get(filename)).findFirst().get().split(",");

		Map<String, List<Double>> map = new HashMap<>();

		for (int i = 4; i < header.length; i++)
		{
			map.put(header[i], new ArrayList<>());
		}


		Files.lines(Paths.get(filename)).skip(1).map(l -> l.split(",")).forEach(t ->
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

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(PATIENT_DIR + "RPPA/Sample1/platform-temp.txt"));
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

	static void writeSamplesInZScores() throws IOException
	{
		String compileDir = PATIENT_DIR + "RPPA/Sample2/";

		// Load TCGA distributions
		Map<String, List<Double>> distMap = loadDistributions(PATIENT_DIR + "RPPA/TCGA-BRCA-L4.csv");

		// Load Patient 101 RPPA data

		String idLine = Files.lines(Paths.get(PATIENT_DIR + RPPAFile)).skip(1).findFirst().get().trim();
		String valLine = Files.lines(Paths.get(PATIENT_DIR + RPPAFile)).skip(11).findFirst().get().trim();

		idLine = idLine.substring(idLine.indexOf("\t") + 1);
		valLine = valLine.substring(valLine.indexOf("Human Tissue") + 13);

		String[] ids = idLine.split("\t");
		Double[] values = Arrays.stream(valLine.split("\t")).map(s -> Double.valueOf(s)).collect(Collectors.toList()).toArray(new Double[0]);

		assert ids.length == values.length;

		// Load RPPA data in the same batch with the patient sample
		Double[][] xMatrix = new Double[10][];
		for (int i = 0; i < xMatrix.length; i++)
		{
			String s = Files.lines(Paths.get(PATIENT_DIR + RPPAFile)).skip(12 + i).findFirst().get().trim();
			s = s.substring(s.indexOf("Human Tissue") + 13);
			xMatrix[i] = Arrays.stream(s.split("\t")).map(ss -> Double.valueOf(ss)).collect(Collectors.toList()).toArray(new Double[0]);
		}

		// Read id mapping between the patient platform and TCGA platform
		Map<String, String> id2tcga = new HashMap<>();
		Files.lines(Paths.get(compileDir + "platform.txt")).skip(1).map(l -> l.split("\t")).filter(t -> t.length < 6 || !t[5].equals("new!")).forEach(t ->
			id2tcga.put(t[1], t[0]));
		Map<String, String> tcga2id = id2tcga.keySet().stream().collect(Collectors.toMap(id2tcga::get, id -> id));


		// Read TCGA Sample IDs in the patient batch
		List<String> sampleIDs = Files.lines(Paths.get(PATIENT_DIR + RPPAFile)).skip(12).map(l -> l.split("\t")[7])
			.collect(Collectors.toList());

		Set<String> sampleIDSet = new HashSet<>(sampleIDs);
		List<String> idsList = Arrays.asList(ids);
		Double[][] yMatrix = new Double[xMatrix.length][xMatrix[0].length];

		// Load corresponding RPPA data from TCGA
		String[] header = Files.lines(Paths.get(DIR + "TCGA-BRCA-L4.csv")).findFirst().get().split(",");
		Files.lines(Paths.get(DIR + "TCGA-BRCA-L4.csv")).skip(1).map(l -> l.split(",")).filter(t -> sampleIDSet.contains(t[0])).forEach(t ->
		{
			int sInd = sampleIDs.indexOf(t[0]);
			for (int i = 4; i < t.length; i++)
			{
				int aInd = idsList.indexOf(tcga2id.get(header[i]));
				if (aInd >= 0)
				{
					if (t[i].equals("NA")) t[i] = "NaN";
					Double v = Double.valueOf(t[i]);
					yMatrix[sInd][aInd] = v;
				}
			}
		});



		// Load patient metastasis #2

		idLine = Files.lines(Paths.get(PATIENT_DIR + RPPAFile2)).skip(1).findFirst().get().trim();
		valLine = Files.lines(Paths.get(PATIENT_DIR + RPPAFile2)).skip(11).findFirst().get().trim();

		idLine = idLine.substring(idLine.indexOf("\t") + 1);
		valLine = valLine.substring(valLine.indexOf("human tissue") + 13);

		String[] ids2 = idLine.split("\t");
		Double[] values2 = Arrays.stream(valLine.split("\t")).map(Double::valueOf).collect(Collectors.toList()).toArray(new Double[0]);

		assert ids2.length == values2.length;

		// Load RPPA data in the same batch with the patient sample
		Double[][] xMatrix2 = new Double[10][];
		for (int i = 0; i < xMatrix2.length; i++)
		{
			String s = Files.lines(Paths.get(PATIENT_DIR + RPPAFile2)).skip(31 + i).findFirst().get().trim();
			s = s.substring(s.indexOf("Breast Cancer") + 14);
			xMatrix2[i] = Arrays.stream(s.split("\t")).map(Double::valueOf).collect(Collectors.toList()).toArray(new Double[0]);
		}

		// print difference in two samples
		CollectionUtil.printVennSets(new HashSet<>(Arrays.asList(ids)), new HashSet<>(Arrays.asList(ids2)));
		if (true) System.exit(0);


//		// Print TCGA and SMMART RPPA Venn sets
//		Set<String> tcga_set = new HashSet<>(Arrays.asList(header).subList(4, Arrays.asList(header).size()));
//		Set<String> patient_set = new HashSet<>(idsList);
//		CollectionUtil.printNameMapping("Patient", "TCGA");
//		CollectionUtil.printVennSets(patient_set, tcga_set, tcga2id);

		for (String id : ids)
		{
			int abInd = idsList.indexOf(id);

			if (id2tcga.containsKey(id))
			{
				TrendLine trendLine = getTheTrendLine(xMatrix, yMatrix, abInd, id);
				values[abInd] = trendLine.predict(values[abInd]);
				distMap.put(id, distMap.get(id2tcga.get(id)));
			}
//			else
//			{
//				List<Double> list = new ArrayList<>();
//				for (Double[] row : xMatrix)
//				{
//					list.add(row[abInd]);
//				}
//				distMap.put(id, list);
//			}
		}

//		hist.plot();

		Map<String, Double> vals = new HashMap<>();
		for (int i = 0; i < ids.length; i++)
		{
			vals.put(ids[i], values[i]);
		}

		Map<String, Double> zScores = ZScore.get(distMap, vals);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(compileDir + "values.txt"));
		writer.write("ID\tValue");
		for (String id : zScores.keySet())
		{
			writer.write("\n" + id + "\t" + zScores.get(id));
		}
		writer.close();
	}

	static Histogram2D hist = new Histogram2D(0.01);

	static TrendLine getTheTrendLine(Double[][] xMatrix, Double[][] yMatrix, int abInd, String name)
	{
		double[] x = new double[xMatrix.length];
		double[] y = new double[x.length];

//		System.out.println("name = " + name);
		for (int i = 0; i < xMatrix.length; i++)
		{
			x[i] = xMatrix[i][abInd];
			y[i] = yMatrix[i][abInd];

//			System.out.println(xMatrix[i][abInd] + "\t" + yMatrix[i][abInd]);
		}

		Tuple corr = Correlation.pearson(x, y);

		if (corr.v < 0.5 || corr.p > 0.05)
		{
			double mX = Summary.mean(x);
			double mY = Summary.mean(y);

			hist.count(mX, mY);
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

	static void extractTCGABRCASubtype() throws IOException
	{
		String subtype = "Basal";

		Set<String> ids = Files.lines(Paths.get(PATIENT_DIR + "RPPA/Tumor-subtypes-for-PanCanPathways-PANCAN.txt"))
			.skip(1).map(l -> l.split("\t")).filter(t -> t[2].equals("BRCA") && t[3].equals(subtype)).map(t -> t[0])
			.collect(Collectors.toSet());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(PATIENT_DIR + "RPPA/TCGA-BRCA-" + subtype + ".csv"));
		writer.write(Files.lines(Paths.get(PATIENT_DIR + "RPPA/TCGA-BRCA-L4.csv")).findFirst().get());

		Files.lines(Paths.get(PATIENT_DIR + "RPPA/TCGA-BRCA-L4.csv")).filter(l -> ids.contains(l.substring(0, 12)))
			.forEach(l -> FileUtil.lnwrite(l, writer));

		writer.close();
	}

	/**
	 * Are they from Primary or Met?
	 */
	static void assessTumorSource() throws IOException
	{
		//read RNA

		double log2 = Math.log(2);
		double rThr = 2;
		double pThr = 2;

		Map<String, Double> rnaMap = new HashMap<>();

		Files.lines(Paths.get(PATIENT_DIR + RNA_FILE)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[1];
			double diff = (Math.log1p(Integer.valueOf(t[2])) - Math.log1p(Integer.valueOf(t[3]))) / log2;

			if (Math.abs(diff) > rThr) rnaMap.put(id, diff);
		});

//		rnaMap.keySet().forEach(k -> System.out.println(k + "\t" + rnaMap.get(k)));
//		System.out.println("\n");

		// evaluate GeneTrails CNA

		Files.lines(Paths.get(PATIENT_DIR + GeneTrails_CNA_FILE)).skip(1).limit(12).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[4];
			String status = t[8];

			if (rnaMap.containsKey(id))
			{
				System.out.println(id + "\t" + status + "\t" + rnaMap.get(id));
			}
		});

		System.out.println("\n");

		// evaluate RPPA

		//read platform
		Map<String, Set<String>> id2genes = new HashMap<>();

		Files.lines(Paths.get(PATIENT_DIR + "RPPA/Sample1/platform.txt")).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[1];
			if (id.contains("_p")) return;
			Set<String> genes = new HashSet<>(Arrays.asList(t[2].split(" ")));
			id2genes.put(id, genes);
		});

		Files.lines(Paths.get(PATIENT_DIR + "RPPA/Sample1/values.txt")).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String id = t[0];
			double val = Double.valueOf(t[1]);

			if (Math.abs(val) > pThr && id2genes.containsKey(id))
			{
				for (String gene : id2genes.get(id))
				{
					if (rnaMap.containsKey(gene))
					{
						System.out.println(gene + "\t" + val + "\t" + rnaMap.get(gene));
					}
				}
			}
		});
	}
}
