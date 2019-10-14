package org.panda.misc.analyses;

import org.panda.resource.proteomics.RPPAIDMapper;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.Summary;
import org.panda.utility.statistics.TTest;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class HillPaper
{
	public static void main(String[] args) throws IOException
	{
		String[] cellLine = new String[]{"BT20", "BT549", "MCF7", "UACC812"};

		for (int i = 0; i < cellLine.length; i++)
		{
//			convertData(cellLine[i]);
//			generateAnalysisFolders(cellLine[i]);
//			generateInhibitorAnalysisFolders(cellLine[i]);
		}

//		writeAllExpectedChanges("/home/ozgun/Analyses/HillPaper/LigandEffects", "/home/ozgun/Analyses/HillPaper/expected.txt");

		testHypotheses();
	}

	public static void convertData(String cellLine) throws IOException
	{
		String inDir = "/home/ozgun/Analyses/HillPaper/DataS1/complete/";
		String outDir = "/home/ozgun/Analyses/HillPaper/DataS1/converted/";

		Map<String, Map<String, Double>> map = load(inDir + cellLine + ".csv");
		RPPAIDMapper.writeAsCPFile(map, outDir + cellLine + ".txt",
			map.values().iterator().next().keySet().stream().sorted().collect(Collectors.toList()));
	}


	public static Map<String, Map<String, Double>> load(String filename) throws IOException
	{
		String[] abs = Files.lines(Paths.get(filename)).filter(l -> l.contains(",Antibody Name,")).findFirst().get().split(",");
		for (int i = 0; i < abs.length; i++)
		{
			abs[i] = abs[i].replaceAll(" ", "-");
		}

		int startCol = Arrays.asList(abs).indexOf("Antibody-Name") + 1;

		boolean halt = false;
		// check if there is missing ab
		for (int i = startCol; i < abs.length; i++)
		{
			String preferredID = RPPAIDMapper.get().getPreferredID(abs[i]);
			if (preferredID == null)
			{
				halt = true;
				System.err.println("Not found: " + abs[i]);
			}
		}

		if (halt)
		{
			System.out.println();
			throw new RuntimeException("There we unrecognized antibodies. Add to the library.");
		}
		//--------

		Map<String, Map<String, Double>> map = new HashMap<>();

		BufferedReader reader = new BufferedReader(new FileReader(filename));
		for (int i = 0; i < 5; i++)
		{
			reader.readLine();
		}

		int j = 0;
		String prevSamp = null;

		for (String line = reader.readLine(); line != null; line = reader.readLine())
		{
			String[] t = line.split(",");

			String currSamp = t[1] + "-" + t[2] + "-" + t[3];

			if (t[1].isEmpty())
			{
				if (!currSamp.equals(prevSamp))
				{
					j = 0;
					prevSamp = currSamp;
				}

				currSamp += "-rep" + (++j);
			}

			for (int i = startCol; i < t.length; i++)
			{
				String ab = RPPAIDMapper.get().getPreferredID(abs[i]);
				Double val = Double.valueOf(t[i]);

				if (!map.containsKey(ab)) map.put(ab, new HashMap<>());
				map.get(ab).put(currSamp, val);
			}
		}

		return map;
	}

	private static void generateAnalysisFolders(String cellLineName) throws IOException
	{
		String[] header = Files.lines(Paths.get("/home/ozgun/Analyses/HillPaper/DataS1/converted/" +
			cellLineName + ".txt")).findFirst().get().split("\t");

		Set<String> controls = new HashSet<>();
		Map<String, Set<String>> ligandToSamples = new HashMap<>();
		Set<String> inhibitors = new HashSet<>();

		for (int i = 4; i < header.length; i++)
		{
			String[] t = header[i].split("-");

			if (t[0].equals("Serum") && t[1].equals("DMSO")) controls.add(header[i]);

			if (!t[0].isEmpty() && !t[0].equals("Serum") && t[1].equals("DMSO"))
			{
				if (!ligandToSamples.containsKey(t[0])) ligandToSamples.put(t[0], new HashSet<>());
				ligandToSamples.get(t[0]).add(header[i]);
			}

			if (!t[1].equals("DMSO") && !t[0].isEmpty() && !t[0].equals("Serum"))
			{
				inhibitors.add(t[1]);
			}
		}

//		inhibitors.stream().sorted().forEach(System.out::println);
//		System.exit(0);

		String outBase = "/home/ozgun/Analyses/HillPaper/LigandEffects/" + cellLineName + "/";

		for (String ligand : ligandToSamples.keySet())
		{
			String dir = outBase + ligand + "/";

			Files.createDirectories(Paths.get(dir));

			BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + "/parameters.txt"));

			writer.write(
				"proteomics-values-file = ../../../DataS1/converted/" + cellLineName + ".txt\n" +
				"id-column = ID\n" +
				"symbols-column = Symbols\n" +
				"sites-column = Sites\n" +
				"effect-column = Effect\n" +
				"\n" +
				"value-transformation = significant-change-of-mean-paired\n" +
				"\n" +
				"fdr-threshold-for-data-significance = 0.1 protein\n" +
				"fdr-threshold-for-data-significance = 0.1 phosphoprotein\n" +
				"pool-proteomics-for-fdr-adjustment = true\n" +
				"\n" +
				"color-saturation-value = 10\n" +
				"show-all-genes-with-proteomic-data = true\n" +
				"\n" +
				"calculate-network-significance = true\n" +
				"permutations-for-significance = 10000\n" +
				"fdr-threshold-for-network-significance = 0.1\n" +
				"use-network-significance-for-causal-reasoning = true\n");

			Map<String, String> pairing = getPairing(controls, ligandToSamples.get(ligand));
			for (String control : pairing.keySet())
			{
				writer.write("\ncontrol-value-column = " + control);
				writer.write("\ntest-value-column = " + pairing.get(control));
			}

			writer.close();
		}
	}

	private static void generateInhibitorAnalysisFolders(String cellLineName) throws IOException
	{
		String[] header = Files.lines(Paths.get("/home/ozgun/Analyses/HillPaper/DataS1/converted/" +
			cellLineName + ".txt")).findFirst().get().split("\t");

		Map<String, Set<String>> ligandToControls = new HashMap<>();
		Map<String, Map<String, Set<String>>> testSamplesMap = new HashMap<>();

		for (int i = 4; i < header.length; i++)
		{
			String[] t = header[i].split("-");

			if (!t[0].isEmpty() && !t[0].equals("Serum") && t[1].equals("DMSO"))
			{
				if (!ligandToControls.containsKey(t[0])) ligandToControls.put(t[0], new HashSet<>());
				ligandToControls.get(t[0]).add(header[i]);
			}

			if (!t[1].equals("DMSO") && !t[0].isEmpty() && !t[0].equals("Serum"))
			{
				if (!testSamplesMap.containsKey(t[0])) testSamplesMap.put(t[0], new HashMap<>());
				if (!testSamplesMap.get(t[0]).containsKey(t[1])) testSamplesMap.get(t[0]).put(t[1], new HashSet<>());
				testSamplesMap.get(t[0]).get(t[1]).add(header[i]);
			}
		}

		Map<String, Set<String>> drugToTargets = readDrugToTargets();

		String outBase = "/home/ozgun/Analyses/HillPaper/InhibitorEffects/" + cellLineName + "/";

		for (String ligand : ligandToControls.keySet())
		{
			for (String inhibitor : testSamplesMap.get(ligand).keySet())
			{
				String dir = outBase + ligand + "/" + inhibitor + "/";

				Files.createDirectories(Paths.get(dir));

				BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + "/parameters.txt"));

				writer.write(
					"proteomics-values-file = ../../../../DataS1/converted/" + cellLineName + ".txt\n" +
					"id-column = ID\n" +
					"symbols-column = Symbols\n" +
					"sites-column = Sites\n" +
					"effect-column = Effect\n" +
					"\n" +
					"value-transformation = significant-change-of-mean-paired\n" +
					"\n" +
					"fdr-threshold-for-data-significance = 0.1 protein\n" +
					"fdr-threshold-for-data-significance = 0.1 phosphoprotein\n" +
					"pool-proteomics-for-fdr-adjustment = true\n" +
					"\n" +
					"color-saturation-value = 10\n" +
					"show-all-genes-with-proteomic-data = true\n" +
					"\n" +
					"calculate-network-significance = true\n" +
					"permutations-for-significance = 10000\n" +
					"fdr-threshold-for-network-significance = 0.1\n" +
					"use-network-significance-for-causal-reasoning = true\n\n");

				for (String target : drugToTargets.get(inhibitor))
				{
					writer.write("gene-activity = " + target + " i\n");
				}

				Map<String, String> pairing = getPairing(ligandToControls.get(ligand), testSamplesMap.get(ligand).get(inhibitor));
				for (String control : pairing.keySet())
				{
					writer.write("\ncontrol-value-column = " + control);
					writer.write("\ntest-value-column = " + pairing.get(control));
				}

				writer.close();
			}
		}
	}

	private static Map<String, String> getPairing(Set<String> controls, Set<String> tests)
	{
		Map<String, String> map = new HashMap<>();

		for (String control : controls)
		{
			String suffix = control.substring(control.lastIndexOf("-") + 1);
			String chosen = null;

			for (String test : tests)
			{
				if (test.substring(test.lastIndexOf("-") + 1).equals(suffix))
				{
					chosen = test;
					break;
				}
			}

			if (chosen != null) map.put(control, chosen);
		}
		return map;
	}

	private static void writeAllExpectedChanges(String baseDir, String outFile) throws IOException
	{
		Map<String, Map<String, Map<String, Map<String, Map<Integer, Set<String>>>>>> map = getAllExpectedChanges(baseDir);
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		for (String cell : map.keySet())
		{
			for (String ligand : map.get(cell).keySet())
			{
				for (String drug : map.get(cell).get(ligand).keySet())
				{
					for (String target : map.get(cell).get(ligand).get(drug).keySet())
					{
						for (Integer expectation : map.get(cell).get(ligand).get(drug).get(target).keySet())
						{
							writer.write(ArrayUtil.merge("\t", new String[]{cell, ligand, drug, target,
								expectation.toString(), CollectionUtil.merge(map.get(cell).get(ligand).get(drug).get(target).get(expectation), " ")}) + "\n");
						}
					}
				}
			}
		}

		writer.close();
	}

	private static Map<String, Map<String, Map<String, Map<String, Map<Integer, Set<String>>>>>> getAllExpectedChanges(String baseDir) throws IOException
	{
		Map<String, Set<String>> drugToTargets = readDrugToTargets();

		// cell-line --> ligand --> drug --> target --> expected-change --> sources
		Map<String, Map<String, Map<String, Map<String, Map<Integer, Set<String>>>>>> map = new HashMap<>();

		for (File cellDir : new File(baseDir).listFiles())
		{
			map.put(cellDir.getName(), new HashMap<>());

			for (File ligandDir : cellDir.listFiles())
			{
				String resultFile = ligandDir.getPath() + "/results.txt";

				if (Files.exists(Paths.get(resultFile)))
				{
					map.get(cellDir.getName()).put(ligandDir.getName(), determineExpectationsAfterInhibitions(
						resultFile, drugToTargets));
				}
			}
		}

		return map;
	}

	private static Map<String, Map<String, Map<Integer, Set<String>>>> determineExpectationsAfterInhibitions(String resultFile, Map<String, Set<String>> drugToTargets) throws IOException
	{
		// drug --> target --> expected-change --> sources
		Map<String, Map<String, Map<Integer, Set<String>>>> map = new HashMap<>();

		Files.lines(Paths.get(resultFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			int relSign = t[1].startsWith("d") ? -1 : 1;

			int actChange = relSign * (int) Math.signum(Double.valueOf(t[8]));

			if (actChange == 1)
			{
				for (String drug : drugToTargets.keySet())
				{
					if (drugToTargets.get(drug).contains(t[0]))
					{
						Integer expChange = (int)(-Math.signum(Double.valueOf(t[8])));

						if (!map.containsKey(drug)) map.put(drug, new HashMap<>());
						if (!map.get(drug).containsKey(t[7])) map.get(drug).put(t[7], new HashMap<>());
						if (!map.get(drug).get(t[7]).containsKey(expChange)) map.get(drug).get(t[7]).put(expChange, new HashSet<>());
						map.get(drug).get(t[7]).get(expChange).add(t[0]);
					}
				}
			}
		});

		return map;
	}

	private static Map<String, Set<String>> readDrugToTargets() throws IOException
	{
		return Files.lines(Paths.get("/home/ozgun/Analyses/HillPaper/drug-targets.txt"))
			.map(l -> l.split("\t")).collect(Collectors.toMap(t -> t[0],
				t -> new HashSet<>(Arrays.asList(t[1].split(" ")))));
	}

	private static void testHypotheses() throws IOException
	{
		Map<Tuple, Double> pvals = new HashMap<>();

		Files.lines(Paths.get("/home/ozgun/Analyses/HillPaper/expected.txt"))
			.map(l -> l.split("\t")).forEach(hy ->
		{
			String cell = hy[0];
			String ligand = hy[1];
			String inh = hy[2];
			String trg = hy[3];
			Integer exp = Integer.valueOf(hy[4]);
			String sources = hy[5];

//			if (inh.equals("GSK690693_GSK1120212")) return;

			String dataFile = "/home/ozgun/Analyses/HillPaper/DataS1/converted/" + cell + ".txt";
			String[] header = loadHeader(dataFile);
			double[] values = loadValues(dataFile, trg);

			String ctrlPrefix = ligand + "-DMSO-";
			String testPrefix = ligand + "-" + inh + "-";

			Set<String> ctrlNames = getColNames(header, ctrlPrefix);
			Set<String> testNames = getColNames(header, testPrefix);

			Map<String, String> pairing = getPairing(ctrlNames, testNames);

			double[][] ct = getSubsets(values, header, pairing);

			Tuple tuple = TTest.testPaired(ct[0], ct[1]);
			tuple.v = getMeanChangePaired(ct[0], ct[1]);
			tuple.v /= Summary.stdev(values);
			tuple.v *= exp;

			if (-Math.log(tuple.p) > 8 && tuple.v > 0 && tuple.v < 0.2)
			{
				System.out.println();
			}

			pvals.put(tuple, tuple.p);

			System.out.println(ArrayUtil.getString("\t", cell, ligand, inh, trg, exp, sources) + "\t" + tuple.v + "\t" + (-Math.log(tuple.p)));

			// comparing effect sizes

//			String prevCtrlPrefix = "Serum-DMSO-";
//			ctrlNames = getColNames(header, prevCtrlPrefix);
//			testNames = getColNames(header, ctrlPrefix);
//
//			pairing = getPairing(ctrlNames, testNames);
//
//			ct = getSubsets(values, header, pairing);
//
//			Tuple tuple2 = TTest.testPaired(ct[0], ct[1]);
//			tuple2.v = getMeanChangePaired(ct[0], ct[1]);
//			tuple2.v /= Summary.stdev(values);
//
//			System.out.println(tuple2.v + "\t" + (-Math.log(tuple.p)));
		});

		List<Tuple> select = FDR.select(pvals, null, 0.05);
		System.out.println("select.size() = " + select.size());
		System.out.println("Agree: " + select.stream().filter(t -> t.v > 0).count());
		System.out.println("Disagree: " + select.stream().filter(t -> t.v < 0).count());

		double pValueThreshold = FDR.getPValueThreshold(pvals, null, 0.05);
		double logP = Math.log(pValueThreshold);
		System.out.println("pValueThreshold = " + pValueThreshold);
		System.out.println("-logP = " + -logP);

		Histogram his = new Histogram(1);
		his.setBorderAtZero(true);
		pvals.keySet().forEach(t -> his.count(t.v));
		his.print();
	}

	private static String[] loadHeader(String dataFile)
	{
		String l = FileUtil.lines(dataFile).findFirst().get();

		for (int i = 0; i < 4; i++)
		{
			l = l.substring(l.indexOf("\t") + 1);
		}

		return l.split("\t");
	}

	private static double[] loadValues(String dataFile, String abID)
	{
		String line = FileUtil.lines(dataFile).filter(l -> l.startsWith(abID + "\t")).findFirst().get();

		for (int i = 0; i < 4; i++)
		{
			line = line.substring(line.indexOf("\t") + 1);
		}

		String[] t = line.split("\t");
		return ArrayUtil.readToDoubleArrayPreserveMissing(t);
	}

	private static double[][] getSubsets(double[] vals, String[] header, Map<String, String> pairing)
	{
		double[] ctrl = new double[pairing.size()];
		double[] test = new double[pairing.size()];

		int j = 0;

		for (String contName : pairing.keySet())
		{
			ctrl[j]   = vals[ArrayUtil.indexOf(header, contName)];
			test[j++] = vals[ArrayUtil.indexOf(header, pairing.get(contName))];
		}

		return new double[][]{ctrl, test};
	}

	private static double getMeanChangePaired(double[] ctrl, double[] test)
	{
		double[] m = new double[test.length];
		for (int i = 0; i < m.length; i++)
		{
			m[i] = test[i] - ctrl[i];
		}
		return Summary.mean(m);
	}

	private static Set<String> getColNames(String[] header, String prefix)
	{
		Set<String> names = new HashSet<>();
		for (int i = 0; i < header.length; i++)
		{
			if (header[i].startsWith(prefix)) names.add(header[i]);
		}
		return names;
	}
}
