package org.panda.misc.analyses;

import org.panda.utility.statistics.KernelDensityPlot;
import org.panda.utility.statistics.TSNE;
import org.panda.utility.statistics.TSNEPlot;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class GordonsExpressionMuts
{
	public static void main(String[] args) throws IOException
	{
		String dir = "/home/babur/Documents/Analyses/Gordon-RPPA-pertubations/";
		String inFile = dir + "MCF10A_Set126_20170913.csv";
		String outFile = dir + "126/data.txt";
//		String inFile = dir + "MCF10A_Set134_20170913.csv";
//		String outFile = dir + "134/data.txt";

//		writeInCPFormat(inFile, outFile);
		generateSampleScatterPlots(inFile);
//		generateAddedGeneExpressionPlots(inFile);
	}

	public static void writeInCPFormat(String inFile, String outFile) throws IOException
	{
		Map<String, String> platMap = loadPlatform("/home/babur/Documents/Analyses/SMMART/Patient1/RPPA/Sample1/platform.txt");
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write("ID\tSymbols\tSites\tEffect");
		String[] abNames = Files.lines(Paths.get(inFile)).skip(6).findFirst().get().split("\t");
		for (int i = 0; i < abNames.length; i++)
		{
			if (abNames[i].length() > 4) abNames[i] = abNames[i].substring(0, abNames[i].length() - 4);
		}

		Map<String, double[]> sampleToVals = new HashMap<>();
		List<String> samples = new ArrayList<>();
		Files.lines(Paths.get(inFile)).skip(11).map(l -> l.split("\t")).forEach(t ->
		{
			String sample = t[7];
			samples.add(sample);
			double[] vals = new double[t.length - 9];
			for (int i = 0; i < vals.length; i++)
			{
				vals[i] = Double.valueOf(t[i + 9]);
			}
			sampleToVals.put(sample, vals);
		});

		for (String sample : samples)
		{
			writer.write("\t" + sample);
		}

		for (int i = 9; i < abNames.length; i++)
		{
			String info = platMap.get(abNames[i].toLowerCase());

			if (info == null)
			{
				System.err.println(abNames[i]);
			}

			writer.write("\n" + info);

			for (String sample : samples)
			{
				writer.write("\t" + sampleToVals.get(sample)[i-9]);
			}
		}

		writer.close();
	}

	public static Map<String, String> loadPlatform(String filename) throws IOException
	{
		Map<String, String> map = new HashMap<>();
		Files.lines(Paths.get(filename)).skip(1).map(l -> l.split("\t")).forEach(t ->
			map.put(t[1].toLowerCase(), t[1] + "\t" + t[2] + "\t" + (t.length > 3 ? t[3] : "") + "\t" +
				(t.length > 4 ? t[4] : "")));
		return map;
	}

	public static void generateAddedGeneExpressionPlots(String inFile) throws IOException
	{
		Map<String, Integer> geneIndices = loadTotalProteinGeneIndices(inFile);
		Map<String, Map<String, double[]>> sampleMap = loadFileRows(inFile);

		for (String gene : sampleMap.keySet())
		{
			Integer ind = geneIndices.get(gene);
			if (ind != null)
			{
				Map<String, double[]> sampleToVals = sampleMap.get(gene);
				Map<String, Set<String>> groups = identifyGroups(sampleToVals.keySet());
				Map<String, double[]> groupsToVals = new LinkedHashMap<>();
				for (String groupName : groups.keySet())
				{
					double[] v = new double[groups.get(groupName).size()];
					int i = 0;
					for (String sample : groups.get(groupName))
					{
						v[i++] = sampleToVals.get(sample)[ind];
					}
					groupsToVals.put(groupName, v);
				}

				KernelDensityPlot.plot(gene, groupsToVals);
			}
		}
	}

	public static Map<String, Integer> loadTotalProteinGeneIndices(String inFile) throws IOException
	{
		String[] abName = Files.lines(Paths.get(inFile)).skip(1).findFirst().get().split("\t");
		String[] geneName = Files.lines(Paths.get(inFile)).skip(3).findFirst().get().split("\t");

		Map<String, Integer> map = new HashMap<>();

		for (int i = 9; i < abName.length; i++)
		{
			if (!abName[i].contains("_p"))
			{
				if (map.containsKey(geneName[i]))
				{
					System.out.println("Already have this gene = " + geneName[i]);
				}
				else
				{
					map.put(geneName[i], i - 9);
				}
			}
		}
		return map;
	}

	public static void generateSampleScatterPlots(String inFile) throws IOException
	{
		Map<String, Map<String, double[]>> sampleMap = loadFileRows(inFile);

		for (String gene : sampleMap.keySet())
		{
			Map<String, double[]> map = sampleMap.get(gene);
			if (map.size() > 10)
			{
				System.out.println("\ngene = " + gene);
				plotSampleTSNE(gene, map);
			}
		}
	}

	private static boolean isSilent(String sample)
	{
		for (String s : sample.split("_"))
		{
			if (s.charAt(0) == s.charAt(s.length() - 1) && Character.isAlphabetic(s.charAt(0)) && Character.isUpperCase(s.charAt(0)))
			{
				boolean allDigits = true;
				for (int i = 1; i < s.length() - 1; i++)
				{
					if (!Character.isDigit(s.charAt(i)))
					{
						allDigits = false;
					}
				}
				if (allDigits) return true;
			}
		}
		return false;
	}

	private static boolean isControl(String sample)
	{
		return  sample.contains("mCherry") || sample.contains("Luc") || sample.contains("pBabe");
	}

	private static Map<String, Map<String, double[]>> loadFileRows(String inFile) throws IOException
	{
		Map<String, Map<String, double[]>> sampleMap = new HashMap<>();
		String controlKey = "Control";

		Files.lines(Paths.get(inFile)).skip(11).map(l -> l.split("\t")).forEach(t ->
		{
			String sample = t[7];
			if (t[7].startsWith("MP")) t[7] = t[7].substring(t[7].indexOf("_") + 1);

			String[] tt = t[7].split("_");
			String gene = tt[0].trim();

			if (isControl(sample))
			{
				gene = controlKey;
			}

			double[] vals = new double[t.length - 9];
			for (int i = 0; i < vals.length; i++)
			{
				vals[i] = Double.valueOf(t[i + 9]);
			}

			if (!sampleMap.containsKey(gene)) sampleMap.put(gene, new HashMap<>());
			sampleMap.get(gene).put(sample, vals);
		});

		Map<String, double[]> control = sampleMap.remove(controlKey);
		if (control != null)
		{
			sampleMap.keySet().stream().forEach(gene -> sampleMap.get(gene).putAll(control));
		}

		return sampleMap;
	}

	private static void plotSampleTSNE(String gene, Map<String, double[]> data)
	{
		Map<String, Set<String>> groups = identifyGroups(data.keySet());
		TSNEPlot.plot(gene, data, new ArrayList<>(groups.values()).toArray(new Set[groups.size()]));
	}
	
	private static Map<String, Set<String>> identifyGroups(Set<String> samples)
	{
		Set<String> wt = samples.stream().filter(s -> s.contains("_WT_")).collect(Collectors.toSet());
		Set<String> silent = samples.stream().filter(GordonsExpressionMuts::isSilent).collect(Collectors.toSet());
		Set<String> control = samples.stream().filter(GordonsExpressionMuts::isControl).collect(Collectors.toSet());
		Set<String> mut = samples.stream().filter(s -> !wt.contains(s) && !silent.contains(s) && !control.contains(s)).collect(Collectors.toSet());

		Map<String, Set<String>> map = new LinkedHashMap<>();
		map.put("Mutated", mut);
		map.put("Silent", silent);
		map.put("Wild-Type", wt);
		map.put("Control", control);

		return map;
	}
}
