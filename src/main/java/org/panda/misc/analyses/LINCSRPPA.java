package org.panda.misc.analyses;

import org.panda.resource.proteomics.MDACCFormatRPPALoader;
import org.panda.resource.proteomics.RPPAIDMapper;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.SIFFileUtil;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.Summary;
import org.panda.utility.statistics.TTest;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

public class LINCSRPPA
{
	public static void main(String[] args) throws IOException
	{
//		convertDataTpCPFormat();
//		checkNormalization();
//		printExpressionRatios();
		printCommonGenesInTwoEGFStudies();
	}

	public static void convertDataTpCPFormat() throws IOException
	{
		Map<String, Map<String, Double>> mapA = addSampleSuffix(removeSamplePrefix(MDACCFormatRPPALoader.load("/home/ozgun/Data/LINCS/repA.csv")), "_A");
		Map<String, Map<String, Double>> mapB = addSampleSuffix(removeSamplePrefix(MDACCFormatRPPALoader.load("/home/ozgun/Data/LINCS/repB.csv")), "_B");
		Map<String, Map<String, Double>> mapC = addSampleSuffix(removeSamplePrefix(MDACCFormatRPPALoader.load("/home/ozgun/Data/LINCS/repC.csv")), "_C");

		Map<String, Map<String, Double>> map = mergeReplicates(mapA, mapB, mapC);
		RPPAIDMapper.writeAsCPFile(map, "/home/ozgun/Data/LINCS/data.txt", null);
	}

	public static void checkNormalization() throws IOException
	{
		Map<String, Map<String, Double>> mapA = removeSamplePrefix(MDACCFormatRPPALoader.load("/home/ozgun/Data/LINCS/repA.csv"));
		Map<String, Map<String, Double>> mapB = removeSamplePrefix(MDACCFormatRPPALoader.load("/home/ozgun/Data/LINCS/repB.csv"));
		Map<String, Map<String, Double>> mapC = removeSamplePrefix(MDACCFormatRPPALoader.load("/home/ozgun/Data/LINCS/repC.csv"));

		Histogram h = new Histogram(0.05);
		for (String ab : mapA.keySet())
		{
			List<Double> diffs1 = new ArrayList<>();
			List<Double> diffs2 = new ArrayList<>();

			for (String sample : mapA.get(ab).keySet())
			{
				diffs1.add((mapB.get(ab).get(sample) - mapA.get(ab).get(sample)) * (mapA.get(ab).get(sample) < 0 ? -1 : 1));
				diffs2.add((mapC.get(ab).get(sample) - mapA.get(ab).get(sample)) * (mapA.get(ab).get(sample) < 0 ? -1 : 1));
			}

			System.out.println("ab = " + ab);
			double[] d1 = ArrayUtil.convertToBasicDoubleArray(diffs1);
			double[] d2 = ArrayUtil.convertToBasicDoubleArray(diffs2);

			Tuple t1 = TTest.test(d1, 0);
			Tuple t2 = TTest.test(d2, 0);

			h.count(t1.p);
			h.count(t2.p);

			System.out.println("t1 = " + t1 + "\t" + Summary.stdev(ArrayUtil.convertToBasicDoubleArray(new ArrayList<>(mapA.get(ab).values()))) + "\t" + Summary.stdev(ArrayUtil.convertToBasicDoubleArray(new ArrayList<>(mapB.get(ab).values()))));
			System.out.println("t2 = " + t2 + "\t" + Summary.stdev(ArrayUtil.convertToBasicDoubleArray(new ArrayList<>(mapC.get(ab).values()))));
		}

		h.print();
	}

	static Map<String, Map<String, Double>> removeSamplePrefix(Map<String, Map<String, Double>> orig)
	{
		Map<String, Map<String, Double>> map = new HashMap<>();

		for (String ab : orig.keySet())
		{
			map.put(ab, new HashMap<>());

			for (String s : orig.get(ab).keySet())
			{
				String sample = s.substring(s.indexOf("_", s.indexOf("_") + 1) + 1);
				map.get(ab).put(sample, orig.get(ab).get(s));
			}
		}
		return map;
	}
	static Map<String, Map<String, Double>> addSampleSuffix(Map<String, Map<String, Double>> orig, String suffix)
	{
		Map<String, Map<String, Double>> map = new HashMap<>();

		for (String ab : orig.keySet())
		{
			map.put(ab, new HashMap<>());

			for (String s : orig.get(ab).keySet())
			{
				String sample = s + suffix;
				map.get(ab).put(sample, orig.get(ab).get(s));
			}
		}
		return map;
	}

	static Map<String, Map<String, Double>> mergeReplicates(Map<String, Map<String, Double>>... reps)
	{
		Map<String, Map<String, Double>> map = new HashMap<>();

		for (Map<String, Map<String, Double>> rep : reps)
		{
			for (String ab : rep.keySet())
			{
				if (!map.containsKey(ab)) map.put(ab, new HashMap<>());
				map.get(ab).putAll(rep.get(ab));
			}
		}
		return map;
	}

	static void printExpressionRatios() throws IOException
	{
		String dir = "/home/ozgun/Data/LINCS/";
		String[] times = new String[]{"1", "4", "8", "24", "48"};
		String[] ligands = new String[]{"EGF", "BMP2", "HGF", "IFNG", "OSM", "TGFB", "PBS"};

		double[][] rats = new double[times.length][ligands.length];

		for (int i = 0; i < ligands.length; i++)
		{
			for (int j = 0; j < times.length; j++)
			{
				String file = dir + File.separator + ligands[i] + File.separator + times[j] + File.separator + "causative.sif";
				long pCnt = Files.lines(Paths.get(file)).filter(l -> l.contains("phospho")).count();
				long eCnt = Files.lines(Paths.get(file)).filter(l -> l.contains("expression")).count();

//				double ratio = (pCnt == 0 && eCnt == 0) ? Double.NaN : eCnt / (double) (eCnt + pCnt);
				double ratio = eCnt + pCnt;

				rats[j][i] = ratio;
			}
		}

		for (String ligand : ligands)
		{
			System.out.print("\t" + ligand);
		}
		for (int i = 0; i < times.length; i++)
		{
			System.out.print("\n" + times[i]);

			for (int j = 0; j < ligands.length; j++)
			{
				System.out.print("\t" + (Double.isNaN(rats[i][j]) ? "" : rats[i][j]));
			}
		}
	}

	static void printCommonGenesInTwoEGFStudies() throws IOException
	{
		Set<String> set1 = SIFFileUtil.getGenesInSIFFile(
			"/home/ozgun/Analyses/CausalPath-paper/EGF-stimulation/run-fdr-0.1-relax/AsSeries/causative.sif");

		Set<String> set2 = SIFFileUtil.getGenesInSIFFile("/home/ozgun/Data/LINCS/EGF/AsSeries/causative.sif");

		CollectionUtil.printVennSets(set1, set2);
	}
}
