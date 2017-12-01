package org.panda.misc.analyses;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.KernelDensityPlot;
import org.panda.utility.statistics.Summary;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class CPTACBreast
{
	public static final String DIR = "/home/babur/Documents/Analyses/CPTACBreastCancer/";
	public static final String PROT_FILE = DIR + "Proteome/TCGA_Breast_BI_Proteome.itraq.tsv";
	public static final String PHOS_FILE = DIR + "Phosphoproteome/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq.tsv";
	public static final String SUBTYPE_FILE = DIR + "Metadata/CPTAC_TCGA_BreastCancer_select_clinical_data_r1.csv";
	public static final String CP_FILE = DIR + "CPTAC-TCGA-BRCA-data.txt";

	public static void main(String[] args) throws IOException
	{
//		convertToCpFormat();
//		printGroups();
		printForSpecificHighLowComparison();
	}

	static void convertToCpFormat() throws IOException
	{
		String[] protH = readHeader(PROT_FILE);
		String[] phosH = readHeader(PHOS_FILE);

		Map<String, Integer> protInd = getIndicesMap(protH);
		Map<String, Integer> phosInd = getIndicesMap(phosH);

		List<String> samples = protInd.keySet().stream().sorted().collect(Collectors.toList());

		Map<String, double[]> protMap = readProteomics();
		Map<String, double[]> phosMap = readPhosphoroteomics();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(CP_FILE));

		writer.write("ID\tSymbols\tSites\tEffect");
		samples.forEach(s -> FileUtil.tab_write(s, writer));

		phosMap.keySet().stream().sorted().forEach(k ->
		{
			String[] t = k.split(" ");
			List<String> sites = splitSites(t[1]);
			String id = t[0] + "-" + CollectionUtil.merge(sites, "-");

			FileUtil.lnwrite(id + "\t" + t[0] + "\t" + CollectionUtil.merge(sites, "|") + "\t", writer);

			for (String sample : samples)
			{
				FileUtil.tab_write(phosMap.get(k)[phosInd.get(sample)], writer);
			}
		});

		protMap.keySet().stream().sorted().forEach(k ->
		{
			FileUtil.lnwrite(k + "\t" + k + "\t\t", writer);

			for (String sample : samples)
			{
				FileUtil.tab_write(protMap.get(k)[protInd.get(sample)], writer);
			}
		});

		writer.close();
	}

	static Map<String, Integer> getIndicesMap(String[] samples)
	{
		Map<String, Integer> map = new HashMap<>();
		for (int i = 0; i < samples.length; i++)
		{
			if (samples[i] != null)
			{
				map.put(samples[i], i);
			}
		}
		return map;
	}

	static Map<String, double[]> readProteomics() throws IOException
	{
		Map<String, double[]> map = new HashMap<>();
		Files.lines(Paths.get(PROT_FILE)).skip(4).map(l -> l.split("\t")).forEach(t ->
		{
			double[] v = new double[t.length];

			for (int i = 0; i < t.length; i++)
			{
				try
				{
					v[i] = Double.valueOf(t[i]);
				}
				catch (Exception e)
				{
					v[i] = Double.NaN;
				}
			}

			map.put(t[0], v);
		});
		return map;
	}

	static Map<String, double[]> readPhosphoroteomics() throws IOException
	{
		Map<String, double[]> map = new HashMap<>();
		Files.lines(Paths.get(PHOS_FILE)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			double[] v = new double[t.length];

			v[0] = Double.NaN;
			for (int i = 0; i < t.length; i++)
			{
				try
				{
					v[i] = Double.valueOf(t[i]);
				}
				catch (Exception e)
				{
					v[i] = Double.NaN;
				}
			}

			map.put(t[113] + " " + t[0].split(":")[1], v);
		});
		return map;
	}

	static List<String> splitSites(String s)
	{
		List<String> sites = new ArrayList<>();

		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < s.length(); i++)
		{
			if (Character.isAlphabetic(s.charAt(i)) && sb.length() > 0)
			{
				sites.add(sb.toString().toUpperCase());
				sb = new StringBuilder();
			}

			sb.append(s.charAt(i));
		}
		sites.add(sb.toString().toUpperCase());

		return sites;
	}

	static String[] readHeader(String file) throws IOException
	{
		String[] header = Files.lines(Paths.get(file)).findFirst().get().split("\t");

		for (int i = 0; i < header.length; i++)
		{
			header[i] = normalizeSampleName(header[i]);
		}
		return header;
	}

	static String normalizeSampleName(String s)
	{
		if (Character.isUpperCase(s.charAt(0)) && !s.contains("Unshared") && s.contains("-"))
		{
			String[] t = s.split(" ")[0].split("-");
			return "TCGA-" + t[0] + "-" + t[1];
		}
		return null;
	}

	static void printGroups() throws IOException
	{
		Map<String, String> subMap = Files.lines(Paths.get(SUBTYPE_FILE)).skip(6).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> t[20].replaceAll(" ", "")));

		List<String> samples = getIndicesMap(readHeader(PROT_FILE)).keySet().stream().sorted()
			.collect(Collectors.toList());

		subMap.values().stream().distinct().forEach(type ->
		{
			System.out.println("\ntype = " + type);
			samples.forEach(s ->
				System.out.println((subMap.get(s).equals(type) ? "test" : "control") + "-value-column = " + s));
		});

		System.out.println("\nCorrelation");
		samples.forEach(s -> System.out.println("value-column = " + s));

		printMarkerGroups("ER", 3, samples);
		printMarkerGroups("PR", 4, samples);
		printMarkerGroups("HER2", 5, samples);
	}

	static void printMarkerGroups(String name, int index, List<String> samples) throws IOException
	{
		Map<String, Boolean> map = Files.lines(Paths.get(SUBTYPE_FILE)).skip(6).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> t[index].equals("Positive")));

		System.out.println("\ngene = " + name);
		samples.stream().filter(map::containsKey).forEach(s ->
			System.out.println((map.get(s) ? "test" : "control") + "-value-column = " + s));
	}

	static void printForSpecificHighLowComparison() throws IOException
	{
		String id = "SHC1-S139";
		String[] header = Files.lines(Paths.get(CP_FILE)).findFirst().get().split("\t");
		String[] row = Files.lines(Paths.get(CP_FILE)).filter(l -> l.startsWith(id + "\t")).findFirst().get().split("\t");
		double[] d = ArrayUtil.readToDoubleArrayPreserveMissing(row);
		double[] v = ArrayUtil.trimNaNs(d);
		Arrays.sort(v);

		double lowBorder = v[(int) Math.round(v.length / 3D)];
		double highBorder = v[(int) Math.round(2 * (v.length / 3D))];

		System.out.println("v[0] = " + v[0]);
		System.out.println("lowBorder = " + lowBorder);
		System.out.println("highBorder = " + highBorder);
		System.out.println("v[" + (v.length - 1) + "] = " + v[v.length - 1]);

		System.out.println("v.length = " + v.length);
//		KernelDensityPlot.plot("", v);

		for (int i = 0; i < header.length; i++)
		{
			if (d[i] <= lowBorder)
			{
				System.out.println("control-value-column = " + header[i]);
			}
			else if (d[i] >= highBorder)
			{
				System.out.println("test-value-column = " + header[i]);
			}
		}
	}
}
