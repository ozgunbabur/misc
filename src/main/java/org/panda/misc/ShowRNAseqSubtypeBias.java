package org.panda.misc;

import org.panda.resource.tcga.ExpressionReader;
import org.panda.utility.Tuple;
import org.panda.utility.statistics.Correlation;
import org.panda.utility.statistics.KernelDensityEstimation;
import org.panda.utility.statistics.KernelDensityPlot;
import org.panda.utility.statistics.TTest;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class ShowRNAseqSubtypeBias
{
	static final String STUDY = "BRCA";
	static final String TCGA_FILE = "/home/ozgun/Data/TCGA/" + STUDY + "/expression.txt";
	static final String SUBTYPES_FILE = "/home/ozgun/Data/TCGA/pancan_samples.txt";
	static final Set<String> CONSIDER_ONLY = new HashSet<>(Arrays.asList("LumA", "LumB", "Her2", "Basal"));
	static String[] genes = new String[]{"ESR1", "KLF4", "EGFR", "CDKN1A", "ERBB2", "ERBB3"};

	public static void main(String[] args) throws IOException
	{
		plot();
	}

	public static void plot() throws IOException
	{
		HashSet<String> geneSet = new HashSet<>(Arrays.asList(genes));
		ExpressionReader er = new ExpressionReader(TCGA_FILE, geneSet, 15);

		Map<String, Set<String>> subtypeMap = loadSubtypes();

		if (!CONSIDER_ONLY.isEmpty())
			new HashSet(subtypeMap.keySet()).stream().filter(type -> !CONSIDER_ONLY.contains(type)).forEach(type ->
				subtypeMap.remove(type));


		Set<String> allSamples = er.getSamples();
		subtypeMap.values().forEach(set -> set.retainAll(allSamples));


		Map<String, String[]> subArrays = subtypeMap.keySet().stream().collect(Collectors.toMap(
			type -> type, type -> subtypeMap.get(type).toArray(new String[0])));

		Map<String, Set<String>> negativeMap = new HashMap<>();
		subtypeMap.keySet().forEach(type -> negativeMap.put(type, new HashSet<>()));
		for (String type1 : subtypeMap.keySet())
		{
			for (String type2 : negativeMap.keySet())
			{
				if (!type1.equals(type2)) negativeMap.get(type2).addAll(subtypeMap.get(type1));
			}
		}

		Map<String, String[]> negArrays = negativeMap.keySet().stream().collect(Collectors.toMap(
			type -> type, type -> negativeMap.get(type).toArray(new String[0])));

		Map<String, Map<String, Tuple>> tMap = new HashMap<>();

		Map<String, Map<String, double[]>> valsMap = new HashMap<>();

		for (String gene : geneSet)
		{
			tMap.put(gene, new HashMap<>());
			valsMap.put(gene, new HashMap<>());

			for (String type : subArrays.keySet())
			{
				double[] posVals = er.getGeneAlterationArray(gene, subArrays.get(type));
				double[] negVals = er.getGeneAlterationArray(gene, negArrays.get(type));

				Tuple result = TTest.test(negVals, posVals);
				tMap.get(gene).put(type, result);

				valsMap.get(gene).put(type, posVals);
			}
		}

		tMap.forEach((gene, map) ->
		{
			System.out.println("\ngene = " + gene);
			map.forEach((type, tuple) -> System.out.println(type + "\t" + tuple.v + "\t" + tuple.p));
		});

		for (String gene : valsMap.keySet())
		{
			Map<String, double[]> map = valsMap.get(gene);
			Map<KernelDensityEstimation, String> kdeMap = new LinkedHashMap<>();

			map.keySet().stream().sorted().forEach(type -> kdeMap.put(new KernelDensityEstimation(map.get(type)), type));
			KernelDensityPlot.plotMap(gene, kdeMap);
		}
	}

	public static Map<String, Set<String>> loadSubtypes() throws IOException
	{
		Map<String, Set<String>> map = new HashMap<>();
		Files.lines(Paths.get(SUBTYPES_FILE)).skip(1).map(l -> l.split("\t")).filter(t -> t[2].equals(STUDY))
			.forEach(t ->
			{
				if (!map.containsKey(t[3])) map.put(t[3], new HashSet<>());
				map.get(t[3]).add(t[1]);
			});
		return map;
	}
}
