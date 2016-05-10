package org.panda.misc.analyses;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.ChEBI;
import org.panda.resource.network.PathwayCommons;
import org.panda.utility.FileUtil;
import org.panda.utility.Tuple;
import org.panda.utility.ValToColor;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.query.QueryExecuter;
import org.panda.utility.graph.query.QueryGraphObject;
import org.panda.utility.statistics.FDR;

import java.awt.*;
import java.io.BufferedWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class Sheep4x4
{
	static final String DIR = "/home/babur/Documents/Analyses/Sheep4x4/";
	static final String METABOLOMICS_FILE = DIR + "Metabolomics_4x4_UrineAndAF_DESeq2results_14Apr2016_withChEBIv2_withNormalized.xlsx";
	static final String PROTEOMICS_FILE = DIR + "Proteomics_4x4_UrineAndAF_DESeq2results_19Apr2016.xlsx";

	static final double FDR_THR = 0.1;

	public static void main(String[] args)
	{
		List<String> comparisons = getComparisonSheetNames(METABOLOMICS_FILE);
		comparisons.retainAll(getComparisonSheetNames(PROTEOMICS_FILE));

		Map<String, String> metToChEBIMap = getMetToChEBIMap();
		Map<String, String> metToNameMap = getMetToNameMap();

		Map<String, String> protToSymMap = getProtIDToSymbol("AD_vs_AC");
		protToSymMap.putAll(getProtIDToSymbol("UD_vs_UC"));

		PathwayCommons pc = new PathwayCommons();
		Graph graph = pc.getGraph(SIFEnum.CONSUMPTION_CONTROLLED_BY, SIFEnum.CONTROLS_PRODUCTION_OF, SIFEnum.USED_TO_PRODUCE);

		for (String comparison : comparisons)
		{
			System.out.println("sheetName = " + comparison);
			Map<String, Tuple> metComp = getMetabolicComparison(comparison);
			Map<String, Tuple> protComp = getProteomicComparison(comparison);

			generatePathsBetween(metComp, protComp, metToChEBIMap, metToNameMap, protToSymMap, comparison, graph);
			generateCommonNeighborhood(metComp, protComp, metToChEBIMap, metToNameMap, protToSymMap, comparison, graph);
		}
	}

	private static void generatePathsBetween(Map<String, Tuple> metComp, Map<String, Tuple> protComp,
		Map<String, String> metToChEBIMap, Map<String, String> metToNameMap, Map<String, String> protToSymMap,
		String comparison, Graph graph) { try
	{
		List<String> metIDs = FDR.select(FDR.extractPval(metComp), null, FDR_THR);
		List<String> protIDs = FDR.select(FDR.extractPval(protComp), null, FDR_THR);

		System.out.println("metIDs = " + metIDs.size());
		System.out.println("protIDs = " + protIDs.size());

		List<String> mets = metIDs.stream().filter(metToChEBIMap::containsKey).map(metToChEBIMap::get)
			.collect(Collectors.toList());

		List<String> prots = protIDs.stream().filter(protToSymMap::containsKey).map(protToSymMap::get)
			.collect(Collectors.toList());

		Set<String> seed = new HashSet<>(mets);
		seed.addAll(prots);

		if (seed.isEmpty() || seed.size() < 2) return;

		Set<QueryGraphObject> result = QueryExecuter.pathsBetween(seed, graph, 1, true, 0, true);
		if (result.isEmpty()) return;

		String sifFile = DIR + comparison + "-paths-between.sif";
		graph.write(sifFile, result);

		// Replace ChEBI Ids with names. Prioritize the names in the data file.

		Map<String, String> mapInFile = metIDs.stream().filter(metToChEBIMap::containsKey)
			.collect(Collectors.toMap(metToChEBIMap::get, metToNameMap:: get));
		Map<String, String> chebiToNameMap = ChEBI.get().getIdToNameMapping();
		chebiToNameMap.putAll(mapInFile);

		FileUtil.replaceNodeNamesInSIFFile(sifFile, chebiToNameMap);

		// Write data colors

		ValToColor vtc = new ValToColor(new double[]{-10, 0, 10}, new Color[]{Color.BLUE, Color.WHITE, Color.RED});
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(sifFile.substring(0, sifFile.lastIndexOf(".")) + ".format"));

		writer.write("node\tall-nodes\tcolor\t255 255 255\n");
		writer.write("node\tall-nodes\tbordercolor\t0 0 0\n");
		for (String metID : metIDs)
		{
			String name = metToNameMap.get(metID);
			writer.write("node\t" + name + "\tcolor\t" + vtc.getColorInString(metComp.get(metID).v) + "\n");
			writer.write("node\t" + name + "\ttooltip\t" + metComp.get(metID).p + "\n");
		}
		for (String protID : protIDs)
		{
			String name = protToSymMap.get(protID);
			writer.write("node\t" + name + "\tcolor\t" + vtc.getColorInString(protComp.get(protID).v) + "\n");
			writer.write("node\t" + name + "\tborderwidth\t2\n");
			writer.write("node\t" + name + "\ttooltip\t" + protComp.get(protID).p + "\n");
		}
		writer.close();
	}
	catch (Exception e){throw new RuntimeException(e);}}


	private static void generateCommonNeighborhood(Map<String, Tuple> metComp, Map<String, Tuple> protComp,
		Map<String, String> metToChEBIMap, Map<String, String> metToNameMap, Map<String, String> protToSymMap,
		String comparison, Graph graph) { try
	{
		List<String> metIDs = FDR.select(FDR.extractPval(metComp), null, FDR_THR);
		List<String> protIDs = FDR.select(FDR.extractPval(protComp), null, FDR_THR);

		System.out.println("metIDs = " + metIDs.size());
		System.out.println("protIDs = " + protIDs.size());

		List<String> mets = metIDs.stream().filter(metToChEBIMap::containsKey).map(metToChEBIMap::get)
			.collect(Collectors.toList());

		List<String> prots = protIDs.stream().filter(protToSymMap::containsKey).map(protToSymMap::get)
			.collect(Collectors.toList());

		Set<String> seed = new HashSet<>(mets);
		seed.addAll(prots);

		if (seed.isEmpty() || seed.size() < 2) return;

		Set<QueryGraphObject> result = QueryExecuter.commonNeighborhood(seed, graph, 2);
		if (result.isEmpty()) return;

		String sifFile = DIR + comparison + "-common-neighborhood.sif";
		graph.write(sifFile, result);
		Set<String> protsInSIF = FileUtil.getNodeNamesInSIFFile(sifFile).stream()
			.filter(name -> !name.startsWith("CHEBI:")).collect(Collectors.toSet());

		// Replace ChEBI Ids with names. Prioritize the names in the data file.

		Map<String, String> mapInFile = metIDs.stream().filter(metToChEBIMap::containsKey)
			.collect(Collectors.toMap(metToChEBIMap::get, metToNameMap:: get));
		Map<String, String> chebiToNameMap = ChEBI.get().getIdToNameMapping();
		chebiToNameMap.putAll(mapInFile);

		FileUtil.replaceNodeNamesInSIFFile(sifFile, chebiToNameMap);

		// Write data colors

		ValToColor vtc = new ValToColor(new double[]{-10, 0, 10}, new Color[]{Color.BLUE, Color.WHITE, Color.RED});
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(sifFile.substring(0, sifFile.lastIndexOf(".")) + ".format"));
		writer.write("node\tall-nodes\tcolor\t255 255 255\n");
		writer.write("node\tall-nodes\tbordercolor\t0 0 0\n");
		for (String metID : metComp.keySet())
		{
			String name = metToNameMap.get(metID);
			writer.write("node\t" + name + "\tcolor\t" + vtc.getColorInString(metComp.get(metID).v) + "\n");
			writer.write("node\t" + name + "\ttooltip\t" + metComp.get(metID).p + "\n");
		}
		for (String protID : protComp.keySet())
		{
			String name = protToSymMap.get(protID);
			writer.write("node\t" + name + "\tcolor\t" + vtc.getColorInString(protComp.get(protID).v) + "\n");
			writer.write("node\t" + name + "\ttooltip\t" + protComp.get(protID).p + "\n");
		}
		for (String prot : protsInSIF)
		{
			writer.write("node\t" + prot + "\tborderwidth\t2\n");
		}
		writer.close();
	}
	catch (Exception e){throw new RuntimeException(e);}}



	private static Map<String, String> getMetToChEBIMap()
	{
		return FileUtil.getXLSXAsStream(METABOLOMICS_FILE, "RawWithChEBI")
			.filter(token -> token[0].startsWith("BBid_") && !token[8].isEmpty())
			.collect(Collectors.toMap(token -> token[0], token -> token[8]));
	}

	private static Map<String, String> getChEBIToNameMap()
	{
		return FileUtil.getXLSXAsStream(METABOLOMICS_FILE, "RawWithChEBI")
			.filter(token -> token[0].startsWith("BBid_") && !token[1].isEmpty() && !token[8].isEmpty())
			.collect(Collectors.toMap(token -> token[8], token -> token[1]));
	}

	private static Map<String, String> getMetToNameMap()
	{
		return FileUtil.getXLSXAsStream(METABOLOMICS_FILE, "RawWithChEBI")
			.filter(token -> token[0].startsWith("BBid_") && !token[1].isEmpty())
			.collect(Collectors.toMap(token -> token[0], token -> token[1]));
	}

	private static Map<String, Tuple> getMetabolicComparison(String sheetName)
	{
		return FileUtil.getXLSXAsStream(METABOLOMICS_FILE, sheetName)
			.filter(token -> token[0].startsWith("BBid_") && !token[4].equals("NA") && !token[5].equals("NA"))
			.collect(Collectors.toMap(token -> token[0],
				token -> new Tuple(Double.parseDouble(token[4]), Double.parseDouble(token[5]))));
	}

	private static Map<String, Tuple> getProteomicComparison(String sheetName)
	{
		return FileUtil.getXLSXAsStream(PROTEOMICS_FILE, sheetName)
			.filter(token -> token[0].startsWith("ENS") && !token[9].isEmpty() && !token[9].equals("0.0") &&
				!token[4].equals("NA") && !token[5].equals("NA"))
			.collect(Collectors.toMap(token -> token[0],
				token -> new Tuple(Double.parseDouble(token[4]), Double.parseDouble(token[5]))));
	}

	private static List<String> getComparisonSheetNames(String filename)
	{
		return FileUtil.getXLSXSheetNames(filename).stream().filter(name -> name.contains("_vs_"))
			.collect(Collectors.toList());
	}

	private static Map<String, String> getProtIDToSymbol(String sheetName)
	{
		return FileUtil.getXLSXAsStream(PROTEOMICS_FILE, sheetName).filter(token -> token[0].startsWith("ENS")
			&& !token[9].isEmpty() && !token[9].equals("0.0"))
			.collect(Collectors.toMap(token -> token[0], token -> token[9]));
	}
}
