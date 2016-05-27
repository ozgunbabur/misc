package org.panda.misc.analyses;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.misc.PathwayEnrichmentSIFGenerator;
import org.panda.resource.ChEBI;
import org.panda.resource.PCPathway;
import org.panda.resource.network.PathwayCommons;
import org.panda.utility.FileUtil;
import org.panda.utility.FormatUtil;
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
import java.text.DecimalFormat;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.toSet;

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
//		Map<String, String> metToNameMap = getMetToNameMap();

		Map<String, String> protToSymMap = getProtIDToSymbol("AD_vs_AC");
		protToSymMap.putAll(getProtIDToSymbol("UD_vs_UC"));

		PathwayCommons pc = new PathwayCommons();
		Graph graph = pc.getGraph(SIFEnum.CONSUMPTION_CONTROLLED_BY, SIFEnum.CONTROLS_PRODUCTION_OF, SIFEnum.USED_TO_PRODUCE);

		for (String comparison : comparisons)
		{
			System.out.println("sheetName = " + comparison);
			Map<String, Tuple> metComp = getMetabolicComparison(comparison);
			Map<String, Tuple> protComp = getProteomicComparison(comparison);

			generatePathsBetween(metComp, protComp, metToChEBIMap, protToSymMap, comparison, graph);
			generateCommonNeighborhood(metComp, protComp, metToChEBIMap, protToSymMap, comparison, graph);
			generateEnrichmentGraphs(metComp, protComp, metToChEBIMap, protToSymMap, comparison);
		}
	}

	private static void generatePathsBetween(Map<String, Tuple> metComp, Map<String, Tuple> protComp,
		Map<String, String> metToChEBIMap, Map<String, String> protToSymMap, String comparison, Graph graph) { try
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

		FileUtil.replaceNodeNamesInSIFFile(sifFile, ChEBI.get().getIdToNameMapping());

		// Write data colors

		ValToColor vtc = new ValToColor(new double[]{-10, 0, 10}, new Color[]{Color.BLUE, Color.WHITE, Color.RED});
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(sifFile.substring(0, sifFile.lastIndexOf(".")) + ".format"));

		writer.write("node\tall-nodes\tcolor\t255 255 255\n");
		writer.write("node\tall-nodes\tbordercolor\t0 0 0\n");
		for (String metID : metIDs)
		{
			String name = ChEBI.get().getName(metToChEBIMap.get(metID));
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
		Map<String, String> metToChEBIMap, Map<String, String> protToSymMap, String comparison, Graph graph) { try
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

		FileUtil.replaceNodeNamesInSIFFile(sifFile, ChEBI.get().getIdToNameMapping());

		// Write data colors

		ValToColor vtc = new ValToColor(new double[]{-10, 0, 10}, new Color[]{Color.BLUE, Color.WHITE, Color.RED});
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(sifFile.substring(0, sifFile.lastIndexOf(".")) + ".format"));
		writer.write("node\tall-nodes\tcolor\t255 255 255\n");
		writer.write("node\tall-nodes\tbordercolor\t0 0 0\n");
		for (String metID : metComp.keySet())
		{
			String name = ChEBI.get().getName(metToChEBIMap.get(metID));
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

	private static void generateEnrichmentGraphs(Map<String, Tuple> metComp, Map<String, Tuple> protComp,
		Map<String, String> metToChEBIMap, Map<String, String> protToSymMap, String comparison) { try
	{
		List<String> metIDs = FDR.select(FDR.extractPval(metComp), null, FDR_THR);
		List<String> protIDs = FDR.select(FDR.extractPval(protComp), null, FDR_THR);

		System.out.println("metIDs = " + metIDs.size());
		System.out.println("protIDs = " + protIDs.size());

		Set<String> mets = metIDs.stream().filter(metToChEBIMap::containsKey).map(metToChEBIMap::get)
			.collect(Collectors.toSet());

		Set<String> prots = protIDs.stream().filter(protToSymMap::containsKey).map(protToSymMap::get)
			.collect(Collectors.toSet());

		Set<String> mols = new HashSet<>(mets);
		mols.addAll(prots);

//		Map<String, Double>[] metEnrichVals = PCPathway.get().getEnrichmentPvals(mets, null, 5, 300);
//		Map<String, Double>[] protEnrichVals = PCPathway.get().getEnrichmentPvals(prots, null, 5, 300);
//		Map<String, Double>[] allEnrichVals = PCPathway.get().getEnrichmentPvals(mols, null, 5, 300);

//		if (!mets.isEmpty()) PCPathway.get().writeEnrichmentResults(mets, 3, 300, DIR + comparison + "-met-enrichment.txt");
//		if (!prots.isEmpty()) PCPathway.get().writeEnrichmentResults(prots, 3, 300, DIR + comparison + "-prot-enrichment.txt");
		if (!mols.isEmpty())
		{
			PCPathway.get().writeEnrichmentResults(mols, 3, 300, DIR + comparison + "-enrichment.txt");

			Map<String, Double>[] pvals = PCPathway.get().getEnrichmentPvals(mols, null, 3, 300);
			Map<String, Double> qvals = FDR.getQVals(pvals[0], pvals[1]);
			OptionalDouble thr = pvals[0].keySet().stream().filter(id -> qvals.get(id) < 0.1).mapToDouble(pvals[0]::get).max();
			if (!thr.isPresent()) return;

			Set<String> pathwayIDs = pvals[0].keySet().stream().filter(id -> pvals[0].get(id) <= thr.getAsDouble()).collect(toSet());

			PathwayEnrichmentSIFGenerator sg = new PathwayEnrichmentSIFGenerator();
			sg.setMolecules(mols);
			sg.setOwlFilename("/home/babur/Documents/PC/PathwayCommons.8.Detailed.BIOPAX.owl");
			sg.setBlacklistFilename("/home/babur/Documents/PC/blacklist.txt");
			sg.setPathwayIDs(pathwayIDs);
			String[][] groupTerms = {
				new String[]{"ECM"},
				new String[]{"Endosomal", "antigen"},
				new String[]{"coagulation", "fibrin", "plasminogen", "clot", "platelet"},
				new String[]{"transport"},
				new String[]{"vitamin"},
				new String[]{"integrin"},
				new String[]{"AP1"},
				new String[]{"integrin"},
				new String[]{"collagen"},
				new String[]{"Ca2+", "calcium"},
			};
			sg.setGroupTerms(groupTerms);
			sg.setTypes(SIFEnum.CONTROLS_STATE_CHANGE_OF, SIFEnum.CONTROLS_EXPRESSION_OF, SIFEnum.USED_TO_PRODUCE,
				SIFEnum.CONTROLS_PRODUCTION_OF, SIFEnum.CONSUMPTION_CONTROLLED_BY, SIFEnum.IN_COMPLEX_WITH);

			ValToColor vtc = new ValToColor(new double[]{-10, 0, 10}, new Color[]{Color.BLUE, Color.WHITE, Color.RED});

			metComp.keySet().stream().filter(metToChEBIMap::containsKey).forEach(id ->
			{
				sg.addNodeColor(metToChEBIMap.get(id), vtc.getColor(metComp.get(id).v));
				sg.addNodeTooltip(metToChEBIMap.get(id),
					String.valueOf(FormatUtil.roundToSignificantDigits(metComp.get(id).p, 4)));
			});

			protComp.keySet().stream().filter(protToSymMap::containsKey).forEach(id ->
			{
				sg.addNodeColor(protToSymMap.get(id), vtc.getColor(protComp.get(id).v));
				sg.addNodeTooltip(protToSymMap.get(id),
					String.valueOf(FormatUtil.roundToSignificantDigits(protComp.get(id).p, 4)));
			});

			String fileNameWithoutExtension = DIR + comparison + "-enriched";
			sg.write(fileNameWithoutExtension);
		}
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
