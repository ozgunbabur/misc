package org.panda.misc.analyses;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.CancerGeneCensus;
import org.panda.resource.OncoKB;
import org.panda.resource.network.IPTMNet;
import org.panda.resource.network.PathwayCommons;
import org.panda.resource.network.PhosphoNetworks;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;
import org.panda.utility.graph.*;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class SMMARTPatient115
{
	public static final String DIR = "/home/babur/Documents/Analyses/SMMART/Patient115/";

	public static final String ACTIVATING_BORDER = "0 180 20";
	public static final String INHIBITING_BORDER = "180 0 20";
	public static final String NEUTRAL_BORDER = "0 0 0";
	public static final String SITE_FEATURE = "rppasite";
	public static final Color UP_COLOR = new Color(255, 100, 100);
	public static final Color DOWN_COLOR = new Color(100, 100, 255);
	public static final String CANCER_GENE_COLOR = "255 255 200";

	public static void main(String[] args) throws IOException
	{
		Map<String, Integer> mutMap = loadMutations();
		System.out.println("mutMap.size() = " + mutMap.size());

		Map<String, Double> cnaMap = loadCNA();
		System.out.println("cnaMap.size() = " + cnaMap.size());

		Map<String, Double> rnaMap = loadRNA();
		System.out.println("rnaMap.size() = " + rnaMap.size());

		Map<String, Integer> tfaMap = loadTFActivity();
		System.out.println("tfaMap.size() = " + tfaMap.size());

		Set<String> cancerGenes = loadCancerGenes();

		writeFormatFile(DIR + "network.format", mutMap, cnaMap, rnaMap, tfaMap, cancerGenes);
		writeRelationsBetweenCancerGenes(DIR + "network-between-cancer-genes.sif", mutMap, cnaMap, rnaMap, tfaMap, cancerGenes);
		writeTFActivityChangesOneByOne(DIR + "network-around-TF-", mutMap, cnaMap, rnaMap, tfaMap, cancerGenes);
		writeAlterationsAroundCancerGenesOnBRCANetwork(DIR + "network-related-to-BRCA.sif", mutMap, cnaMap, rnaMap, tfaMap, cancerGenes);
	}

	private static Map<String, Integer> loadMutations() throws IOException
	{
		Map<String, Integer> map = Files.lines(Paths.get(DIR + "Galaxy156-[OncotatorMAF].maf"))
			.filter(l -> !l.startsWith("#")).skip(1).map(l -> l.split("\t"))
			.filter(t -> t.length > 8 && (t[8].equals("Missense_Mutation") || t[8].equals("Nonsense_Mutation") || t[8].equals("Splice_Site") || t[8].startsWith("In_Frame") || t[8].startsWith("Frame_Shift")))
			.collect(Collectors.toMap(t -> t[0], t -> t[8].equals("Missense_Mutation") || t[8].startsWith("In_Frame") ? 0 : -1, Math::min));

		Files.lines(Paths.get(DIR + "SMMART16113-115_mutation_list.tsv")).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			if (!map.containsKey(t[8]) || map.get(t[8]) == 0)
			{
				map.put(t[8], t[7].contains("*") ? -1 : 0);
			}
		});

		return map;
	}

	private static Map<String, Double> loadCNA() throws IOException
	{
		return Files.lines(Paths.get(DIR + "SMMART16113-115_cna.tsv"))
			.skip(1).map(l -> l.split("\t")).collect(Collectors.toMap(t -> t[4], t ->
			{
				double v = Double.valueOf(t[7]);
				if (v > 1) return v;
				else return -(1 / v);
			}));
	}

	private static Map<String, Double> loadRNA() throws IOException
	{
		return Files.lines(Paths.get(DIR + "RNAseq-zscores.txt"))
			.skip(1).map(l -> l.split("\t")).filter(t -> Math.abs(Double.valueOf(t[1])) > 2)
			.collect(Collectors.toMap(t -> t[0], t -> Double.valueOf(t[1])));
	}

	private static Map<String, Integer> loadTFActivity() throws IOException
	{
		return Files.lines(Paths.get(DIR + "tf-enrich-signed-tcga-consensus-thr2/TF-activity-results.txt"))
			.skip(1).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> t[1].startsWith("a") ? 1 : -1));
	}

	private static Set<String> loadCancerGenes()
	{
		Set<String> set = new HashSet<>(OncoKB.get().getAllSymbols());
		set.addAll(CancerGeneCensus.get().getAllSymbols());
		return set;
	}

	private static DirectedGraph loadStateChangeNetwork()
	{
		DirectedGraph graph = new DirectedGraph("State change", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		graph.merge((DirectedGraph) PathwayCommons.get().getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF));
		graph.merge(PhosphoNetworks.get().getGraph());
		graph.merge(IPTMNet.get().getGraph(SignedType.PHOSPHORYLATES));
		return graph;
	}

	private static UndirectedGraph loadPPI()
	{
		UndirectedGraph graph = new UndirectedGraph("PPI", SIFEnum.IN_COMPLEX_WITH.getTag());
		graph.merge((UndirectedGraph) PathwayCommons.get().getGraph(SIFEnum.IN_COMPLEX_WITH));
		graph.merge((UndirectedGraph) PathwayCommons.get().getGraph(SIFEnum.INTERACTS_WITH));
		return graph;
	}

	private static GraphList loadSignedExpressionTCGAConsensus()
	{
		DirectedGraph posG = new DirectedGraph("TCGA consensus upregulates expression", SignedType.UPREGULATES_EXPRESSION.getTag());
		DirectedGraph negG = new DirectedGraph("TCGA consensus downregulates expression", SignedType.DOWNREGULATES_EXPRESSION.getTag());

		String filename = "/home/babur/Documents/PC/SignedByTCGAConsensusFiltered.sif";
		posG.load(filename, Collections.singleton(SignedType.UPREGULATES_EXPRESSION.getTag()));
		negG.load(filename, Collections.singleton(SignedType.DOWNREGULATES_EXPRESSION.getTag()));

		GraphList gl = new GraphList("Expression TCGA consensus");
		gl.addGraph(posG);
		gl.addGraph(negG);
		return gl;
	}

	private static GraphList loadTCGABRCACPNetwork()
	{
		SiteSpecificGraph posP = new SiteSpecificGraph("TCGA BRCA CP phosphorylation", SignedType.PHOSPHORYLATES.getTag());
		SiteSpecificGraph negP = new SiteSpecificGraph("TCGA BRCA CP dephosphorylation", SignedType.DEPHOSPHORYLATES.getTag());
		DirectedGraph posE = new DirectedGraph("TCGA BRCA CP upregulates expression", SignedType.UPREGULATES_EXPRESSION.getTag());
		DirectedGraph negE = new DirectedGraph("TCGA BRCA CP downregulates expression", SignedType.DOWNREGULATES_EXPRESSION.getTag());

		String filename = "/home/babur/Documents/Analyses/CPTACBreastCancer/correlation-based-phospho-0.1/causative.sif";
		posP.load(filename, Collections.singleton(SignedType.PHOSPHORYLATES.getTag()));
		negP.load(filename, Collections.singleton(SignedType.DEPHOSPHORYLATES.getTag()));

		filename = "/home/babur/Documents/Analyses/CPTACBreastCancer/correlation-based-expression-using-rnaseq/causative.sif";
		posE.load(filename, Collections.singleton(SignedType.UPREGULATES_EXPRESSION.getTag()));
		negE.load(filename, Collections.singleton(SignedType.DOWNREGULATES_EXPRESSION.getTag()));

		GraphList gl = new GraphList("Expression TCGA consensus");
		gl.addGraph(posP);
		gl.addGraph(negP);
		gl.addGraph(posE);
		gl.addGraph(negE);
		return gl;
	}

	private static void writeAlterationsAroundCancerGenesOnBRCANetwork(String filename, Map<String, Integer> mutMap,
		Map<String, Double> cnaMap, Map<String, Double> rnaMap, Map<String, Integer> tfaMap, Set<String> cancerGenes)
		throws IOException
	{
		GraphList gl = loadTCGABRCACPNetwork();
		SiteSpecificGraph posP = (SiteSpecificGraph) gl.getGraph(SignedType.PHOSPHORYLATES.getTag());
		SiteSpecificGraph negP = (SiteSpecificGraph) gl.getGraph(SignedType.DEPHOSPHORYLATES.getTag());
		DirectedGraph posE = (DirectedGraph) gl.getGraph(SignedType.UPREGULATES_EXPRESSION.getTag());
		DirectedGraph negE = (DirectedGraph) gl.getGraph(SignedType.DOWNREGULATES_EXPRESSION.getTag());

		Set<String> gaGenes = new HashSet<>(mutMap.keySet());
		gaGenes.addAll(cnaMap.keySet());

		Set<String> cgOrTFA = new HashSet<>(cancerGenes);
		cgOrTFA.addAll(tfaMap.keySet());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));

		gaGenes.forEach(s -> posP.getDownstream(s).stream().filter(cgOrTFA::contains).forEach(t ->
			FileUtil.writeln(s + "\t" + posP.getEdgeType() + "\t" + t + "\t" + posP.getMediatorsInString(s, t), writer)));
		gaGenes.forEach(t -> posP.getUpstream(t).stream().filter(cgOrTFA::contains).forEach(s ->
			FileUtil.writeln(s + "\t" + posP.getEdgeType() + "\t" + t + "\t" + posP.getMediatorsInString(s, t), writer)));

		gaGenes.forEach(s -> negP.getDownstream(s).stream().filter(cgOrTFA::contains).forEach(t ->
			FileUtil.writeln(s + "\t" + negP.getEdgeType() + "\t" + t + "\t" + negP.getMediatorsInString(s, t), writer)));
		gaGenes.forEach(t -> negP.getUpstream(t).stream().filter(cgOrTFA::contains).forEach(s ->
			FileUtil.writeln(s + "\t" + negP.getEdgeType() + "\t" + t + "\t" + negP.getMediatorsInString(s, t), writer)));

		gaGenes.addAll(tfaMap.keySet());

		gaGenes.forEach(s -> posE.getDownstream(s).stream().filter(rnaMap::containsKey).forEach(t ->
			FileUtil.writeln(s + "\t" + posE.getEdgeType() + "\t" + t + "\t" + posE.getMediatorsInString(s, t), writer)));
		gaGenes.forEach(s -> negE.getDownstream(s).stream().filter(rnaMap::containsKey).forEach(t ->
			FileUtil.writeln(s + "\t" + negE.getEdgeType() + "\t" + t + "\t" + negE.getMediatorsInString(s, t), writer)));

		writer.close();
	}

	private static void writeTFActivityChangesOneByOne(String filename, Map<String, Integer> mutMap, Map<String, Double> cnaMap,
		Map<String, Double> rnaMap, Map<String, Integer> tfaMap, Set<String> cancerGenes) throws IOException
	{
		Set<String> altGenes = new HashSet<>(mutMap.keySet());
		altGenes.addAll(cnaMap.keySet());

		DirectedGraph scGraph = loadStateChangeNetwork();
		UndirectedGraph ppiGraph = loadPPI();

		GraphList expGraphs = loadSignedExpressionTCGAConsensus();
		DirectedGraph posE = (DirectedGraph) expGraphs.getGraph(SignedType.UPREGULATES_EXPRESSION.getTag());
		DirectedGraph negE = (DirectedGraph) expGraphs.getGraph(SignedType.DOWNREGULATES_EXPRESSION.getTag());

		for (String tf : tfaMap.keySet())
		{
			BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename + tf + ".sif"));

			scGraph.getUpstream(tf).stream()
				.filter(altGenes::contains).forEach(s ->
					FileUtil.writeln(s + "\t" + scGraph.getEdgeType() + "\t" + tf + "\t" +
						scGraph.getMediatorsInString(s, tf), writer));

			ppiGraph.getNeighbors(tf).stream()
				.filter(altGenes::contains)
				.filter(s -> !scGraph.hasRelation(s, tf))
				.forEach(s -> FileUtil.writeln(s + "\t" + ppiGraph.getEdgeType() + "\t" + tf + "\t" +
					ppiGraph.getMediatorsInString(s, tf), writer));

			posE.getDownstream(tf).stream().filter(rnaMap.keySet()::contains)
				.filter(t -> expressionEdgeSignCompatible(tfaMap.get(tf), 1, rnaMap.get(t)))
				.forEach(t -> FileUtil.writeln(tf + "\t" + posE.getEdgeType() + "\t" + t + "\t" +
					posE.getMediatorsInString(tf, t), writer));

			negE.getDownstream(tf).stream().filter(rnaMap.keySet()::contains)
				.filter(t -> expressionEdgeSignCompatible(tfaMap.get(tf), -1, rnaMap.get(t)))
				.forEach(t -> FileUtil.writeln(tf + "\t" + negE.getEdgeType() + "\t" + t + "\t" +
					negE.getMediatorsInString(tf, t), writer));

			writer.close();
		}
	}

	private static void writeTFActivityChangesInOneGraph(String filename, Map<String, Integer> mutMap, Map<String, Double> cnaMap,
		Map<String, Double> rnaMap, Map<String, Integer> tfaMap, Set<String> cancerGenes) throws IOException
	{
		Set<String> altGenes = new HashSet<>(mutMap.keySet());
		altGenes.addAll(cnaMap.keySet());

		DirectedGraph scGraph = loadStateChangeNetwork();
		UndirectedGraph ppiGraph = loadPPI();

		GraphList expGraphs = loadSignedExpressionTCGAConsensus();
		DirectedGraph posE = (DirectedGraph) expGraphs.getGraph(SignedType.UPREGULATES_EXPRESSION.getTag());
		DirectedGraph negE = (DirectedGraph) expGraphs.getGraph(SignedType.DOWNREGULATES_EXPRESSION.getTag());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));

		altGenes.forEach(s -> scGraph.getDownstream(s).stream()
			.filter(tfaMap.keySet()::contains).forEach(t ->
				FileUtil.writeln(s + "\t" + scGraph.getEdgeType() + "\t" + t + "\t" +
					scGraph.getMediatorsInString(s, t), writer)));

		altGenes.forEach(s -> ppiGraph.getNeighbors(s).stream()
			.filter(tfaMap.keySet()::contains)
			.filter(t -> !scGraph.hasRelation(s, t))
			.forEach(t -> FileUtil.writeln(s + "\t" + ppiGraph.getEdgeType() + "\t" + t + "\t" +
					ppiGraph.getMediatorsInString(s, t), writer)));

		tfaMap.keySet().forEach(s -> posE.getDownstream(s).stream().filter(rnaMap.keySet()::contains)
			.filter(t -> expressionEdgeSignCompatible(tfaMap.get(s), 1, rnaMap.get(t)))
			.forEach(t ->
				FileUtil.writeln(s + "\t" + posE.getEdgeType() + "\t" + t + "\t" + posE.getMediatorsInString(s, t), writer)));

		tfaMap.keySet().forEach(s -> negE.getDownstream(s).stream().filter(rnaMap.keySet()::contains)
			.filter(t -> expressionEdgeSignCompatible(tfaMap.get(s), -1, rnaMap.get(t)))
			.forEach(t ->
				FileUtil.writeln(s + "\t" + negE.getEdgeType() + "\t" + t + "\t" + negE.getMediatorsInString(s, t), writer)));

		writer.close();
	}

	private static boolean expressionEdgeSignCompatible(int tfActChange, int edgeSign, double rnaChange)
	{
		return tfActChange * edgeSign * rnaChange > 0;
	}

	private static void writeRelationsBetweenCancerGenes(String filename, Map<String, Integer> mutMap, Map<String, Double> cnaMap,
		Map<String, Double> rnaMap, Map<String, Integer> tfaMap, Set<String> cancerGenes) throws IOException
	{
		Set<String> usedRels = new HashSet<>();
		Set<String> usedGenes = new HashSet<>();

		Set<String> genes = new HashSet<>(mutMap.keySet());
		genes.addAll(cnaMap.keySet());
//		genes.addAll(rnaMap.keySet());
		genes.addAll(tfaMap.keySet());

		DirectedGraph scGraph = loadStateChangeNetwork();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));

		genes.stream().filter(cancerGenes::contains).forEach(s -> scGraph.getDownstream(s).stream()
			.filter(cancerGenes::contains).filter(genes::contains).forEach(t ->
			{
				String key = s + " " + t;
				usedRels.add(key);
				usedGenes.add(s);
				usedGenes.add(t);
				FileUtil.writeln(s + "\t" + scGraph.getEdgeType() + "\t" + t + "\t" +
					scGraph.getMediatorsInString(s, t), writer);
			}));

		UndirectedGraph ppiGraph = (UndirectedGraph) PathwayCommons.get().getGraph(SIFEnum.IN_COMPLEX_WITH);

		genes.stream().filter(cancerGenes::contains).forEach(s -> ppiGraph.getNeighbors(s).stream()
			.filter(cancerGenes::contains).filter(genes::contains).forEach(t ->
			{
				String key = t + " " + s;
				if (usedRels.contains(key)) return;
				key = s + " " + t;
				if (usedRels.contains(key)) return;
				usedRels.add(key);
				usedGenes.add(s);
				usedGenes.add(t);

				FileUtil.writeln(s + "\t" + ppiGraph.getEdgeType() + "\t" + t + "\t" +
					ppiGraph.getMediatorsInString(s, t), writer);
			}));

		genes.stream().filter(cancerGenes::contains).filter(g -> !usedGenes.contains(g)).forEach(g ->
			FileUtil.writeln(g, writer));

		writer.close();
	}

	private static void writeFormatFile(String filename, Map<String, Integer> mutMap, Map<String, Double> cnaMap,
		Map<String, Double> rnaMap, Map<String, Integer> tfaMap, Set<String> cancerGenes) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(filename));
		writer.write("node\tall-nodes\tcolor\t255 255 255");
		writer.write("\nnode\tall-nodes\tbordercolor\t0 0 0");

		mutMap.keySet().stream().forEach(g -> FileUtil.lnwrite("node\t" + g + "\t" + SITE_FEATURE + "\t" +
			(mutMap.get(g) < 0 ? "inhibiting " : mutMap.get(g) > 0 ? "activating " : "")  + "mutation|m|255 255 255|" +
			(mutMap.get(g) < 0 ? INHIBITING_BORDER : mutMap.get(g) > 0 ? ACTIVATING_BORDER : NEUTRAL_BORDER) + "|", writer));

		ValToColor cnaVTC = new ValToColor(new double[]{-5, 0, 5}, new Color[]{DOWN_COLOR, Color.WHITE, UP_COLOR});
		cnaMap.keySet().stream().forEach(g -> FileUtil.lnwrite("node\t" + g + "\t" + SITE_FEATURE + "\tcopy number " +
			(cnaMap.get(g) < 0 ? "loss" : "gain") + "|c|" + cnaVTC.getColorInString(cnaMap.get(g)) + "|" +
			NEUTRAL_BORDER + "|" + cnaMap.get(g), writer));

		ValToColor rnaVTC = new ValToColor(new double[]{-5, 0, 5}, new Color[]{DOWN_COLOR, Color.WHITE, UP_COLOR});
		rnaMap.keySet().stream().forEach(g -> FileUtil.lnwrite("node\t" + g + "\t" + SITE_FEATURE + "\trna expression " +
			(rnaMap.get(g) < 0 ? "decrease" : "increase") + "|e|" + rnaVTC.getColorInString(rnaMap.get(g)) + "|" +
			NEUTRAL_BORDER + "|" + rnaMap.get(g), writer));

		ValToColor tfaVTC = new ValToColor(new double[]{-3, 0, 3}, new Color[]{DOWN_COLOR, Color.WHITE, UP_COLOR});
		tfaMap.keySet().stream().forEach(g -> FileUtil.lnwrite("node\t" + g + "\t" + SITE_FEATURE + "\tTF activity " +
			(tfaMap.get(g) < 0 ? "decreased" : "increased") + "|a|" + tfaVTC.getColorInString(tfaMap.get(g)) + "|" +
			NEUTRAL_BORDER + "|", writer));

		cancerGenes.forEach(g -> FileUtil.lnwrite("node\t" + g + "\tcolor\t" + CANCER_GENE_COLOR, writer));

		writer.close();
	}

}
