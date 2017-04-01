package org.panda.misc.pancan;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.misc.MutexReader;
import org.panda.misc.MutexReader.*;
import org.panda.resource.network.PathwayCommons;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class MutexToGraphAligner
{
	static final Map<String, ValToColor> edgeColorsMap = new HashMap<>();
	final static double[] range = {-1, -Math.log(0.001)};
	static
	{
		edgeColorsMap.put(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag(),
			new ValToColor(range, new Color[]{Color.WHITE, new Color(50, 100, 150).darker()}));
		edgeColorsMap.put(SIFEnum.CONTROLS_EXPRESSION_OF.getTag(),
			new ValToColor(range, new Color[]{Color.WHITE, new Color(50, 150, 50).darker()}));
		edgeColorsMap.put(SIFEnum.IN_COMPLEX_WITH.getTag(),
			new ValToColor(range, new Color[]{Color.WHITE, new Color(100, 100, 100)}));
		edgeColorsMap.put("in-same-group",
			new ValToColor(range, new Color[]{Color.WHITE, new Color(50, 50, 50)}));
	}
	public static final ValToColor pvalVTC = new ValToColor(range, new Color[]{Color.WHITE, Color.BLACK});
	public static final ValToColor covVTC = new ValToColor(new double[]{0, 0.2}, new Color[]{Color.WHITE, Color.RED});

	static void writeNetworkSupportingMutexResults(Set<Group> results, double scoreThr, List<DirectedGraph> graphs,
		String sifFileNameWithoutExtension) throws IOException
	{
		Map<Pair, Double> pairMap = preparePairMap(results);
		Map<String, Double> geneScoreMap = prepareGeneScoreMap(results);
		Map<String, Double> coverage = readMutationCoverages(geneScoreMap.keySet());
		Set<String> nodesInGraph = new HashSet<>();
		Set<Pair> edgesInGraph = new HashSet<>();

		BufferedWriter sifWriter = Files.newBufferedWriter(Paths.get(sifFileNameWithoutExtension + ".sif"));
		BufferedWriter formatWriter = Files.newBufferedWriter(Paths.get(sifFileNameWithoutExtension + ".format"));
		formatWriter.write("node\tall-nodes\tcolor\t255 255 255");
		formatWriter.write("\nnode\tall-nodes\tborderwidth\t2");
		formatWriter.write("\nedge\tall-edges\twidth\t2");

		pairMap.keySet().stream().filter(pair -> pairMap.get(pair) <= scoreThr).forEach(pair ->
		{
			for (DirectedGraph graph : graphs)
			{
				if (graph.isDirected())
				{
					if (graph.getDownstream(pair.gene1).contains(pair.gene2))
					{
						addEdge(pair.gene1, pair.gene2, pairMap.get(pair), graph.getEdgeType(), nodesInGraph, edgesInGraph, pair, sifWriter, formatWriter);
					}
					if (graph.getDownstream(pair.gene2).contains(pair.gene1))
					{
						addEdge(pair.gene2, pair.gene1, pairMap.get(pair), graph.getEdgeType(), nodesInGraph, edgesInGraph, pair, sifWriter, formatWriter);
					}
				}
				else
				{
					if (graph.getNeighbors(pair.gene1).contains(pair.gene2))
					{
						addEdge(pair.gene1, pair.gene2, pairMap.get(pair), graph.getEdgeType(), nodesInGraph, edgesInGraph, pair, sifWriter, formatWriter);
					}
				}
			}
		});

		NumberFormat fmt = new DecimalFormat("0.####");

		nodesInGraph.forEach(gene -> FileUtil.lnwrite(
			"node\t" + gene + "\tbordercolor\t" + pvalVTC.getColorInString(-Math.log(geneScoreMap.get(gene))) +
			"\nnode\t" + gene + "\tcolor\t" + covVTC.getColorInString(coverage.get(gene)) +
				"\nnode\t" + gene + "\ttooltip\tp=" + fmt.format(geneScoreMap.get(gene)) + ", c=" + fmt.format(coverage.get(gene)), formatWriter));

		sifWriter.close();
		formatWriter.close();

		Set<String> genes = geneScoreMap.keySet().stream().filter(g -> geneScoreMap.get(g) <= scoreThr).collect(Collectors.toSet());
		genes.removeAll(nodesInGraph);
		System.out.println("nodesInGraph = " + nodesInGraph.size());
		System.out.println("genes.size() = " + genes.size());
		System.out.println("genes = " + genes);
	}

	static void writeInSameGroupNetwork(Set<Group> results, double scoreThr, String sifFileNameWithoutExtension) throws IOException
	{
		Map<Pair, Double> pairMap = preparePairMap(results);
		Map<String, Double> geneScoreMap = prepareGeneScoreMap(results);
		Map<String, Double> coverage = readMutationCoverages(geneScoreMap.keySet());
		Set<String> nodesInGraph = new HashSet<>();
		Set<Pair> edgesInGraph = new HashSet<>();

		BufferedWriter sifWriter = Files.newBufferedWriter(Paths.get(sifFileNameWithoutExtension + ".sif"));
		BufferedWriter formatWriter = Files.newBufferedWriter(Paths.get(sifFileNameWithoutExtension + ".format"));
		formatWriter.write("node\tall-nodes\tcolor\t255 255 255");
		formatWriter.write("\nnode\tall-nodes\tborderwidth\t2");
//		formatWriter.write("\nedge\tall-edges\twidth\t2");

		pairMap.keySet().stream().filter(pair -> pairMap.get(pair) <= scoreThr).forEach(pair ->
			addEdge(pair.gene1, pair.gene2, pairMap.get(pair), "in-same-group", nodesInGraph, edgesInGraph, pair, sifWriter, formatWriter));

		NumberFormat fmt = new DecimalFormat("0.####");

		nodesInGraph.forEach(gene -> FileUtil.lnwrite(
			"node\t" + gene + "\tbordercolor\t" + pvalVTC.getColorInString(-Math.log(geneScoreMap.get(gene))) +
			"\nnode\t" + gene + "\tcolor\t" + covVTC.getColorInString(coverage.get(gene)) +
				"\nnode\t" + gene + "\ttooltip\tp=" + fmt.format(geneScoreMap.get(gene)) + ", c=" + fmt.format(coverage.get(gene)), formatWriter));

		sifWriter.close();
		formatWriter.close();
	}

	static void addEdge(String g1, String g2, double score, String edgeType,
		Set<String> nodesSoFar, Set<Pair> edgesSoFar, Pair pair, BufferedWriter sifWriter, BufferedWriter formatWriter)
	{
		if (edgesSoFar.contains(pair)) return;
		FileUtil.writeln(g1 + "\t" + edgeType + "\t" + g2, sifWriter);
		FileUtil.lnwrite("edge\t" + g1 + " " + edgeType + " " + g2 + "\tcolor\t" +
			edgeColorsMap.get(edgeType).getColorInString(-Math.log(score)), formatWriter);
		nodesSoFar.add(g1);
		nodesSoFar.add(g2);
		edgesSoFar.add(pair);
	}

	static List<DirectedGraph> loadGraphs()
	{
		List<DirectedGraph> graphs = new ArrayList<>();
		DirectedGraph csc = (DirectedGraph) PathwayCommons.get().getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF);
//		Graph g = new Graph("Reach", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
//		g.load("/home/babur/Documents/mutex/networks/REACH.sif", Collections.emptySet(), Collections.singleton(g.getEdgeType()));
//		csc.merge(g);
		graphs.add(csc);
		graphs.add((DirectedGraph) PathwayCommons.get().getGraph(SIFEnum.CONTROLS_EXPRESSION_OF));

//		Graph ppi = PathwayCommons.get().getGraph(SIFEnum.IN_COMPLEX_WITH);
//		ppi.merge(PathwayCommons.get().getGraph(SIFEnum.INTERACTS_WITH));
//		graphs.add(ppi);

		return graphs;
	}

	static Map<String, Double> prepareGeneScoreMap(Set<Group> results)
	{
		Map<String, Double> scores = new HashMap<>();

		for (Group group : results)
		{
			for (String gene : group.genes)
			{
				if (!scores.containsKey(gene) || scores.get(gene) > group.score)
				{
					scores.put(gene, group.score);
				}
			}
		}
		return scores;
	}
	static Map<Pair, Double> preparePairMap(Set<Group> results)
	{
		Map<Pair, Double> pairs = new HashMap<>();

		for (Group group : results)
		{
			for (String g1 : group.genes)
			{
				for (String g2 : group.genes)
				{
					if (g1.compareTo(g2) < 0)
					{
						Pair p = new Pair(g1, g2);
						if (!pairs.containsKey(p) || pairs.get(p) > group.score) pairs.put(p, group.score);
					}
				}
			}
		}
		return pairs;
	}

	public static Map<String, Double> readMutationCoverages(Set<String> genes) throws IOException
	{
		Map<String, Double> covMap = new HashMap<>();
		Files.lines(Paths.get("/home/babur/Documents/mutex/TCGA/PanCan/1/1/DataMatrix.txt")).skip(1)
			.filter(l -> genes.contains(l.substring(0, l.indexOf("\t"))))
			.map(l -> l.split("\t")).forEach(t ->
		{
			int c = 0;
			for (int i = 1; i < t.length; i++)
			{
				if (!t[i].equals("0")) c++;
			}
			covMap.put(t[0], c / (double) (t.length - 1));
		});
		return covMap;
	}

	static class Pair
	{
		String gene1;
		String gene2;

		public Pair(String gene1, String gene2)
		{
			if (gene1.compareTo(gene2) < 0)
			{
				this.gene1 = gene1;
				this.gene2 = gene2;
			}
			else if (gene1.equals(gene2))
			{
				throw new IllegalArgumentException("Same genes!");
			}
			else
			{
				this.gene2 = gene1;
				this.gene1 = gene2;
			}
		}

		@Override
		public boolean equals(Object obj)
		{
			return obj instanceof Pair && gene1.equals(((Pair) obj).gene1) && gene2.equals(((Pair) obj).gene2);
		}

		@Override
		public int hashCode()
		{
			return gene1.hashCode() + gene2.hashCode();
		}
	}

	public static void main(String[] args) throws IOException
	{
		String base = "/home/babur/Documents/PanCan/";
		Object[] o = PanCanResultLoader.readGroupsWithFlattenedControl(true, base + "PanCan-results",
			base + "PanCan-shuffled-?-results", f -> useDir(f.getName()));
		Set<Group> groups = (Set<Group>) o[0];

//		MutexReader.readMutexResultsRecursive(base + "PanCan-results", new HashSet<>(),
//			f -> useDir(f.getName()));

		writeNetworkSupportingMutexResults(groups, 0.1475, loadGraphs(), base + "aligned-network-PC");
//		writeInSameGroupNetwork(groups, 0.01, base + "mutex-pairings");
	}

	static boolean useDir(String dirName)
	{
		return Character.isDigit(dirName.charAt(0)) || dirName.equals("no-network") || dirName.equals("REACH-PC2v8") || dirName.equals("fries290K-PC2v8");
	}
}
