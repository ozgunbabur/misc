package org.panda.misc;

import org.panda.resource.network.SignedPC;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.CollectionUtil;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.PhosphoGraph;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author Ozgun Babur
 */
public class SignedSIFQuery
{
	public static void pathsBetweenFromPC(Set<String> genes, String outFile) throws IOException
	{
		Map<SignedType, DirectedGraph> graphs = SignedPC.get().getAllGraphs();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		for (SignedType type : graphs.keySet())
		{
			DirectedGraph graph = graphs.get(type);

			for (String source : genes)
			{
				Set<String> down = new HashSet<>(graph.getDownstream(source));
				down.retainAll(genes);

				for (String target : down)
				{
					writer.write(source + "\t" + type.getTag() + "\t" + target + "\t" +
						graph.getMediatorsInString(source, target));

					if (graph instanceof PhosphoGraph && ((PhosphoGraph) graph).hasSites(source, target))
					{
						writer.write("\t" + CollectionUtil.merge(((PhosphoGraph) graph).getSites(source, target), ";"));
					}

					writer.write("\n");
				}
			}
		}

		writer.close();
	}

	public static void neighborhoodFromFile(String inFile, Set<String> genes, String outFile) throws IOException
	{
		Set<String> rows = Files.lines(Paths.get(inFile)).filter(l -> l.contains("\t")).collect(Collectors.toSet());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		for (String row : rows)
		{
			String[] t = row.split("\t");
			if (t.length > 2 && (genes.contains(t[0]) || genes.contains(t[1])))
			{
				writer.write(row + "\n");
			}
		}

		writer.close();
	}

	public static void pathsBetweenFromFile(String inFile, Set<String> genes, String outFile) throws IOException
	{
		Set<String> rows = Files.lines(Paths.get(inFile)).filter(l -> l.contains("\t")).collect(Collectors.toSet());

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		for (String row : rows)
		{
			String[] t = row.split("\t");
			if (t.length > 2 && (genes.contains(t[0]) && genes.contains(t[2])))
			{
				writer.write(row + "\n");
			}
		}

		writer.close();
	}

	public static void streamOverGenes(String inFile, Set<String> genes, String outFile, int dwDist, int upDist) throws IOException
	{
		Map<String[], String> rows = Files.lines(Paths.get(inFile)).filter(l -> l.contains("\t"))
			.collect(Collectors.toMap(l -> l.split("\t"), l -> l));

		Set<String> upBorder = new HashSet<>(genes);
		Set<String> dwBorder = new HashSet<>(genes);
		Set<String> visitedUp = new HashSet<>();
		Set<String> visitedDw = new HashSet<>();

		Set<String[]> result = new HashSet<>();

		while (dwDist > 0 && !dwBorder.isEmpty())
		{
			Set<String> nextRound = new HashSet<>();
			for (String[] t : rows.keySet())
			{
				if (dwBorder.contains(t[0]) && !visitedUp.contains(t[2]) && !visitedDw.contains(t[2]) && !dwBorder.contains(t[2]))
				{
					result.add(t);
					nextRound.add(t[2]);
				}
			}
			visitedDw.addAll(dwBorder);
			nextRound.removeAll(visitedDw);
			dwBorder = nextRound;
			dwDist--;
		}

		while (upDist > 0 && !upBorder.isEmpty())
		{
			Set<String> nextRound = new HashSet<>();
			for (String[] t : rows.keySet())
			{
				if (upBorder.contains(t[2]) && !visitedUp.contains(t[0]) && !visitedDw.contains(t[0]) && !upBorder.contains(t[0]))
				{
					result.add(t);
					nextRound.add(t[0]);
				}
			}
			visitedUp.addAll(upBorder);
			nextRound.removeAll(visitedUp);
			upBorder = nextRound;
			upDist--;
		}

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		for (String[] t : result)
		{
			writer.write(rows.get(t) + "\n");
		}

		writer.close();
	}

	public static Set<String> getNodesNameStarting(String sifFile, String... prefix) throws IOException
	{
		return getNodeStreamInSIFFile(sifFile).filter(g -> startsWithOneOfThese(g, prefix)).collect(Collectors.toSet());
	}

	public static boolean startsWithOneOfThese(String s, String... q)
	{
		return Arrays.stream(q).filter(s::startsWith).findAny().isPresent();
	}

	public static Stream<String> getNodeStreamInSIFFile(String file) throws IOException
	{
		return Files.lines(Paths.get(file)).map(l -> l.split("\t")).filter(t -> t.length > 2)
			.map(t -> new String[]{t[0], t[2]}).flatMap(Arrays::stream);
	}

	public static void main(String[] args) throws IOException
	{
		// For Anil
//		String dir = "/home/babur/Documents/Analyses/JQ1/";
//		Set<String> genes = Files.lines(Paths.get(dir + "platform.txt")).skip(1)
//			.map(l -> l.split("\t")[2].split(" ")).flatMap(Arrays::stream).collect(Collectors.toSet());
//		pathsBetweenFromPC(genes, dir + "base-graph.sif");

//		String dir = "/home/babur/Documents/Analyses/CPTACBreastCancer/correlation-based/";
//		Set<String> genes = Files.lines(Paths.get(dir + "significance-pvals.txt")).skip(2).map(l -> l.split("\t"))
//			.filter(t -> t.length > 1 && Double.valueOf(t[1]) < 0.01).map(t -> t[0]).collect(Collectors.toSet());
//		neighborhoodFromFile(dir + "causative.sif", genes, dir + "causative-subset.sif");

//		String dir = "/home/babur/Documents/Analyses/CPTACBreastCancer/correlation-based/";
//		streamOverGenes(dir + "causative.sif", Collections.singleton("SHC1"), dir + "causative-subset.sif", 0, 2);

		String dir = "/home/babur/Documents/Analyses/CPTACBreastCancer/correlation-based/";
		String inFile = dir + "causative-data-centric.sif";
		pathsBetweenFromFile(inFile, getNodesNameStarting(inFile, "BRAF", "AKT", "PPP"), dir + "causative-data-centric-subset.sif");

	}
}
