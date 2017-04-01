package org.panda.misc;

import org.panda.resource.network.SignedPC;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.CollectionUtil;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.PhosphoGraph;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class SignedSIFPathsBetween
{
	public static void query(Set<String> genes, String outFile) throws IOException
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

	public static void main(String[] args) throws IOException
	{
		String dir = "/home/babur/Documents/Analyses/JQ1/";
		Set<String> genes = Files.lines(Paths.get(dir + "platform.txt")).skip(1)
			.map(l -> l.split("\t")[2].split(" ")).flatMap(Arrays::stream).collect(Collectors.toSet());

		query(genes, dir + "base-graph.sif");
	}
}
