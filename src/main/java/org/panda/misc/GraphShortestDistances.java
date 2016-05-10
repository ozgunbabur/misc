package org.panda.misc;

import org.panda.resource.network.SignedPC;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.graph.Graph;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

/**
 * Calculates and writes shortest distances in a given graph.
 *
 * @author Ozgun Babur
 */
public class GraphShortestDistances
{
	public static void main(String[] args) throws IOException
	{
		SignedPC spc = new SignedPC();
		Graph graph = spc.getGraph(SignedType.values());
		graph.printStats();

		writeShortestDistances(graph, 4, "/home/babur/Documents/Temp/distances.txt");
	}

	public static void writeShortestDistances(Graph graph, int limit, String filename) throws IOException
	{
		Map<String, Map<String, Integer>> dist = graph.getShortestDistances(limit);
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		for (String source : dist.keySet())
		{
			for (String target : dist.get(source).keySet())
			{
				writer.write(source + "\t" + target + "\t" + dist.get(source).get(target) + "\n");
			}
		}

		writer.close();
	}

}
