package org.panda.misc;

import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.UndirectedGraphWithEdgeWeights;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Set;

/**
 * Created by Owner on 12/23/2019.
 */
public class SIFDownstreamSimilarity
{
	public static void main(String[] args)
	{
		processRecursively("/Users/ozgun/Documents/Analyses/CPTAC-LSCC/clusters");
	}

	public static void processRecursively(String dir)
	{
		try
		{
			if (Files.exists(Paths.get(dir + "/causative-sig-neigh.sif")))
			{
				generateDownstreamSimilarityGraph(dir + "/causative-sig-neigh.sif", dir + "/causative-sig-neigh-similarity.sif");
			}

			Files.list(Paths.get(dir)).filter(Files::isDirectory).forEach(p -> processRecursively(p.toString()));
		} catch (IOException e) {e.printStackTrace();}
	}

	public static void generateDownstreamSimilarityGraph(String sifFile, String outFile) throws IOException
	{
		DirectedGraph sif = new DirectedGraph("SIF", "directed");
		sif.load(sifFile, FileUtil.getTermsInTabDelimitedColumn(sifFile, 1, 0));

		UndirectedGraphWithEdgeWeights simGraph = sif.getDownstreamSimilarityGraph();
		simGraph.setEdgeType("interacts-with");

		BufferedWriter writer1 = Files.newBufferedWriter(Paths.get(outFile));
		simGraph.write(writer1);
		Set<String> syms = sif.getOneSideSymbols(true);
		syms.removeAll(simGraph.getSymbols());
		syms.forEach(s -> FileUtil.writeln(s, writer1));
		writer1.close();

		outFile = outFile.substring(0, outFile.lastIndexOf(".")) + ".format";
		BufferedWriter writer2 = Files.newBufferedWriter(Paths.get(outFile));
		writer2.write("edge\tall-edges\twidth\t2\n");

		ValToColor vtc = new ValToColor(new double[]{0, 1}, new Color[]{new Color(200, 200, 200), Color.BLACK});

		for (String s1 : simGraph.getSymbols())
		{
			for (String s2 : simGraph.getNeighbors(s1))
			{
				if (s1.compareTo(s2) < 0)
				{
					writer2.write("edge\t"+s1 + " " + simGraph.getEdgeType() + " " + s2 + "\tcolor\t" +
						vtc.getColorInString(simGraph.getWeight(s1, s2)) + "\n");
				}
			}
		}

		writer2.close();
	}
}
