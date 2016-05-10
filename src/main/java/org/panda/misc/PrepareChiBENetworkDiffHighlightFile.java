package org.panda.misc;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.utility.graph.Graph;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by babur on 3/10/16.
 */
public class PrepareChiBENetworkDiffHighlightFile
{
	public static void main(String[] args) throws IOException
	{
		String base = "/home/babur/Documents/mutex/networks/";
		Graph g1 = new Graph("PC", "pc-edge");
		g1.load(base + "PC2v7.sif", Collections.emptySet(), new HashSet<>(Arrays.asList(
			SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag(), SIFEnum.CONTROLS_EXPRESSION_OF.getTag())));
		Graph g2 = new Graph("Fries", "fries-edge");
		g2.load(base + "fries290K.sif", Collections.emptySet(), new HashSet<>(Arrays.asList(
			SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag(), SIFEnum.CONTROLS_EXPRESSION_OF.getTag())));

		prepare(g1, g2, base + "../../../Projects/chibe/fries-diff-edge-only.highlight");
	}

	public static void prepare(Graph g1, Graph g2, String outFile) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));

		Set<String> diffGenes = new HashSet<>(g2.getSymbols());
		diffGenes.removeAll(g1.getSymbols());

//		for (String gene : diffGenes)
//		{
//			writer.write("node\t" + gene + "\n");
//		}

		for (String source : g2.getOneSideSymbols(true))
		{
			Set<String> existing = g1.getDownstream(source);

			for (String target : g2.getDownstream(source))
			{
				if (existing.contains(target)) continue;

				writer.write("edge\t" + source + "\t" + target + '\n');
			}
		}

		writer.close();
	}
}
