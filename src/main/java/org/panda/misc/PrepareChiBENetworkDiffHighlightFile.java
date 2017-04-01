package org.panda.misc;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.utility.graph.DirectedGraph;
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
		drawVenn();
		if (true) return;

		String base = "/media/babur/6TB1/REACH-mutex/network/";
		DirectedGraph g1 = new DirectedGraph("PC", "pc-edge");
		g1.load(base + "PC2v8.sif", new HashSet<>(Arrays.asList(
			SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag(), SIFEnum.CONTROLS_EXPRESSION_OF.getTag())));
		DirectedGraph g2 = new DirectedGraph("Fries", "fries-edge");
		g2.load(base + "REACH.sif", new HashSet<>(Arrays.asList(
			SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag(), SIFEnum.CONTROLS_EXPRESSION_OF.getTag())));

		prepare(g1, g2, "/home/babur/Projects/chibe/fries-diff.highlight");
	}

	public static void drawVenn() throws IOException
	{
		String base = "/media/babur/6TB1/REACH-mutex/network/";
		DirectedGraph g1 = new DirectedGraph("PC", "pc-edge");
		g1.load(base + "PC2v8.sif", Collections.singleton(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag()));
		DirectedGraph g2 = new DirectedGraph("Fries", "fries-edge");
		g2.load(base + "REACH.sif", Collections.singleton(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag()));

		g1.printVennIntersections(g2);
	}

	/**
	 * Prepares a .highlight file to use in ChiBE sif views. This file highlights the components that exists in the
	 * second graph but not the first graph.
	 */
	public static void prepare(DirectedGraph g1, DirectedGraph g2, String outFile) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));

		Set<String> diffGenes = new HashSet<>(g2.getSymbols());
		diffGenes.removeAll(g1.getSymbols());

		for (String gene : diffGenes)
		{
			writer.write("node\t" + gene + "\n");
		}

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
