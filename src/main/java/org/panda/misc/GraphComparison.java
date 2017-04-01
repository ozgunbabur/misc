package org.panda.misc;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

/**
 * Created by babur on 3/3/16.
 */
public class GraphComparison
{
	public static void main(String[] args)
	{
		DirectedGraph g1 = new DirectedGraph("PC", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		g1.load("/home/babur/Documents/mutex/networks/PC2v7.sif", new HashSet<>
			(Arrays.asList(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag(), SIFEnum.CONTROLS_EXPRESSION_OF.getTag())));
		DirectedGraph g2 = new DirectedGraph("Fries290K", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		g2.load("/home/babur/Documents/mutex/networks/fries290K.sif", new HashSet<>
			(Arrays.asList(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag(), SIFEnum.CONTROLS_EXPRESSION_OF.getTag())));
		DirectedGraph g3 = new DirectedGraph("Leidos", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		g3.load("/home/babur/Documents/mutex/networks/leidos.sif", new HashSet<>
			(Arrays.asList(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag(), SIFEnum.CONTROLS_EXPRESSION_OF.getTag())));
		DirectedGraph g4 = new DirectedGraph("fries-Ras", "is-upstream-of");
		g4.load("/home/babur/Documents/mutex/networks/fries1K.sif", new HashSet<>
			(Arrays.asList("is-upstream-of")));

		g1.printVennIntersections(g2, g3, g4);
//		g2.printVennIntersections(g4);
	}
}
